import argparse
import subprocess
import os
import sys

def run_command(cmd, description, log_file):
    """Runs a shell command and logs the output."""
    print(f"\nRunning: {description}")
    with open(log_file, "a") as log:
        log.write(f"\nRunning: {description}\n")
        log.write(f"Command: {cmd}\n")
        
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        log.write(stdout + "\n")
        if stderr:
            log.write(f"ERROR in {description}: {stderr}\n")
        
        if process.returncode != 0:
            log.write("Exiting pipeline due to error.\n")
            print(f"ERROR in {description}: {stderr}")
            sys.exit(1)
        else:
            log.write(f"{description} completed successfully.\n")
            print(f"{description} completed successfully.")

def check_and_build_bowtie2_index(genome_fasta, genome_index_prefix, log_file):
    """Checks if the Bowtie2 index exists; if not, builds it."""
    expected_files = [f"{genome_index_prefix}.{i}.bt2" for i in range(1, 5)]
    
    if all(os.path.exists(f) for f in expected_files):
        print("✅ Bowtie2 index found. Skipping indexing step.")
        return
    else:
        print("⚠️ Bowtie2 index not found. Building index now...")
        build_cmd = f"bowtie2-build {genome_fasta} {genome_index_prefix}"
        run_command(build_cmd, "Building Bowtie2 index", log_file)

def main():
    parser = argparse.ArgumentParser(description="Optimized Bulk RNA-Seq Pipeline for Multiple Samples")
    parser.add_argument("--input_dir", required=True, help="Directory containing raw FASTQ files")
    parser.add_argument("--output_dir", required=True, help="Directory to store results")
    parser.add_argument("--genome_fasta", required=True, help="Path to genome FASTA file")
    parser.add_argument("--genome_index", required=True, help="Prefix for Bowtie2 genome index")
    parser.add_argument("--annotation", required=True, help="Path to GTF annotation file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    log_file = os.path.join(args.output_dir, "pipeline.log")
    
    with open(log_file, "w") as log:
        log.write("RNA-Seq Pipeline Log\n")
    
    # ✅ Step 0: Check and Build Bowtie2 Index if Needed
    check_and_build_bowtie2_index(args.genome_fasta, args.genome_index, log_file)
    
    fastq_files = sorted([f for f in os.listdir(args.input_dir) if f.endswith(".fastq.gz")])
    sample_pairs = {}
    
    for fastq in fastq_files:
        sample_name = fastq.replace("_R1.fastq.gz", "").replace("_R2.fastq.gz", "")
        if sample_name not in sample_pairs:
            sample_pairs[sample_name] = {"R1": None, "R2": None}
        if "_R1.fastq.gz" in fastq:
            sample_pairs[sample_name]["R1"] = os.path.join(args.input_dir, fastq)
        elif "_R2.fastq.gz" in fastq:
            sample_pairs[sample_name]["R2"] = os.path.join(args.input_dir, fastq)
    
    for sample, files in sample_pairs.items():
        if not files["R1"] or not files["R2"]:
            print(f"Skipping {sample} as it lacks a complete R1/R2 pair.")
            continue
        
        sample_output = os.path.join(args.output_dir, sample)
        os.makedirs(sample_output, exist_ok=True)
        
        # Step 1: Run FastQC
        raw_fastqc_dir = os.path.join(sample_output, "raw_fastqc")
        os.makedirs(raw_fastqc_dir, exist_ok=True)
        run_command(f"fastqc {files['R1']} {files['R2']} -o {raw_fastqc_dir}", f"FastQC on raw reads for {sample}", log_file)
        
        # Step 2: Trim Reads
        trimmed_fastq1 = os.path.join(sample_output, "trimmed_1.fastq.gz")
        trimmed_fastq2 = os.path.join(sample_output, "trimmed_2.fastq.gz")
        trimmomatic_cmd = (
            f"trimmomatic PE -threads {args.threads} {files['R1']} {files['R2']} "
            f"{trimmed_fastq1} /dev/null {trimmed_fastq2} /dev/null "
            "SLIDINGWINDOW:4:20 MINLEN:36"
        )
        run_command(trimmomatic_cmd, f"Trimming low-quality reads for {sample}", log_file)
        
        # Step 3: Run FastQC Again
        trimmed_fastqc_dir = os.path.join(sample_output, "trimmed_fastqc")
        os.makedirs(trimmed_fastqc_dir, exist_ok=True)
        run_command(f"fastqc {trimmed_fastq1} {trimmed_fastq2} -o {trimmed_fastqc_dir}", f"FastQC after trimming for {sample}", log_file)
        
        # Step 4: Align Reads using Bowtie2
        aligned_sam = os.path.join(sample_output, "aligned.sam")
        bowtie2_cmd = (
            f"bowtie2 -x {args.genome_index} -1 {trimmed_fastq1} -2 {trimmed_fastq2} "
            f"-S {aligned_sam} --threads {args.threads}")
        run_command(bowtie2_cmd, f"Aligning reads for {sample} using Bowtie2", log_file)

        # Convert SAM to BAM and Sort
        aligned_bam = os.path.join(sample_output, "aligned.bam")
        run_command(f"samtools view -bS {aligned_sam} | samtools sort -o {aligned_bam}", 
                    f"Converting and sorting BAM file for {sample}", log_file)
        
        # Step 5: Generate Count Matrix
        count_matrix = os.path.join(sample_output, "counts.txt")
        featurecounts_cmd = (
            f"featureCounts -T {args.threads} -p -a {args.annotation} -o {count_matrix} {aligned_bam}"
            )

        run_command(featurecounts_cmd, f"Generating count matrix for {sample}", log_file)
        
    print("\nPipeline completed successfully for all samples! 🎉")
    with open(log_file, "a") as log:
        log.write("\nPipeline completed successfully for all samples!\n")
    
if __name__ == "__main__":
    main()
