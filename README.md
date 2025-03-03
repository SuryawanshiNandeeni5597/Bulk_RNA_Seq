# Bulk RNA-Seq Analysis Pipeline

## Overview
This is the Bulk RNA-Seq pipeline designed to process multiple RNA-Seq samples efficiently. The pipeline automates the following steps:

1. **FastQC**: Quality control of raw FASTQ files.
2. **Trimmomatic**: Trimming low-quality reads.
3. **FastQC (Post-trimming)**: Quality control check after trimming.
4. **Bowtie2 Alignment**: Aligning reads to a reference genome.
5. **SAM to BAM Conversion**: Converting aligned SAM files to sorted BAM format.
6. **FeatureCounts**: Generating a gene expression count matrix.

## Requirements
Ensure the following tools are installed:
- Python 3.x
- FastQC
- Trimmomatic
- Bowtie2
- SAMtools
- FeatureCounts

## Installation
1. Clone this repository:
   ```bash
   git clone <repository_link>
   cd bulk_RNA_seq_pipeline
   ```
2. Install dependencies:
   ```bash
   conda create -n rnaseq_pipeline python=3.8 -y
   conda activate rnaseq_pipeline
   conda install -c bioconda fastqc trimmomatic bowtie2 samtools subread
   ```

## Usage
Run the pipeline using the following command:
```bash
python bulkRNA_seq_Analysis_pipeline.py \
    --input_dir /path/to/input \
    --output_dir /path/to/output \
    --genome_fasta /path/to/genome.fasta \
    --genome_index /path/to/genome_index \
    --annotation /path/to/annotation.gtf \
    --threads 4
```

### Parameters
- `--input_dir` : Directory containing raw FASTQ files.
- `--output_dir` : Directory to store results.
- `--genome_fasta` : Path to the genome FASTA file.
- `--genome_index` : Prefix for the Bowtie2 genome index.
- `--annotation` : Path to the GTF annotation file.
- `--threads` : Number of CPU threads to use (default: 4).

## Output
- `raw_fastqc/` : FastQC reports for raw reads.
- `trimmed_fastqc/` : FastQC reports after trimming.
- `trimmed_1.fastq.gz`, `trimmed_2.fastq.gz` : Trimmed paired-end reads.
- `aligned.sam` : Aligned reads in SAM format.
- `aligned.bam` : Sorted BAM file.
- `counts.txt` : Gene expression count matrix.
- `pipeline.log` : Log file recording pipeline progress.

## Author
**Nandeeni Suryawanshi**

 
