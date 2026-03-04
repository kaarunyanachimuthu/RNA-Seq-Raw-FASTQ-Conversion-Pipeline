# RNA-Seq-Raw-FASTQ-Conversion-Pipeline

A bash pipeline designed sequentially for processing paired-end bulk RNA-seq data on HPC cluster (LSF + Docker). Takes raw FASTQ files through quality control, alignment, and quantification to produce a gene counts matrix ready for downstream analysis (e.g. DESeq2, edgeR).

---

## Table of Contents

- [Overview](#overview)
- [Samples](#samples)
- [Pipeline Steps](#pipeline-steps)
- [Requirements](#requirements)
- [Directory Structure](#directory-structure)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output Files](#output-files)
- [Tools & Versions](#tools--versions)

---

## Overview

```
Raw FASTQ (4 lanes)
      │
      ▼
  [Step 1] fastp — Quality trimming + adapter removal (per lane → merged)
      │
      ▼
  [Step 2] bowtie2 — Paired-end alignment to hg38
      │
      ▼
  [Step 3] samtools — SAM → sorted BAM → index
      │
      ▼
  [Step 4] featureCounts — Gene-level read quantification
      │
      ▼
  all_samples_counts.txt  ← ready for DESeq2 / edgeR
```

---

## Samples

| Sample ID   | Description   |
|-------------|---------------|
| Proband_1   | Proband replicate 1 |
| Proband_2   | Proband replicate 2 |
| Dad_1       | Father replicate 1  |
| Dad_2       | Father replicate 2  |
| Mum_1       | Mother replicate 1  |
| Mum_2       | Mother replicate 2  |

Each sample has **4 sequencing lanes** (L001–L004) that are trimmed individually and then merged before alignment.

---

## Pipeline Steps

### Step 1 — fastp (QC & Trimming)

Trims adapters and low-quality bases from paired-end reads for each sample across all 4 lanes, then concatenates lanes into a single merged FASTQ per sample.

**Parameters used:**
- `--detect_adapter_for_pe` — auto-detect adapters for paired-end reads
- `--qualified_quality_phred 20` — filter bases below Q20
- `--length_required 35` — discard reads shorter than 35 bp
- `--thread 4`

**Docker image:** `staphb/fastp`

---

### Step 2 — bowtie2 (Alignment)

Aligns merged paired-end reads to the **hg38** human reference genome using bowtie2.

**Parameters used:**
- `-p 8` — 8 threads
- `-x` — path to bowtie2 hg38 index
- `-1` / `-2` — R1 and R2 merged clean FASTQs
- `-S` — output SAM file

**Docker image:** `quay.io/biocontainers/subread:2.0.1--hed695b0_0`

> **Note:** The hg38 bowtie2 index must be pre-built. See [Configuration](#configuration).

---

### Step 3 — samtools (Sort & Index)

Converts each SAM to a coordinate-sorted BAM, indexes it, and removes the intermediate SAM to save disk space.

```bash
samtools view -b sample.sam | samtools sort -T sample.tmp -o sample.sorted.bam
samtools index sample.sorted.bam
```

**Docker image:** `quay.io/biocontainers/samtools:1.2-0`

---

### Step 4 — featureCounts (Quantification)

Counts reads mapping to each gene across all 6 samples simultaneously, using a hg38 GTF annotation. Outputs a single counts matrix.

**Parameters used:**
- `-p --countReadPairs` — paired-end mode
- `-s 0` — unstranded (adjust to `1` or `2` if your library is stranded)
- `-T 8` — 8 threads
- `-a` — GTF annotation file

**Docker image:** `quay.io/biocontainers/subread:2.0.1--hed695b0_0`

---

## Requirements

- WashU RIS HPC cluster access (LSF job scheduler)
- Docker-enabled compute nodes
- The following Docker images (pulled automatically on first run):
  - `staphb/fastp`
  - `quay.io/biocontainers/samtools:1.2-0`
  - `quay.io/biocontainers/subread:2.0.1--hed695b0_0`
- Pre-built **hg38 bowtie2 index** (see below)
- **hg38 GTF annotation file** (see below)

### Building the hg38 bowtie2 index

```bash
# Download hg38 reference genome
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
gunzip hg38.fa.gz

# Build index (run inside bowtie2 container)
bowtie2-build hg38.fa hg38_index
```

### Downloading the hg38 GTF annotation

```bash
# From Gencode (recommended)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
gunzip gencode.v45.annotation.gtf.gz
```

---

## Directory Structure

```
kaarunya_files/
├── FASTQ_uploads/                  # Raw input FASTQs (one per sample/lane/read)
│   ├── 180622_18_07769_A_S1_L001_R1_001.fastq.gz
│   └── ...
├── fastp_output/                   # Trimmed FASTQs + QC reports
│   ├── Proband_1/                  # Per-lane clean FASTQs + HTML/JSON reports
│   ├── Proband_1_R1.fastq.gz       # Merged clean R1 (used for alignment)
│   ├── Proband_1_R2.fastq.gz       # Merged clean R2
│   └── ...
├── alignments_bt2/                 # Alignment outputs
│   ├── Proband_1.sorted.bam
│   ├── Proband_1.sorted.bam.bai
│   └── ...
├── featurecounts_output/           # Final counts
│   └── all_samples_counts.txt
├── pipeline_logs/                  # Log files
└── hg38_index.*                    # bowtie2 index files
```

---

## Configuration

Before running, open `bulk_rnaseq_pipeline.sh` and update the following variables from the beginning of the script:

```bash
# Your HPC paths
export SCRATCH1=/scratch1/fs1/bigley
export HOME=/home/kaarunya
export STORAGE1=/storage1/fs1/bigley/Active

# Path to pre-built bowtie2 hg38 index prefix
INDEX="${BASE_DIR}/hg38_index"

# ⚠️  UPDATE THIS: Path to your hg38 GTF annotation file
GTF="/path/to/hg38.gtf"
```

Also check the `SAMPLE_PREFIX` mapping if your raw FASTQ filenames differ from those listed in the script.

---

## Usage

```bash
# Clone the repository
git clone https://github.com/your-username/your-repo.git
cd your-repo

# Make the script executable
chmod +x bulk_rnaseq_pipeline.sh

# Run the pipeline
bash bulk_rnaseq_pipeline.sh
```

Each step is submitted as an interactive LSF job. The pipeline will not advance to the next step until the previous one completes successfully — if any step fails, the pipeline exits immediately (`set -euo pipefail`).

---

## Output Files

| File | Description |
|------|-------------|
| `fastp_output/<sample>/<sample>_<lane>_fastp.html` | Per-lane QC report (open in browser) |
| `fastp_output/<sample>/<sample>_<lane>_fastp.json` | Per-lane QC metrics (machine-readable) |
| `fastp_output/<sample>_R1/R2.fastq.gz` | Merged trimmed reads per sample |
| `alignments_bt2/<sample>.sorted.bam` | Coordinate-sorted alignment |
| `alignments_bt2/<sample>.sorted.bam.bai` | BAM index |
| `alignments_bt2/<sample>.bowtie2.log` | Alignment summary statistics |
| `featurecounts_output/all_samples_counts.txt` | **Final gene counts matrix** |

The counts matrix can be loaded directly into R for differential expression analysis:

```r
counts <- read.table("all_samples_counts.txt", header=TRUE, skip=1, row.names=1)
```

---

## Tools & Versions

| Tool | Version | Purpose |
|------|---------|---------|
| [fastp](https://github.com/OpenGene/fastp) | 1.0.1 | QC & adapter trimming |
| [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | 2.5.4 | Read alignment |
| [samtools](http://www.htslib.org/) | 1.2 | BAM processing |
| [featureCounts](https://subread.sourceforge.net/) | 2.0.1 | Read quantification |

---

## Notes

- The intermediate `.sam` files are automatically deleted after sorting to conserve disk space.
- Each `bsub` job requests 90 GB RAM. Adjust `MEM` in the script if needed.
