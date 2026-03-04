#!/bin/bash
# =============================================================================
# Bulk RNA-seq Pipeline: FASTQ → Counts
# Samples: Proband_1, Proband_2, Dad_1, Dad_2, Mum_1, Mum_2
# =============================================================================

# =============================================================================
# CONFIGURE PATHS — update these before running
# =============================================================================

export SCRATCH1=/path/to/scratch/directory
export HOME=/home/kaarunya
export STORAGE1=/path/to/storage/directory
export LSF_DOCKER_VOLUMES="$HOME:$HOME $STORAGE1:$STORAGE1 $SCRATCH1:$SCRATCH1"

BASE_DIR="${STORAGE1}/scRNAseq_DataSets/kaarunya_files"
RAW_DIR="${BASE_DIR}/FASTQ_uploads"
FASTP_DIR="${BASE_DIR}/fastp_output"
ALIGN_DIR="${BASE_DIR}/alignments_bt2"
COUNTS_DIR="${BASE_DIR}/featurecounts_output"
LOG_DIR="${BASE_DIR}/pipeline_logs"

INDEX="${BASE_DIR}/hg38_index"           
GTF="/path/to/hg38.gtf"                  

SAMPLES=(Proband_1 Proband_2 Dad_1 Dad_2 Mum_1 Mum_2)
LANES=(L001 L002 L003 L004)

declare -A SAMPLE_PREFIX=(
    [Proband_1]="180622_18_07769_A_S1"
    [Proband_2]="180622_18_07769_B_S2"
    [Dad_1]="180724_18_11842_b_1_S1"
    [Dad_2]="180724_18_11842_b_2_S2"
    [Mum_1]="180724_18_11844_b_1_S3"
    [Mum_2]="180724_18_11844_b_2_S4"
)

THREADS=8
MEM=90000  

# =============================================================================
# SETUP — create output directories
# =============================================================================

echo "=============================================="
echo " Bulk RNA-seq Pipeline Starting"
echo " $(date)"
echo "=============================================="

mkdir -p "${FASTP_DIR}"/{Proband_1,Proband_2,Dad_1,Dad_2,Mum_1,Mum_2}
mkdir -p "${ALIGN_DIR}"
mkdir -p "${COUNTS_DIR}"
mkdir -p "${LOG_DIR}"

echo "[INFO] Output directories ready."

# =============================================================================
# STEP 1: fastp — Quality Control & Adapter Trimming (per lane, then merge)
# =============================================================================

echo ""
echo "----------------------------------------------"
echo " STEP 1: fastp — Trimming & QC"
echo "----------------------------------------------"

bsub -Is -q general-interactive -M ${MEM} -R "rusage[mem=${MEM}]" \
     -a 'docker(staphb/fastp)' /bin/bash << 'FASTP_BLOCK'

set -euo pipefail

RAW_DIR="${STORAGE1}/scRNAseq_DataSets/kaarunya_files/FASTQ_uploads"
FASTP_DIR="${STORAGE1}/scRNAseq_DataSets/kaarunya_files/fastp_output"

declare -A SAMPLE_PREFIX=(
    [Proband_1]="180622_18_07769_A_S1"
    [Proband_2]="180622_18_07769_B_S2"
    [Dad_1]="180724_18_11842_b_1_S1"
    [Dad_2]="180724_18_11842_b_2_S2"
    [Mum_1]="180724_18_11844_b_1_S3"
    [Mum_2]="180724_18_11844_b_2_S4"
)

SAMPLES=(Proband_1 Proband_2 Dad_1 Dad_2 Mum_1 Mum_2)
LANES=(L001 L002 L003 L004)

for SAMPLE in "${SAMPLES[@]}"; do
    PREFIX="${SAMPLE_PREFIX[$SAMPLE]}"
    echo "[fastp] Processing ${SAMPLE}..."

    # Per-lane trimming
    for LANE in "${LANES[@]}"; do
        echo "  [fastp] ${SAMPLE} ${LANE}..."
        fastp \
            --in1  "${RAW_DIR}/${PREFIX}_${LANE}_R1_001.fastq.gz" \
            --in2  "${RAW_DIR}/${PREFIX}_${LANE}_R2_001.fastq.gz" \
            --out1 "${FASTP_DIR}/${SAMPLE}/${SAMPLE}_${LANE}_R1.clean.fastq.gz" \
            --out2 "${FASTP_DIR}/${SAMPLE}/${SAMPLE}_${LANE}_R2.clean.fastq.gz" \
            --html "${FASTP_DIR}/${SAMPLE}/${SAMPLE}_${LANE}_fastp.html" \
            --json "${FASTP_DIR}/${SAMPLE}/${SAMPLE}_${LANE}_fastp.json" \
            --thread 4 \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 35
    done

    # Merge lanes into a single clean FASTQ per sample
    echo "  [fastp] Merging lanes for ${SAMPLE}..."
    cat "${FASTP_DIR}/${SAMPLE}"/*_R1.clean.fastq.gz > "${FASTP_DIR}/${SAMPLE}_R1.fastq.gz"
    cat "${FASTP_DIR}/${SAMPLE}"/*_R2.clean.fastq.gz > "${FASTP_DIR}/${SAMPLE}_R2.fastq.gz"

    echo "  [fastp] ${SAMPLE} done."
done

echo "[STEP 1 COMPLETE] fastp trimming and lane merging finished."
FASTP_BLOCK

# =============================================================================
# STEP 2: bowtie2 — Alignment to hg38
# =============================================================================

echo ""
echo "----------------------------------------------"
echo " STEP 2: bowtie2 — Alignment to hg38"
echo "----------------------------------------------"

bsub -Is -q general-interactive -M ${MEM} -R "rusage[mem=${MEM}]" \
     -a 'docker(quay.io/biocontainers/subread:2.0.1--hed695b0_0)' /bin/bash << 'BOWTIE_BLOCK'

set -euo pipefail

BASE_DIR="${STORAGE1}/scRNAseq_DataSets/kaarunya_files"
FASTP_DIR="${BASE_DIR}/fastp_output"
ALIGN_DIR="${BASE_DIR}/alignments_bt2"
INDEX="${BASE_DIR}/hg38_index"
SAMPLES=(Proband_1 Proband_2 Dad_1 Dad_2 Mum_1 Mum_2)

mkdir -p "${ALIGN_DIR}"
cd "${BASE_DIR}"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "[bowtie2] Aligning ${SAMPLE}..."

    bowtie2 -p 8 \
        -x "${INDEX}" \
        -1 "${FASTP_DIR}/${SAMPLE}_R1.fastq.gz" \
        -2 "${FASTP_DIR}/${SAMPLE}_R2.fastq.gz" \
        -S "${ALIGN_DIR}/${SAMPLE}.bt2.sam" \
        2> "${ALIGN_DIR}/${SAMPLE}.bowtie2.log"

    echo "[bowtie2] ${SAMPLE} alignment complete."
done

echo "[STEP 2 COMPLETE] bowtie2 alignment finished."
BOWTIE_BLOCK

# =============================================================================
# STEP 3: samtools — SAM → Sorted BAM → Index
# =============================================================================

echo ""
echo "----------------------------------------------"
echo " STEP 3: samtools — Sort & Index BAMs"
echo "----------------------------------------------"

bsub -Is -q general-interactive -M ${MEM} -R "rusage[mem=${MEM}]" \
     -a 'docker(quay.io/biocontainers/samtools:1.2-0)' /bin/bash << 'SAMTOOLS_BLOCK'

set -euo pipefail

ALIGN_DIR="${STORAGE1}/scRNAseq_DataSets/kaarunya_files/alignments_bt2"
SAMPLES=(Proband_1 Proband_2 Dad_1 Dad_2 Mum_1 Mum_2)

cd "${ALIGN_DIR}"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "[samtools] Processing ${SAMPLE}..."

    # Convert SAM to sorted BAM
    samtools view -b "${SAMPLE}.bt2.sam" \
        | samtools sort -T "${SAMPLE}.tmp" -o "${SAMPLE}.sorted.bam"

    # Index the sorted BAM
    samtools index "${SAMPLE}.sorted.bam"

    # Remove intermediate SAM to save disk space
    rm -f "${SAMPLE}.bt2.sam"

    echo "[samtools] ${SAMPLE} sorted and indexed."
done

echo "[STEP 3 COMPLETE] samtools sorting and indexing finished."
SAMTOOLS_BLOCK

# =============================================================================
# STEP 4: featureCounts — Generate Gene Counts Matrix
# =============================================================================

echo ""
echo "----------------------------------------------"
echo " STEP 4: featureCounts — Count Reads per Gene"
echo "----------------------------------------------"

bsub -Is -q general-interactive -M ${MEM} -R "rusage[mem=${MEM}]" \
     -a 'docker(quay.io/biocontainers/subread:2.0.1--hed695b0_0)' /bin/bash << 'FEATURECOUNTS_BLOCK'

set -euo pipefail

BASE_DIR="${STORAGE1}/scRNAseq_DataSets/kaarunya_files"
ALIGN_DIR="${BASE_DIR}/alignments_bt2"
COUNTS_DIR="${BASE_DIR}/featurecounts_output"
GTF="/path/to/hg38.gtf"    # UPDATE THIS PATH to your hg38 GTF file
SAMPLES=(Proband_1 Proband_2 Dad_1 Dad_2 Mum_1 Mum_2)

mkdir -p "${COUNTS_DIR}"

# Build list of all sorted BAM files
BAM_FILES=()
for SAMPLE in "${SAMPLES[@]}"; do
    BAM_FILES+=("${ALIGN_DIR}/${SAMPLE}.sorted.bam")
done

echo "[featureCounts] Running on all samples together..."

featureCounts \
    -T 8 \
    -p \
    --countReadPairs \
    -s 0 \
    -a "${GTF}" \
    -o "${COUNTS_DIR}/all_samples_counts.txt" \
    "${BAM_FILES[@]}"

echo "[featureCounts] Counts matrix saved to: ${COUNTS_DIR}/all_samples_counts.txt"
echo "[STEP 4 COMPLETE] featureCounts finished."
FEATURECOUNTS_BLOCK

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

echo ""
echo "=============================================="
echo " PIPELINE COMPLETE"
echo " $(date)"
echo " Final counts matrix:"
echo "   ${COUNTS_DIR}/all_samples_counts.txt"
echo "=============================================="
