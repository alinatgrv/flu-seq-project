#!/usr/bin/env bash
set -euo pipefail

REF="data/raw/KF848938.1.fasta"
FASTQ="data/raw/SRR1705851.fastq.gz"
OUT_BAM="data/processed/roommate.sorted.bam"

echo "Indexing reference (if not indexed yet)..."
bwa index "${REF}"

echo "Aligning roommate sample..."
bwa mem "${REF}" "${FASTQ}" \
  | samtools view -S -b - \
  | samtools sort -o "${OUT_BAM}" -

echo "Indexing BAM..."
samtools index "${OUT_BAM}"

echo "Calculating mapping statistics..."
TOTAL=$(samtools view -c "${OUT_BAM}")
UNMAPPED=$(samtools view -c -f4 "${OUT_BAM}")
MAPPED=$((TOTAL - UNMAPPED))

echo "Total reads in BAM: ${TOTAL}"
echo "Unmapped reads:    ${UNMAPPED}"
echo "Mapped reads:      ${MAPPED}"

