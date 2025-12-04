#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="data/raw"
mkdir -p "${RAW_DIR}"

echo "Downloading roommate sample (SRR1705851)..."
wget -O "${RAW_DIR}/SRR1705851.fastq.gz" \
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/SRR1705851.fastq.gz"

echo "Downloading control 1 (SRR1705858)..."
wget -O "${RAW_DIR}/SRR1705858.fastq.gz" \
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz"

echo "Downloading control 2 (SRR1705859)..."
wget -O "${RAW_DIR}/SRR1705859.fastq.gz" \
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz"

echo "Downloading control 3 (SRR1705860)..."
wget -O "${RAW_DIR}/SRR1705860.fastq.gz" \
  "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz"

echo "Downloading reference HA sequence KF848938.1..."
curl \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KF848938.1&rettype=fasta&retmode=text" \
  -o "${RAW_DIR}/KF848938.1.fasta"

echo "All downloads finished."
