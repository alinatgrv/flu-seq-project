# flu-seq-project

# Influenza HA Deep Sequencing Analysis

This repository contains the analysis workflow for a practical assignment on influenza viral
evolution, sequencing error control, and detection of antigenically relevant mutations in the
hemagglutinin (HA) gene.  
All steps follow the instructions provided in the course handout.

## Repository structure

```
flu-seq-project/
├── data/
│ ├── raw/ # Reference sequence and raw fastq files (roommate + controls)
│ └── processed/ # BAM files, mpileup files
├── results/
│ ├── roommate/ # Variant calls for roommate sample
│ └── controls/ # Variant calls for three control samples
├── scripts/ # Bash and Python scripts used in analysis
└── report/
└── lab_report.md # Full lab report with methods, results, figures, and conclusions
```


## Overview of the analysis

1. **Download data**  
   Roommate sample and three control samples were downloaded from NCBI SRA (fastq format).  
   Reference sequence KF848938.1 was obtained from NCBI GenBank.

2. **Alignment**  
   All samples were aligned to the reference using `bwa mem`, sorted and indexed with `samtools`.

3. **Variant calling**  
   - High-frequency variants (`≥95%`) were detected using VarScan (`--min-var-freq 0.95`).  
   - Rare variants were detected using VarScan (`--min-var-freq 0.001`).

4. **Error estimation**  
   Variant frequencies from control samples were used to estimate the background sequencing
   error rate (mean + 3 × SD). A Python script (`scripts/error_stats.py`) identified
   high-confidence mutations in the roommate sample.

5. **Epitope mapping**  
   High-confidence mutations were compared with known HA epitope regions (from Muñoz et al.).
   Two mutations were found within antigenic epitopes.

6. **Report**  
   The complete analysis, explanation of all steps, figures, and conclusions are included in  
   `report/lab_report.md`.

## Scripts included

- `scripts/align_roommate.sh` – alignment pipeline for the roommate sample  
- `scripts/align_controls.sh` – alignment pipeline for three controls (optional)  
- `scripts/error_stats.py` – analysis of background error rates and identification of
  high-confidence mutations

## Requirements

- `bwa`  
- `samtools`  
- `VarScan`  
- Python 3 + pandas, numpy  
- IGV (optional, for visualization)

## How to reproduce the analysis

Clone the repository and run:

```bash
bash scripts/align_roommate.sh
bash scripts/align_controls.sh
python scripts/error_stats.py
```

The resulting variant tables will appear in the results/ directory.

## Author

A. Tagirova 