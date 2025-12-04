# Flu Deep Sequencing Project — Lab Notebook  
### Progress Log (Steps 1–3)

---

## 1. Downloading and inspecting raw data

### 1.1. Samples and files

| Sample | SRA ID | Description |
|--------|--------|--------------|
| Roommate sample | **SRR1705851** | Viral H3N2 sample from sick roommate; used to detect mutations |
| Control 1 | SRR1705858 | Isogenic reference control (should match reference exactly) |
| Control 2 | SRR1705859 | Isogenic reference control |
| Control 3 | SRR1705860 | Isogenic reference control |
| Reference | **KF848938.1** | HA gene sequence of A/Hong Kong/4801/2014 (H3N2) |

All files were downloaded into `data/raw/` using a custom script.

---

### 1.2. Inspecting the roommate FASTQ

Checked the FASTQ structure:

`gzcat data/raw/SRR1705851.fastq.gz | head`


FASTQ format is correct (4-line entries).

Number of reads:

`expr $(gzcat SRR1705851.fastq.gz | wc -l) / 4`


### 1.3. Basic observations
- single-end Illumina reads  
- high coverage  
- good apparent base quality  

---

## 2. Mapping to the reference sequence

### 2.1. Reference indexing

`bwa index data/raw/KF848938.1.fasta`


---

### 2.2. Aligning the roommate sample

`bwa mem data/raw/KF848938.1.fasta data/raw/SRR1705851.fastq.gz | samtools view -S -b - | samtools sort -o data/processed/roommate.sorted.bam -`


Indexing BAM:

`samtools index data/processed/roommate.sorted.bam`


---

### 2.3. Mapping statistics

Total aligned reads:

`samtools view -c roommate.sorted.bam`
361349


Unmapped reads:

`samtools view -c -f4 roommate.sorted.bam`

233


Summary:

- **Mapped:** 361116  
- **Total:** 361349  
- **Mapping rate:** **99.94%**

Conclusion: reference is correct; HA gene matches sample well.

---

## 3. Generating mpileup

`samtools mpileup -d 100000 -f data/raw/KF848938.1.fasta data/processed/roommate.sorted.bam > data/processed/roommate.mpileup `



First lines confirm correct format and coverage (300–1100×).

---

## 4. Detecting common (high-frequency) variants

### 4.1. Running VarScan (min freq = 0.95)

`varscan mpileup2snp data/processed/roommate.mpileup --min-var-freq 0.95 --output-vcf 1 > results/roommate/roommate_common.vcf` 


VarScan summary:

1665 bases in pileup file
5 variant positions (5 SNP, 0 indel)
5 SNPs reported


---

### 4.2. Counting known variants

`grep -v '^#' roommate_common.vcf | wc -l`

5


### 4.3. Extracting simplified variant table

`grep -v '^#' roommate_common.vcf | awk '{print $1, $2, $4, $5}' > roommate_common_simple.tsv`


This produced a 4-column table (CHROM, POS, REF, ALT) for 5 high-frequency SNPs.

---

## Current summary

- Raw data successfully downloaded  
- Roommate sample aligned with **99.94% mapping**  
- mpileup generated  
- VarScan detected **5 high-frequency SNPs**  
- Variant table extracted and saved  

Next steps:
- Rare variant calling (min freq = 0.001)  
- Processing three control samples  
- Error modeling and identifying real rare variants  
- Epitope mapping  

---

## 5. Detection of rare variants in the roommate sample (min freq = 0.001)

After identifying high-frequency variants, the next task was to search for rare mutations that may represent minor quasispecies within the viral population.

### 5.1. Running VarScan with a rare-variant frequency threshold (0.001)

The following command was used to call rare SNPs:

```bash
varscan mpileup2snp data/processed/roommate.mpileup --min-var-freq 0.001 --output-vcf 1 > results/roommate/roommate_rare.vcf
```

VarScan output summary:

1665 bases in pileup file
23 variant positions (21 SNP, 2 indel)
0 were failed by the strand-filter
21 variant positions reported (21 SNP, 0 indel)

This result shows that lowering the frequency threshold significantly increases the number of detected positions compared to the 5 high-frequency variants.

### 5.2. Counting rare variants

`grep -v '^#' results/roommate/roommate_rare.vcf | wc -l`

21

So, 21 rare SNPs were detected at a 0.1% frequency threshold.

### 5.3. Inspecting raw variant calls

The first several lines of the rare-variant VCF:

`grep -v '^#' results/roommate/roommate_rare.vcf | head`

Frequencies are stored in the sixth field of the sample column (GT:GQ:DP:RD:AD:FREQ:...).

### 5.4. Extracting CHROM, POS, REF, ALT, FREQ into a clean table

A simplified rare-variant table was generated using:

`grep -v '^#' results/roommate/roommate_rare.vcf | awk -F '\t' '{ split($10,a,":"); print $1, $2, $4, $5, a[6] }' > results/roommate/roommate_rare_simple.tsv`

This produced a table containing 21 rare variants, with their chromosome, position, reference base, alternate base, and estimated frequency.


---

## 5. Inspection and alignment of control samples

To distinguish true rare mutations in the roommate sample from sequencing or PCR errors, three technical control datasets (isogenic reference virus) were analyzed. These samples should match the reference sequence exactly; thus any detected “mutations” represent background error rates of the sequencing workflow.

### 6.1. Read counts in control FASTQ files

The number of lines in each FASTQ file was counted using `gzcat` and `wc -l`:

```bash
gzcat data/raw/SRR1705858.fastq.gz | wc -l   # 1026344
gzcat data/raw/SRR1705859.fastq.gz | wc -l   # 933308
gzcat data/raw/SRR1705860.fastq.gz | wc -l   # 999856
```

| Sample ID  | Total lines | Reads             |
| ---------- | ----------- | ----------------- |
| SRR1705858 | 1,026,344   | **256,586 reads** |
| SRR1705859 | 933,308     | **233,327 reads** |
| SRR1705860 | 999,856     | **249,964 reads** |


### 6.2. Reference sequence length

The length of the HA reference (KF848938.1) was determined as:

```grep -v '^>' data/raw/KF848938.1.fasta | tr -d '\n' | wc -c```

1665

The HA gene is 1665 bp long.

### 6.3. Rough coverage estimate

Coverage was approximated using:

```
coverage ≈ number of reads × read length / 1665
```

Given Illumina single-end reads (~150–250 bp), each control dataset provides coverage on the order of 20,000–40,000×, sufficient for detecting low-frequency sequencing errors.

### 6.4. Alignment of control samples

Each control FASTQ file was aligned to the HA reference using bwa mem, piped directly into samtools for BAM conversion and sorting:

``` bwa mem data/raw/KF848938.1.fasta data/raw/SRR1705858.fastq.gz | samtools view -S -b - | samtools sort -o data/processed/control1.sorted.bam -```

```bwa mem data/raw/KF848938.1.fasta data/raw/SRR1705859.fastq.gz | samtools view -S -b - | samtools sort -o data/processed/control2.sorted.bam - ```

``` bwa mem data/raw/KF848938.1.fasta data/raw/SRR1705860.fastq.gz | samtools view -S -b - | samtools sort -o data/processed/control3.sorted.bam - ```

Each resulting BAM file was indexed:

```
samtools index data/processed/control1.sorted.bam
samtools index data/processed/control2.sorted.bam
samtools index data/processed/control3.sorted.bam
```

### 6.5. Mapping statistics for control samples

Mapping statistics were collected using samtools view -c:

Control 1 (SRR1705858)

```
Total reads:     256744
Unmapped reads:      86
Mapped reads:    256658
Mapping rate:    99.97%
```

Control 2 (SRR1705859)

```
Total reads:     233451
Unmapped reads:      76
Mapped reads:    233375
Mapping rate:    99.97%
```

Control 3 (SRR1705860)

```
Total reads:     250184
Unmapped reads:      76
Mapped reads:    250108
Mapping rate:    99.97%
```


All three samples show >99.9% mapping efficiency, confirming that the sequencing reads originate from a viral strain identical to the reference sequence. This provides a clean dataset for estimating the background sequencing error profile.

---

## 7. Rare variant calling in control samples using VarScan

To quantify the background sequencing error rate, rare variants were called in each of the three control (isogenic reference) samples using the same HA reference sequence (KF848938.1) and the same parameters as for the roommate sample, but with a very low frequency cutoff (0.1%).

### 7.1. Generating mpileup files for control BAMs

For each control alignment, an mpileup file was generated with a high depth limit:

```bash
samtools mpileup -d 100000 \
  -f data/raw/KF848938.1.fasta \
  data/processed/control1.sorted.bam \
  > data/processed/control1.mpileup

samtools mpileup -d 100000 \
  -f data/raw/KF848938.1.fasta \
  data/processed/control2.sorted.bam \
  > data/processed/control2.mpileup

samtools mpileup -d 100000 \
  -f data/raw/KF848938.1.fasta \
  data/processed/control3.sorted.bam \
  > data/processed/control3.mpileup
```

Each mpileup contains coverage and base information for all 1665 positions of the HA gene.

### 7.2. Calling rare variants with VarScan (min freq = 0.001)

VarScan was run on each mpileup file with a minimum variant frequency of 0.001 (0.1%), requesting VCF output:

```bash
varscan mpileup2snp data/processed/control1.mpileup \
  --min-var-freq 0.001 \
  --output-vcf 1 \
  > results/controls/control1_rare.vcf

varscan mpileup2snp data/processed/control2.mpileup \
  --min-var-freq 0.001 \
  --output-vcf 1 \
  > results/controls/control2_rare.vcf

varscan mpileup2snp data/processed/control3.mpileup \
  --min-var-freq 0.001 \
  --output-vcf 1 \
  > results/controls/control3_rare.vcf
```

VarScan summaries:
```
Control 1:
1665 bases in pileup file
58 variant positions (58 SNP, 0 indel)
1 were failed by the strand-filter
57 variant positions reported (57 SNP, 0 indel)

Control 2:
1665 bases in pileup file
54 variant positions (54 SNP, 0 indel)
2 were failed by the strand-filter
52 variant positions reported (52 SNP, 0 indel)

Control 3:
1665 bases in pileup file
61 variant positions (61 SNP, 0 indel)
0 were failed by the strand-filter
61 variant positions reported (61 SNP, 0 indel)
```

Counting non-header lines in each VCF confirms the number of SNPs:

```bash
grep -v '^#' results/controls/control1_rare.vcf | wc -l   
grep -v '^#' results/controls/control2_rare.vcf | wc -l  
grep -v '^#' results/controls/control3_rare.vcf | wc -l
```
Thus, the three control samples contain 57, 52, and 61 rare SNP calls at a 0.1% threshold, respectively. These calls mainly represent sequencing / PCR errors.

### 7.3. Parsing VCF files into simple variant tables

VarScan reports the detailed genotype information in the sample column using the format:

```
GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:...
```

Here FREQ (6th field) is the estimated variant allele frequency.
For each control sample, the VCF file was parsed to extract chromosome, position, reference base, alternative base, and frequency:

```bash
grep -v '^#' results/controls/control1_rare.vcf \
  | awk -F '\t' '{ split($10,a,":"); print $1, $2, $4, $5, a[6] }' \
  > results/controls/control1_rare_simple.tsv

grep -v '^#' results/controls/control2_rare.vcf \
  | awk -F '\t' '{ split($10,a,":"); print $1, $2, $4, $5, a[6] }' \
  > results/controls/control2_rare_simple.tsv

grep -v '^#' results/controls/control3_rare.vcf \
  | awk -F '\t' '{ split($10,a,":"); print $1, $2, $4, $5, a[6] }' \
  > results/controls/control3_rare_simple.tsv
```

Each resulting *_rare_simple.tsv file contains:
### Parsed rare variants

| CHROM       | POS  | REF | ALT | FREQ    |
|-------------|------|-----|-----|---------|
| KF848938.1  | 254  | A   | G   | 0.17%   |
| KF848938.1  | 307  | C   | T   | 0.94%   |
| KF848938.1  | 389  | T   | C   | 0.22%   |
| KF848938.1  | 691  | A   | G   | 0.17%   |
| KF848938.1  | 722  | A   | G   | 0.20%   |
| ...         | ...  | ... | ... | ...     |


## 8. Comparison of rare variants between controls and roommate sample

### 8.1 Processing control variant frequencies

Rare variants were extracted from all three control samples using VarScan  
(`--min-var-freq 0.001`) and parsed into tables containing:

- reference base (`REF`)
- position (`POS`)
- alternative base (`ALT`)
- variant frequency (`FREQ`)

Frequencies were converted from percentages to fractions.

### 8.2 Error statistics in control samples

To estimate the background sequencing error rate, the mean and standard deviation of
variant frequencies were computed for each control sample.

Using the script `scripts/error_stats.py`, the following values were obtained:

Per-control error statistics (fraction units):

```
sample mean_error sd_error
control1 0.002565 0.000717
control2 0.002369 0.000524
control3 0.002503 0.000780
```


Thresholds for error outliers were defined as:

```
mean_error + 3 × sd_error
```

```
control1: 0.004717
control2: 0.003941
control3: 0.004844
```


We used the **maximum of these three thresholds** as the conservative cutoff:
```
cutoff = 0.004844
```


### 8.3 Identification of high-confidence mutations in roommate sample

Rare variants from the roommate sample were compared against the control-based
cutoff. Any variant with frequency > 0.004844 was considered a true mutation
rather than sequencing noise.

The following variants exceeded the cutoff:

### High-confidence variants in roommate sample

| CHROM       | POS  | REF | ALT | FREQ     |
|-------------|------|-----|-----|----------|
| KF848938.1  | 72   | A   | G   | 0.9996   |
| KF848938.1  | 117  | C   | T   | 0.9982   |
| KF848938.1  | 307  | C   | T   | 0.0094   |
| KF848938.1  | 774  | T   | C   | 0.9996   |
| KF848938.1  | 999  | C   | T   | 0.9986   |
| KF848938.1  | 1260 | A   | C   | 0.9994   |
| KF848938.1  | 1458 | T   | C   | 0.0084   |


### 8.4 Interpretation

- Several mutations (72, 117, 774, 999, 1260) are present at ~100% frequency,
  indicating that the patient's virus is genetically distinct from the reference strain.
- Two low-frequency variants (positions 307 and 1458) also exceed the 3×SD threshold,
  suggesting possible within-host viral evolution or mixed viral populations.
- These mutations should be examined in IGV to determine:
  - the codon affected
  - the amino acid change
  - the position in the HA protein sequence


## 9. Epitope mapping

To determine whether any of the high-confidence mutations in the roommate’s influenza
sample occur within known antigenic epitopes of hemagglutinin (HA), we used epitope
definitions from Muñoz & Deem (2004).  

### 9.1 Epitope definitions (H3 numbering)

- **Epitope A:** 122, 124, 126, 130–133, 135, 137, 138, 140, 142–146, 150, 152, 168  
- **Epitope B:** 128, 129, 155–160, 163, 165, 186–190, 192–194, 196–198  
- **Epitope C:** 44–48, 50, 51, 53, 54, 273, 275, 276, 278–280, 294, 297, 299–300,  
  304–305, **307–312**  
- **Epitope D:** 96, 102, 103, **117**, 121, 167, 170–177, 179, 182, 201, 203,  
  207–209, 212–219, 226–230, 238, 240, 242, 244, 246–248  
- **Epitope E:** 57, 59, 62–63, 67, 75, 78, 80–83, 86–88, 91–92, 94, 109, 260–262, 265

### 9.2 High-confidence mutations in the roommate sample

Variants considered high-confidence were those with frequencies exceeding  
the control-derived threshold (**mean + 3×SD**):

| CHROM       | POS  | REF | ALT | FREQ     |
|-------------|------|-----|-----|----------|
| KF848938.1  | 72   | A   | G   | 0.9996   |
| KF848938.1  | 117  | C   | T   | 0.9982   |
| KF848938.1  | 307  | C   | T   | 0.0094   |
| KF848938.1  | 774  | T   | C   | 0.9996   |
| KF848938.1  | 999  | C   | T   | 0.9986   |
| KF848938.1  | 1260 | A   | C   | 0.9994   |
| KF848938.1  | 1458 | T   | C   | 0.0084   |

### 9.3 Epitope-mapped mutations

Two mutations fall within known epitope regions of hemagglutinin:

| POS  | REF→ALT | Epitope | Interpretation |
|------|----------|---------|----------------|
| **117** | C→T | **D** | Clear epitope mutation; position explicitly listed in epitope D. |
| **307** | C→T | **C** | Within 307–312 block of epitope C. |

### 9.4 Interpretation

- Mutations at **positions 117 (epitope D)** and **307 (epitope C)** occur within
  immunodominant regions of HA and may influence antibody recognition.
- All other high-confidence mutations occur outside epitope regions and therefore
  are unlikely to affect epitope-mediated immune escape.
- These findings suggest the roommate’s viral strain contains at least two mutations
  in antigenically relevant sites.

---

## 10. Conclusion

In this project, we used deep sequencing data from one infected individual (“roommate”) and
three reference control samples to identify true mutations in the influenza hemagglutinin (HA)
segment. Our goal was to understand why a vaccinated person could still become infected with
the circulating viral strain.

### 10.1 Why vaccination did not prevent infection

By comparing the roommate's viral population with reference sequences, we identified several
high-confidence mutations. Two of these substitutions occurred within well-characterized
antigenic epitope regions of the hemagglutinin (HA) protein:

- **Position 117 → epitope D**
- **Position 307 → epitope C**

Mutations within HA epitopes can disrupt antibody binding by modifying the shape, charge,
or accessibility of antigenic surfaces. Modern studies on influenza immunity show that even
single amino-acid substitutions in key antigenic regions can substantially alter antibody
recognition and promote immune escape (Koel et al., 2013). Such antigenic drift is one of
the main reasons influenza vaccines vary in effectiveness across seasons and viral subtypes
(Belongia et al., 2016). Reviews of the human antibody response further emphasize that HA is
under continuous immune pressure, and the accumulation of epitope-directed mutations enables
circulating strains to evade pre-existing vaccine-induced immunity (Krammer, 2019).

Therefore, the most likely explanation for the infection is that the roommate’s virus carried
epitope-altering mutations that diminished the effectiveness of the prior influenza vaccine.


### 10.2 Which mutations are real?

To distinguish real biological variants from sequencing noise, we quantified the background
error rate using the three reference samples. For each, we calculated:

- mean variant frequency  
- standard deviation  

We then defined a conservative threshold:

**variant frequency > mean + 3 × SD (across all controls)**

Any mutation exceeding all three thresholds was considered high-confidence. This approach
allows us to filter out:

- random base calling errors  
- polymerase (PCR) errors  
- low-frequency sequencing artifacts  

Using this method, we identified **7 real mutations**, two of which fall within HA epitopes.

### 10.3 Why error control matters

Deep sequencing errors can easily appear as low-frequency variants. Without filtering, these
artifacts could be misinterpreted as evidence of viral evolution, mixed infections, or drug
escape. Accurate error control is therefore essential for:

- identifying true viral mutations  
- quantifying variant frequencies  
- understanding within-host viral diversity  

### 10.4 Ways to improve error control

Below are additional strategies that would strengthen confidence in variant calls:

#### **(1) Use Unique Molecular Identifiers (UMIs)**
UMIs tag each original RNA molecule before amplification. Reads sharing the same UMI can be
collapsed into a consensus sequence, eliminating PCR-introduced mutations.  
This substantially reduces artificial variants, especially at low frequencies.

#### **(2) Use higher-fidelity polymerases**
Switching to a polymerase with 50–100× lower error rates reduces the number of
PCR-generated substitutions before sequencing even begins. This decreases background noise and
makes true low-frequency variants easier to detect.

#### **(3) Apply bioinformatic error-correction tools**
Software such as **LoFreq**, **iVar**, or **DeepVariant** models position-specific sequencing
errors and produces statistically robust low-frequency variant calls. These approaches outperform simple frequency cutoffs.

#### **(4) Increase sequencing depth**
More coverage reduces stochastic sampling error, improving confidence in true rare variants.

#### **(5) Compare multiple independent replicates**
True biological variants will appear in all replicates at similar frequencies, while
sequencing artifacts will not reproduce consistently.

Each of these steps would reduce noise and improve the reliability of rare-variant detection.

---

### References

Belongia, E. A., Skowronski, D. M., McLean, H. Q., Chambers, C., Sundaram, M. E., & De Serres, G. (2016).
**Variable influenza vaccine effectiveness by subtype: a systematic review and meta-analysis of test-negative design studies.**
*The Lancet Infectious Diseases*, 16(8), 942–951. https://doi.org/10.1016/S1473-3099(16)00129-8

Krammer, F. (2019).
**The human antibody response to influenza A virus.**
*Nature Reviews Immunology*, 19, 383–397. https://doi.org/10.1038/s41577-019-0143-6

Koel, B. F., et al. (2013).
**Substitutions near the receptor-binding site determine major antigenic change during influenza virus evolution.**
*Science*, 342(6161), 976–979. https://doi.org/10.1126/science.1244730



