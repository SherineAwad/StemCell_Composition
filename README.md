Snakemake Workflow for Stem Cell Project
==========================================================================


We use mixture of cells to estimate cell fractions. 
The project includes variant calling and stem cell fraction estimators. 

In this version, we call variants from REFERENCE using a simple samtools pipeline. 
We find out  common region list of variants between cells and reference using bam-count. 
Then we apply the estimator to get cells' fraction 

# Workflow Summary: From VCF to Donor Proportion Estimates

## Step 1: Extract SNP Positions from VCF
- **Input:** VCF file containing SNPs for multiple donors.  
- **Process:**  
  - Extract chromosome (`CHROM`) and position (`POS`) of each SNP.  
  - Add an `end` column equal to `POS` (since SNPs are single-base).  
  - Only positions where **at least one donor has a variant** are kept.  
- **Output:** Regions file with SNP positions:

| CHROM | start | end |
|-------|-------|-----|
| chr1  | 100   | 100 |
| chr1  | 200   | 200 |

**Purpose:** Create a list of variable positions to target for read-level analysis.

---

## Step 2: Count Reads at SNP Positions (`bam_count` rule)
- **Input:**  
  - Donor BAM file (`{sample}.sorted.bam`)  
  - Regions file from Step 1 (`{cohort}_regions.txt`)  
- **Process:**  
  - For each SNP position, count how many sequencing reads contain **A, C, G, T**.  
  - Collect additional metrics such as total reads, base quality, strand bias, deletions (`DEL`), and ambiguous bases (`N`).  
- **Output:** `bam-readcount` file (per donor):

| CHROM | POS | REF | TOTAL_READS | A | C | G | T | DEL | N |
|-------|-----|-----|-------------|---|---|---|---|-----|---|
| chr1  | 100 | A   | 14          | 2 | 3 | 4 | 5 | 0   | 0 |
| chr1  | 200 | G   | 10          | 0 | 1 | 7 | 2 | 0   | 0 |

**Purpose:** Measure the actual nucleotide composition at each SNP for the donor.

---

## Step 3: Estimate Minor Allele Frequency (`b_estimate`)
- **Input:**  
  - `bam-readcount` file (donor-specific)  
  - Original VCF file with SNPs and genotypes  
- **Process:**  
  - Match SNPs between bam-readcount and VCF.  
  - Convert donor genotypes in the VCF into **minor allele dosages**:  
    - `0` → homozygous reference  
    - `0.5` → heterozygous  
    - `1` → homozygous minor allele  
  - Calculate **minor allele frequency** for each SNP in the donor:

\[
b\_estimate = \frac{\text{reads supporting ALT}}{\text{total reads at that position}}
\]

- **Output:**  
  - `b_estimate.csv` → donor-specific minor allele frequency per SNP

| CHROM | POS | REF | TOTAL_READS | A | C | G | T | ALT | b_estimate |
|-------|-----|-----|-------------|---|---|---|---|-----|------------|
| chr1  | 100 | A   | 14          | 2 | 3 | 4 | 5 | G   | 0.286      |
| chr1  | 200 | G   | 10          | 0 | 1 | 7 | 2 | T   | 0.2        |

- **Minor Allele Dosage Table:** contains **all donors** with 0, 0.5, 1 values per SNP:

| SNP        | Donor1 | Donor2 | Donor3 | … |
|------------|--------|--------|--------|---|
| chr1_100_A | 0.5    | 1      | 0      | … |
| chr1_200_G | 0      | 0.5    | 1      | … |

---

## Step 4: Estimate Donor Proportions (`w_estimate`)
- **Input:**  
  - `b_estimate.csv` → minor allele frequencies for the sample/donor  
  - Minor allele dosage table (`genotype_minor_allele_dosage.csv`)  
- **Process:**  
  - Solve the linear system \( A \cdot w \approx b \) where:  
    - \( A \) = minor allele dosage matrix (SNPs × donors)  
    - \( w \) = vector of donor proportions to estimate  
    - \( b \) = observed minor allele frequencies from bam-readcount  
  - Constraints:  
    - \( 0 \le w_i \le 1 \) for each donor  
    - \( \sum w_i = 1 \)  
  - SNPs where all donors are identical or not observed in the reads are removed.  
- **Output:**  
  - `w_estimate.txt` → estimated proportion of each donor in the sample

| Donor  | Estimated Proportion |
|--------|--------------------|
| Donor1 | 0.25               |
| Donor2 | 0.50               |
| Donor3 | 0.25               |

**Purpose:** Determine the fraction of each donor contributing to the sample based on observed minor allele frequencies.

---

### **Summary of Data Flow**
1. **VCF → regions file** → list of SNP positions.  
2. **BAM + regions → bam-readcount** → per-donor read counts at each SNP.  
3. **Bam-readcount + VCF → b_estimate & minor allele dosages** → minor allele frequencies for donor and dosages for all donors.  
4. **b_estimate + dosages → w_estimate** → estimated donor proportions in the sample.




