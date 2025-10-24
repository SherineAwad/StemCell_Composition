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
  - `b_estimate.csv` → observed minor allele frequencies in the sample (from sequencing reads)
  - Minor allele dosage table (`genotype_minor_allele_dosage.csv`) → known genotypes for all donors

- **Process (conceptually):**
  - For each SNP, the estimator compares the observed minor allele frequency in the sample to the known minor allele dosages of each donor.
  - It determines a set of donor proportions that best explains the observed frequencies across all SNPs.
  - SNPs where all donors have the same genotype or positions not observed in the sequencing data are ignored.
  - The solution ensures that:
    - Each donor proportion is between 0% and 100%
    - The total of all donor proportions sums to 100%

- **Output:**
  - `w_estimate.txt` → estimated contribution of each donor to the sample


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



# Conceptual Explanation of the Donor Proportion Estimator

---

### What the estimator does

Imagine you have a **mixed sample** that contains DNA from several donors. You want to figure out **how much each donor contributed** to that sample.

To do this, the estimator uses:

1. **Observed data from the sequencing reads (`b_estimate`)**:

   * At each variable position (SNP), it knows what fraction of reads in the mixed sample carry the minor allele.
   * This tells you, for each SNP, how common the minor allele is in the mixture.

2. **Donor genotypes (`genotype_minor_allele_dosage`)**:

   * For every donor, you know whether they have 0, 1, or 2 copies of the minor allele at each SNP.
   * This is essentially the **donor’s “signature”** for each SNP.

---

### How it works conceptually

* For each SNP, the estimator asks:  
  *“Given that I see this fraction of minor alleles in the sample, and knowing each donor’s SNP profile, what combination of donors could produce this observation?”*

* It looks across all SNPs and tries to find a **set of donor proportions** that best explains the observed minor allele frequencies.

* It ensures that the proportions:

  * Are between 0 and 100% for each donor
  * Add up to 100% in total

---

### What it produces

* A simple **list of donor contributions** (weights/proportions) for the sample:

| Donor  | Estimated Contribution |
| ------ | ---------------------- |
| Donor1 | 25%                    |
| Donor2 | 50%                    |
| Donor3 | 25%                    |

* This tells you **how much DNA each donor contributed** to the sample.

---

✅ **In short:**  
The estimator takes the observed DNA mixture (from sequencing reads) and compares it to known donor genotypes to **work out the fraction of each donor in the mixture**. It’s like figuring out the recipe of a cake by tasting it and knowing the possible ingredients.



