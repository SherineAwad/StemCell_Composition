#This is a slight update to code written by Marta Perez-Alcantara
# Extract from vcf the chr, start column, and add an end column
# Filter VCF to only the variants that are different across all donors
# Input: VCF of non-identical SNPs for all donors
# Output: for input SNPs, just chr, pos and pos (start and end are the same, it's 1bp)

library(vcfR)
library(data.table)
library(tidyverse)
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
vcf_path = args[1] #"input_for_create_region_list_bam-readcount/test_vcf.vcf"
output_path = args[2] # "output_for_create_region_list_bam-readcount/test_vcf.txt"
  

# Function to convert gt vcf to minor allele dosages

message("...Reading VCF...")


vcf = vcfR::read.vcfR(vcf_path)

message("...Getting right columns...")

subset = cbind(vcf@fix[,"CHROM"],vcf@fix[,"POS"])
subset=as.data.frame(subset)
print(head(subset))
subset$V2 = as.numeric(subset$V2)
subset$V3 = subset$V2
subset = as.data.table(subset)
fwrite(subset,file = output_path, 
       sep = "\t", col.names = FALSE)
