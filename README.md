Snakemake Workflow for Stem Cell Project
==========================================================================


We use mixture of cells to estimate cell fractions. 
The project includes variant calling and stem cell fraction estimators. 

In this version, we call variants from REFERENCE using a simple samtools pipeline. 
We find out  common region list of variants between cells and reference using bam-count. 
Then we apply the estimator to get cells' fraction 





