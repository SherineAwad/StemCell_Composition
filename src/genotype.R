args <- commandArgs(trailingOnly = TRUE)
source("/home/ec2-user/stemcell/pipeline/src/getProp.R")


# allele specific count table at 6453 SNPs
Y=as.matrix(read.table(args[1], as.is=T, header=T))

# genotype dosage of 11 donors at the 6453 SNPs
G=as.matrix(read.table(args[2], as.is=T, header=T))

# fitting the model
res=getProp(Y,G)

# result
res

# true propotions 
ptruth=scan(args[3])
pdf(args[4])
# visualising the estimated proportions with confidence intervals (y-axis) aginst the truth (x-axis)
plot(ptruth, res$P[,1], xlim=c(0,0.2), ylim=c(0,0.2))
segments(ptruth, res$P[,2], ptruth, res$P[,3])


