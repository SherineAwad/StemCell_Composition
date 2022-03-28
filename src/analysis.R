library(edgeR)
#Read in raw counts from 3 lanes
data<-read.delim("raw_count.tsv",row.names=1)
data2<-read.delim("raw_count_lane2.tsv",row.names=1)
data3<-read.delim("raw_count_lane3.tsv",row.names=1)

#Create edgeR DGE object
y<-DGEList(counts=data+data2+data3)

#Sample info
meta<-read.delim("~/genehunter/fastq/SLX-16082.HWMJKBBXX.s_7.contents.csv",sep=",")
meta2<-read.delim("samples.txt")
samples<-merge(meta,meta2,by.x="Sample.name",by.y="Sample",sort=T)
samples<-samples[order(samples$Barcode),]

#Groups
y$samples$group<-relevel(as.factor(samples$Time.vs.fast),ref="Pre")
fast<-relevel(as.factor(samples$Time.vs.fast),ref="Pre")
gender<-as.factor(samples$Gender)
subject<-as.factor(samples$Participant)

#Remove low expressed genes
keep<-rowSums(cpm(y)>=5)>=13
y<-y[keep,]
y$samples$lib.size<-colSums(y$counts)

#TMM Normalisation
y<-calcNormFactors(y)

#Model dispersions
design<-model.matrix(~subject+fast)
y<-estimateGLMCommonDisp(y,design,verbose=T)
y<-estimateGLMTagwiseDisp(y,design)
y<-estimateGLMTrendedDisp(y,design)

#Fit
fit<-glmFit(y,design)

#LRTs
noGender_post_vs_pre_lrt<-glmLRT(fit,coef=14)
#write.table(topTags(noGender_post_vs_pre_lrt,n="inf"),"noGender_post_vs_pre_lrt.txt",sep="\t")

#interaction with sex 

#Create new model for interaction
design2<-model.matrix(~gender*fast)
y<-estimateGLMCommonDisp(y,design2,verbose=T)
y<-estimateGLMTrendedDisp(y,design2)
y<-estimateGLMTagwiseDisp(y,design2)

#Fit and LRT
interaction_fit<-glmFit(y,design2)
interaction_lrt<-glmLRT(interaction_fit,coef=4)

#write.table(topTags(interaction_lrt,n="inf"),"sex_difference.txt",sep="\t")
