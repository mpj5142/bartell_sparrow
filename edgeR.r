library("edgeR")
#Change to working directory with R Studio interface

#Load data
biol_reps_data<-read.table("kallisto_master_liver.tsv",stringsAsFactors=FALSE) #Or heart

#Change column and row names to samples and genes, so that data frame cells are just count values
colnames(biol_reps_data)<-biol_reps_data[1,]
rownames(biol_reps_data)<-biol_reps_data[,1]
biol_reps_data<-biol_reps_data[,-1]
biol_reps_data<-biol_reps_data[-1,]

#Convert to matrix and convert characters to numeric type
biol_reps_matrix<-as.matrix(biol_reps_data)
class(biol_reps_matrix)<-"numeric"

#Load data matrix into DGEList format to begin edgeR pipepline
#Groups are as follows: 5 each of non-migratory AM, migratory AM, non-migratory PM, migratory PM; in numerical order
groups<-c("non-mig-am", "mig-pm", "non-mig-am", "non-mig-pm", "mig-pm", "mig-pm", "mig-am", "mig-am", "mig-pm", "mig-am", "mig-am", "mig-am", "mig-pm", "non-mig-pm", "non-mig-pm", "non-mig-pm", "non-mig-am", "non-mig-am", "non-mig-am", "non-mig-pm")
y<-DGEList(biol_reps_matrix,group=groups)

#Filter variants: Recommended defaults here are 5 counts in smallest library, translated to counts/million (here ~25 million = 0.2 cpm), 
#and in at least 5 libraries (or expressed in all members of one group).
#For liver, use 0.15 instead of 0.2 (larger min. library size of 30M)
keep<-rowSums(cpm(y)>0.2)>=5
y<-y[keep, , keep.lib.sizes=FALSE] #Recalculate library sizes afterwards

#Calculate normalization factors
y<-calcNormFactors(y)

#Construct multi-dimensional graph for clustering of samples
colors<-c("red","blue",rep("red",2), rep("blue",9), rep("red",7)) #Color by migratory/non-migratory
points<-c(16, 17, 16, rep(17,3), rep(16,2), 17, rep(16,3), rep(17,4), rep(16,3), 17) #Point shape by morning/evening
plotMDS(y, col=colors, pch=points)

#Legend will need some work
legend("topright", legend=c("Non-migratory day","Migratory day","Non-migratory night","Migratory night"),col=c("red","blue","red","blue"),pch=c(16,16,17,17))
title(main="Multi Dimensional Plot of Sparrow Liver Transcriptome") #Or heart...

#Create the design matrix for all samples
design<-model.matrix(~ 0 + groups)
#Estimate dispersion of variances along model
y <- estimateDisp(y, design, robust=TRUE)
#Model the quasi-lielihood dispersions of variation within model
fit <- glmQLFit(y, design, robust=TRUE)

#Differential expression analysis
diffexp_all<-glmQLFTest(fit, contrast=c(1,1,-1,-1)) #All differences between migrating and non-migrating
diffexp_am<-glmQLFTest(fit, contrast=c(1,0,-1,0)) #Differences between migrating and non-migrating in AM
diffexp_pm<-glmQLFTest(fit, contrast=c(0,1,0,-1)) #Differences between migrating and non-migrating in PM
#Note: we can do other contrasts later as necessary...

#Write as output all (unfiltered) differentialy expressed genes, sorted by FDR. File suffix is either heart or liver
write.csv(topTags(diffexp_all, n=nrow(diffexp_all$table))$table, file="diffexp_all_liver.csv")
write.csv(topTags(diffexp_am, n=nrow(diffexp_am$table))$table, file="diffexp_am_liver.csv")
write.csv(topTags(diffexp_pm, n=nrow(diffexp_pm$table))$table, file="diffexp_pm_liver.csv")
