library("DESeq2")
#Change to working directory with R Studio interface

#Load data
biol_reps_data<-read.table("kallisto_master_heart.tsv",stringsAsFactors=FALSE)

#Change column and row names to samples and genes, so that data frame cells are just count values
colnames(biol_reps_data)<-0:20
rownames(biol_reps_data)<-biol_reps_data[,1]
biol_reps_data<-biol_reps_data[,-1]
biol_reps_data<-biol_reps_data[-1,]

#Convert to matrix and convert characters to numeric type/round to integers
biol_reps_numeric<-as.matrix(biol_reps_data)
class(biol_reps_numeric)<-"numeric"
biol_reps_numeric<-round(biol_reps_numeric)

#Design data frame with condition data for the experiment
#Two separate runs will be one: One (dds_sep) that separates AM and PM samples; one (dds_all) with combined AM/PM samples
#These cannot be done in the same run due to limitations with how DESeq2 factors the desgin matrix; it cannot differentiate between the conditions
condition<-c("non-mig-am", "mig-pm", "non-mig-am", "non-mig-pm", "mig-pm", "mig-pm", "mig-am", "mig-am", "mig-pm", "mig-am", "mig-am", "mig-am", "mig-pm", "non-mig-pm", "non-mig-pm", "non-mig-pm", "non-mig-am", "non-mig-am", "non-mig-am", "non-mig-pm") 
migration<-c("non-mig", "mig", rep("non-mig", 2), rep("mig", 9), rep("non-mig", 7))
sample<-c(1:20)
sample.data.sep<-data.frame(sample,condition)
sample.data.all<-data.frame(sample, migration)

#Load data matrix into DESeq format to begin pipepline. Design will account for infection and time
dds_sep<-DESeqDataSetFromMatrix(countData=biol_reps_numeric, colData=sample.data.sep, design= ~ condition)
dds_all<-DESeqDataSetFromMatrix(countData=biol_reps_numeric, colData=sample.data.all, design= ~ migration)

#Pre-filter to remove rows with <2 total counts; additional filtering will be applied later
dds_sep<-dds_sep[rowSums(counts(dds_sep))>1,]
dds_all<-dds_all[rowSums(counts(dds_all))>1,]

#Perform differential expression analysis
#These steps take a while...
dds_sep<-DESeq(dds_sep)
dds_all<-DESeq(dds_all)

#Generate pairwise comparisions. P-value is set to 0.05, default in DESeq is 0.1
results_all<-results(dds_all, contrast=c("migration", "mig", "non-mig"), alpha=0.05) 	#All infected or non-infected samples
results_am<-results(dds_sep, contrast=c("condition", "mig-am", "non-mig-am"), alpha=0.05) 	#Infected/non-infected in AM
results_pm<-results(dds_sep, contrast=c("condition", "mig-pm", "non-mig-pm"), alpha=0.05)	#Infected/non-infected in PM

summary(results_all)
summary(results_am)
summary(results_pm)

#MA Plot comparing mean expression to log-fold change (quality check)
plotMA(results_all, main="Migratory versus Non-Migratory--All")
plotMA(results_am, main="Migratory versus Non-Migratory--AM")
plotMA(results_pm, main="Migratory versus Non-Migratory--PM")

##Can also get plot counts for individual genes--once annotation occurs, this may be useful. Check documentation

#Export results to text files
write.csv(as.data.frame(results_all), file="results_all_deseq2.csv")
write.csv(as.data.frame(results_am), file="results_am_deseq2.csv")
write.csv(as.data.frame(results_pm), file="results_pm_deseq2.csv")

#Export normalized counts to files
write.csv(counts(dds_all, normalized=TRUE), file="counts_all_norm_deseq2.csv")
write.csv(counts(dds_sep, normalized=TRUE), file="counts_sep_norm_deseq2.csv")

#######Working on visualization later...
