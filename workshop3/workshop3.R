#workshop3.R
#R source file for workshop 3

#Load DESeq2 and pheatmap
require(DESeq2)
require(pheatmap)

#Reads in the matrix of log-transformed gene expression counts, and prints
log_counts = read.table("counts.txt")
print(head(log_counts))
print(dim(log_counts)) #prints the amount of genes in the dataset (51414)
#If there's a significant difference between gene expression in healthy and diseased individuals, then we should see that samples of similar types are more correlated with one another than with samples of different types.
# Compare counts in two normal samples
plot(log_counts[,1], log_counts[,3])
abline(0,1)
print(cor(log_counts[,1], log_counts[,3]))

dev.new()

#Compare counts in a healthy sample and a diseased sample
plot(log_counts[,1], log_counts[,12])
abline(0,1)
print(cor(log_counts[,1], log_counts[,12]))

dev.new()

#Compute the correlation between all different samples
print(cor(log_counts))

#View matrix produced above as heatmap
log_cor <- cor(log_counts)
diag(log_cor) <- NA
pheatmap(log_cor, show_rownames = TRUE)
dev.new()

##Dealing with batch effects -- patterns in sequencing data caused by way sample is handled before/during sequencing (eg. the environment the samples are left in, different sequencing runs). Can prove to be problematic if they overpower the actual biological signals we're looking for. Options to handle this / minimize risk:
#1. Loading covariates
covariates = read.table("covariates.txt")
#in linux: run 
#view covariates.txt

#2. Principal Component Analysis -- an algorithm that detects systematic changes/variation in the sequencing data to detect batch effects. Works on the premise that the same set of genes will be affected similarly across affected samples. PCA identifies these genes that vary, and outputs 'principal components' that represent the presence of varying genes. 
pca = prcomp(log_counts)
#Plotting by disease status
# A simple function for color-coding points
color_code <- function(x, vals, colors)
{
  i = which(x == vals)
  return(colors[i])
}
# Plot PC1 and PC2
color_status = sapply(substring(covariates$samples, 1, 3), FUN=color_code,
                      vals=c("IPF", "Nor"), colors=c("blue", "red"))
plot(pca$rotation[,1], pca$rotation[,2], col=color_status, pch=16)

#Plotting by sequencing batch
color_seqbatch = sapply(covariates$seq.batch, FUN=color_code,
                        vals=c(1, 2), colors=c("blue", "red"))
plot(pca$rotation[,1], pca$rotation[,2], col=color_seqbatch, pch=16)
dev.new()
##The above show that there is a strong clustering effect caused by the sequencing batch, as shown by principal component #2

# Plot PC3 and PC4
plot(pca$rotation[,3], pca$rotation[,4], col=color_status, pch=16)
dev.new()

#Re-plot heatmap after removing first 2 principal components --> why do we have to remove the first one?? 
## ASK STEPHEN about understanding principal components -- what exactly they are, how do we decide to get rid of which
pca_mat = cor(t(as.matrix(pca$rotation[,3:6])))
diag(pca_mat) = NA
pheatmap(pca_mat)

##DIFFERENTIAL EXPRESSION -- DESeq2 looks across samples for differentially expressed genes; the BIOS201 class has already given us the original read count data with simulated (above) sequencing batch effects removed. 
#Question: How would we remove the batch effects once we've identified them? (we only identified them in prev section)
counts = read.table("counts_denoised.txt")
rownames(covariates) = covariates$samples
dds <- DESeqDataSetFromMatrix(countData <- counts,
                              colData <- covariates,
                              design = ~ factor(status, levels=c("Norm","IPF")))
dds <- DESeq(dds, betaPrior=FALSE)
res_nosex = results(dds)
plotMA(dds)

#ASK STEPHEN to explain the DESeq2 graph
dev.new()
#Correct for individual's sex (covariate)
dds <- DESeqDataSetFromMatrix(countData <- counts,
                              colData <- covariates,
                              design = ~ sex + factor(status, levels=c("Norm","IPF")))
dds <- DESeq(dds, betaPrior=FALSE)
res = results(dds)
plotMA(dds)

#Compare the strength of results between the sex-corrected DESeq2 run and the run with no covariates corrected.
plot(log(res_nosex$pvalue, base=10), log(res$pvalue, base=10), xlim=c(-30,0), ylim=c(-30,0),
            xlab="p-value, No sex correction", ylab="p-value, sex correction")
abline(a=0,b=1)
abline(a=-5,b=0,lty=3,col="red")
abline(v=-5,lty=3,col="red")

#see how many genes are differentially expressed
print(sum(res_nosex$padj < 0.05, na.rm=TRUE))
print(sum(res$padj < 0.05, na.rm=TRUE))

print(head(res[order(res$pvalue),]))
status <- c(rep("Norm", 7), rep("IPF", 8))
boxplot(as.numeric(log_counts[rownames(log_counts)=="ENSG00000170962",]) ~ status)

#MUC5B Gene - ENSG00000117983
# DSP Gene - ENSG00000096696

