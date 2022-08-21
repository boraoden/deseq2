setwd("C:\\Users\\borao\\Desktop\\SciLifeLab\\DESeq")

library(readxl)
library(dplyr)
library(tximport)
library(DESeq2)
library(ggplot2)

# read count data 
hss_ac_raw <- read.delim("hss_all_Count.txt", header = TRUE)

# set rownames of count data dataframe 
hss_ac <- read.delim("hss_all_Count.txt", header = TRUE, row.names = hss_ac_raw$Ensemble)

# remove a unnecessary row
hss_ac$Ensemble <- NULL

# read metadata
hss_md <- read_excel("hss_metadata.xlsx", col_names = TRUE)

# move two columns in metadata
hss_md_2 <- hss_md %>% relocate(PatientNo, .before = PatientENAno, .after = NULL) %>% relocate(Sarcopenia, .before = PatientENAno, .after = NULL)

# write count data and metadata to a file
# write.table(hss_ac, file = "hss_countdata.csv" , sep = ",")
# write.table(hss_md_moved, file = "hss_metadata_moved.csv" , sep = ",")

# call rows in metadata (trash code below)
# hss_md_trimmed <- hss_md_moved[-(c(1:320) - c(1,9,17)), ]
# hss_md_trimmed <- hss_md_moved[-c((sapply(c(1:320),"-",c(1,9,17)))), ]
# hss_md_trimmed <- hss_md_moved[ %in% 1:320, ]
# hss_md_trimmed <- hss_md_moved[intersect(c(1,9,17), c(1:320)), ]

# call rows in metadata (my for loop)
called_rows <- NULL
for (i in 0:39){
  row_numbers <- 8*i + 1
  called_rows <- c(called_rows, row_numbers)
  hss_md_3 <- hss_md_2[called_rows, ]
}

# set row names of metadata as PatientNo
hss_md_4 <- as.data.frame(hss_md_3)
row.names(hss_md_4) <- hss_md_4$PatientNo

# easy way
# hss_md_called_easy <- hss_md_moved[seq(1,313,by=8), ]
# View(hss_md_called_easy)

#DESeq 
dds <- DESeqDataSetFromMatrix(countData=round(hss_ac), colData=hss_md_4, design= ~Sarcopenia, tidy = FALSE)
dds <- DESeq(dds)

# Results & Summary

dds_results <- as.data.frame(results(dds))
head(results(dds, tidy=TRUE))
summary(dds_results)

# sort results by adjusted p value
dds_results_sorted_padj <- dds_results[order(dds_results$padj),]
head(dds_results_sorted_padj)

# sort results by p value
dds_results_sorted_pvalue <- dds_results[order(dds_results$pvalue),]
head(dds_results_sorted_pvalue)

# write results
write.table(as.data.frame(dds_results), file = "C:\\Users\\borao\\Desktop\\SciLifeLab\\DESeq\\DESeqResults.txt", quote = F, sep = '\t')

# plot counts

# compare the normalized counts between treated and control groups for our top 6 genes (adjusted p value)
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000170624", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000162688", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000081189", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000115221", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000163171", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000182253", intgroup="Sarcopenia")

# compare the normalized counts between treated and control groups for our top 6 genes (p value)
par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000285343", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000170624", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000188517", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000143502", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000188452", intgroup="Sarcopenia")
plotCounts(dds, gene="ENSG00000162688", intgroup="Sarcopenia")

# other plots
plotCounts(dds, gene="ENSG00000170624", intgroup="Sarcopenia", returnData=TRUE) %>% 
  ggplot(aes(Sarcopenia, count)) + geom_boxplot(aes(fill=Sarcopenia)) + scale_y_log10() + ggtitle("SGCD")

# volcano plot (pvalue)

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(dds_results_sorted_pvalue, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(dds_results_sorted_pvalue, pvalue<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_results_sorted_pvalue, pvalue<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# PCA

#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Sarcopenia") #using the DESEQ2 plotPCA fxn we can


# adjusted p value check
dds_results_2_padj <- as.data.frame(dds_results) %>% dplyr::mutate(sig=padj<0.05)
# How many of each?
dds_results_2_padj %>% 
  group_by(sig) %>% 
  summarize(n=n())