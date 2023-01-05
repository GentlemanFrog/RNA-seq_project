
# RNA-seq project to perform gene expression analysis with the DESeq2 library

# The data comes from the paper:
# Human Airway Smooth Muscle Transcriptome Changes in Response to Asthma
# Medications (GSE52778)

# Design of the study:

# 1. RNA-seq was performed on 4 primary human airway smooth muscle cell lines treated with
# 1 micromolar dexamethasone (drug with a strong and effective anti-inflammatory, anti-allergic and immunosuppressive effect)
# for 18  hours

# 2. For each of the 4 cell lines, study has a treated and untreated sample

# 3. The aim of the study: It is to understand the transcriptional changes occurring due to the treatment with used drug

# Data download and csv preperation:
# To download library needed for this experiment the following command must be used:
# BiocManager::install("airway")

setwd("C:/Users/Alek/Desktop/RNA-seq_project")
library(DESeq2)
library(tidyverse)
library(airway)

data(airway)
airway

# Changing the data to data_frame
sample_info <- as.data.frame(colData(airway))

# Restricting the data only to interesting columns 
sample_info <- sample_info[,c(2,3)]

# Replacing shortcuts with proper name for column dex
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)

# Renaming the column names: cell, dex to more descriptive names 
names(sample_info) <- c('cellLine', 'dexamethasone')

# Writing chosen data info into csv file 
write.table(sample_info, file = "sample_info.csv", sep = ',',
            col.names = T, row.names = T, quote = F)

# Preparing data of experiment where rows are genes, columns are counts for each sample and gene
countsData <- assay(airway)

# Writing count data info into csv file
write.table(countsData, file = "counts_data.csv", sep = ',',
            col.names = T, row.names = T, quote = F)



# Chapter 1: Preparing the data:

# Reading in counts data from csv:
counts_data<- read.csv('counts_data.csv')

# we see the rows are gene IDs (Ensemble), columns are the sample name but by this names
# we don't know which samples are treated and which one are untreated which we will set up
head(counts_data)

# Reading the sample_info:
column_data<-read.csv('sample_info.csv')
head(column_data)

# To perform calculations with DESeq we have to check the row names of column_data
# matches the names in counts_data.
# We performing this by checking if all colnames in counts_data are present in column_data
# We get the information in the boolean if there are all of them or not
all(colnames(counts_data) %in% rownames(column_data))

# In addition we check if they are in the same order:
all(colnames(counts_data) == rownames(column_data))



# Chapter 2: The construction of DESeqDataSet object:

# First argument countData refers to counts genes accross samples, colData refers to data about cell line of each sample
# that we are checking and information about it being treated with drug or not, design is a formula which express
# how the counts for each gene depend on the variables in colData, so for our example we want to check the differences
# between treated and untreated samples so design refers to the dexamethasone column
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = column_data,
                       design = ~ dexamethasone)
dds

# Pre-filtering: allows us to remove the rows with low gene count across all the samples
# We are removing rows that have at least 10 reads total

# Keeped_rows variable will contain all genes with the true or false vale of our argument
# which allows to filter the dds object

# This step allows us to reduce the memory size of the dds object, and increase the speed
# of the transformation and testing functions with DESeq2. It can also improve visualization,
# as features with no information for differential expression will not be plotted
keeped_rows <- rowSums(counts(dds)) >= 10
dds <- dds[keeped_rows,] # we drooped ~ 40k rows with lower than 10 counts
dds

# Setting the factor level: This steps is for telling the DESeq which sample is treated and which is not
# In our example we have 2 levels: treated and untreated. So we have to tell the DESeq that we want to use the
# untreated samples as reference level, to compare treated sample with the untreated. 

dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")
dds$dexamethasone

# If this part will not be set the DESeq will alphabetically choose which level is reference

# Its good to remember if the samples have technical replicates we should collapse them before running DESeq
# Technical replicates are repeated measurements of the same sample that represent independent measures of
# the random noise associated with protocols or equipment (https://www.nature.com/articles/nmeth.3091). 

# Chapter 3: Running the DESeq functions

# Now the script will perform the standard differential expression analysis, which all steps are wrapped in the single 
# function DESeq. 
dds <- DESeq(dds)

# Saving results from DESeqData object:
res <- results(dds)

# As we can see in results, values of log2 FC are calculated as we set them in our design factor,
# so all values we see are based on the treated samples. 

# What the different columns mean?
# baseMean - avg. of the normalized counts taken over all the samples

# log2FC - is the FC value of current gene (from row) compared with the untreated value, so the positive values
#         means that gene is up-regulated in treated sample, and negative values means that they are down-regulated in the terated conditions

# lfcSE - provides the standard error estimates for the log2FC

# stat - this column represent the Wald test values for checked genes. This test ia a parametric 
#       statistical measure to confirm whether a set of independent variables are collecively 'significant' 
#       for model or not. 
# pvalue - is basically the p_value of the the test statistics

# padj - this column contains adjusted p_values for multiple testing. We have to correct the p_values
#       for multiple testing because whenever the statistical test is performed we use a p value of 0.05.
#       5% of our differential expressed genes are not really differential expressed and they are there due to random chance.
#       To reduce problem of this false positives we perform methods to calculate the adjusted p_values to avoid this detection
#       of false positive genes. 
res



# Chapter 4: Exploration of the results

# In the summary we can see how many genes are up-regulated, down-regulated, how much of them are outlires
# and how many of them have low counts. 
summary(res)

# Also we can change the threshold of adjusted p_values we look at
res_p_0.01 <- results(dds, alpha = 0.01)
summary(res_p_0.01)

# MA plot:

# This plot is the plot of log2FC vs the mean of the normalized counts. The genes that are differentially
# expressed are colored in blue. This genes are significantly differentially expressed genes in our example
# also there are visible triangles to the top of the plot, this are genes with the higher foldchanges, and the direction 
# triangles tell the direction of the fold change. In this type of plots we want to see genes in the upper right 
# part of the plot, which means that these will have high mean of normalized counts and high log2FC values
# which will make these genes candidates to be further look into. 
plotMA(res)

