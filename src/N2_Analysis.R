library(dplyr)
library(tidyverse)
library(DESeq2)
library(tximport)
#library(biomaRt)
library(ggplot2)
library(SummarizedExperiment)
library(stringr)
library(ggfortify)
library(data.table)
setwd('~/repos/MalePheromoneRNAseq/src')
source('utils.R')

# get tph1 + n2 data
rds = readRDS('~/repos/MalePheromoneRNAseq/nxt2/salmon/salmon.merged.gene_counts.rds')
rds$Names = rds$names %>% gsub('\\.[0-9].*', '', .)

# annotate metadata:
colData(rds)['Treatment'] = sapply(rds$names, treatment_getter) %>% as.factor
colData(rds)['Genotype'] = sapply(rds$Names, genotype_getter) %>% as.factor
colData(rds)['Age'] = sapply(rds$Names, time_getter) %>% as.factor
colData(rds)['Batch'] = sapply(rds$names, batch_getter) %>% as.factor

# put into DESeq
dds50 = loadDDS(rds[,rds$Age == '50'], which='wt', design=~Treatment)
res50 = get_results(dds50, 'Treatment', 'ascr', 'cnt', filter=FALSE)
s50 = lfcShrink(dds50, 'Treatment_ascr_vs_cnt')

# at time 58
dds58 = loadDDS(rds[,rds$Age == '58'], which='wt', design=~Treatment + Batch)
# dds58 = loadDDS(rds[,(rds$Age == '58') & (rds$Batch == 'batch2')], which='wt',
#                design=~Treatment)
res58 = get_results(dds58, 'Treatment', 'ascr', 'cnt', filter=FALSE)
s58 = lfcShrink(dds58, 'Treatment_ascr_vs_cnt')


# compare gene expression thru time
# ddsTime = loadDDS(rds[,rds$Genotype == 'wt'], which='wt', design=~Age)
# resTime = get_results(ddsTime, 'Age', '58', '50', filter=FALSE)


res50 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_N250.csv', quote = FALSE)
res58 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_N258.csv', quote = FALSE)
s50 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_N250_shrunken.csv', quote = FALSE)
s58 %>% as.data.frame %>% write.csv('../data/diff_exp/DE_N258_shrunken.csv', quote = FALSE)


dds = loadDDS(rds, which='', design=~Age + Treatment + Batch)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Age", "Treatment"))