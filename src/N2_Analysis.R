library(dplyr)
library(tidyverse)
library(DESeq2)
library(tximport)
library(biomaRt)
library(ggplot2)
library(SummarizedExperiment)
library(stringr)
library(ggfortify)
library(data.table)
source('src/utils.R')

# get tph1 + n2 data
rds = readRDS('nxt/salmon/salmon_merged_gene_counts.rds')
rds$Names = rds$names %>% gsub('\\.[0-9].*', '', .)

# annotate metadata:
colData(rds)['Treatment'] = sapply(rds$Names, treatment_getter) %>% as.factor
colData(rds)['Genotype'] = sapply(rds$Names, genotype_getter) %>% as.factor
colData(rds)['Age'] = sapply(rds$Names, time_getter) %>% as.factor

# put into DESeq
dds50 = loadDDS(rds[,rds$Age == '50' & rds$Genotype == 'wt'], which='wt',
                design=~Treatment)
res50 = get_results(dds50, 'Treatment', 'ascr', 'cnt', filter=FALSE)

dds58 = loadDDS(rds[,rds$Age == '58' & rds$Genotype == 'wt'], which='wt',
                design=~Treatment)
res58 = get_results(dds58, 'Treatment', 'ascr', 'cnt', filter=FALSE)

res50 %>% as.data.frame %>% write.csv('../DE_N250.csv', quote = FALSE)
res58 %>% as.data.frame %>% write.csv('../DE_N258.csv', quote = FALSE)