library(DESeq2)
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(mosaic)
library(ggh4x)
metacyc_abundance<-read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv",sep="\t",header=T)%>%
  select(-description)

ec_abundance<-read.delim("picrust2/EC_pathways_out/path_abun_unstrat.tsv.gz",sep="\t",header=T)%>%
  select(-description)

# metacyc_abundance<-metacyc_abundance%>%

key<-read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv",sep="\t",header=T)%>%
  select(pathway,description)%>%distinct_all()
metadata <- read_delim("IRB_Human_Metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
metadata

metacyc_abundance
key
m<-metacyc_abundance %>% column_to_rownames("pathway")
dim(m)

m<-m%>%mutate(across(.cols=everything(),as.integer))
m<-m%>%mutate(keep=rowSums(m)>10)%>%filter(keep==TRUE)%>%select(-keep)

dds <- DESeqDataSetFromMatrix(countData = m, colData = metadata, design = ~Diabetes_Status) 
assay(dds)
keep <- rowSums(assay(dds)) >= 10
keep
dds <- dds[keep,]

dds$Diabetes_Status <- relevel(dds$Diabetes_Status, ref = "HC")
# dds$Diabetes_Status <- factor(dds$Diabetes_Status, levels = c("T2D","HC"))
dds <- DESeq(dds,parallel = T)
res <- results(dds,)
write.table(data.frame(res),"T2D_DESeq2_results_full.tsv",sep = "\t")

res2<-data.frame(res)%>%
  filter(padj<0.05)%>%
  filter(abs(log2FoldChange)>0.25)%>%
  filter(!is.na(padj))
dim(res2)

MA_plot<-plotMA(res, ylim=c(-7,7))
# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]
##############################################################################
res2



