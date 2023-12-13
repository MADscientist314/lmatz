library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(mosaic)
library(ggh4x)
getwd()
# setwd("/home/jochum00/lmatz/LM_RawSequences")
# If you want to analyze KEGG pathway abundance instead of KO within the pathway, turn ko_to_kegg to TRUE.
# KEGG pathways typically have more explainable descriptions.

# Load metadata as a tibble
# data(metadata)
metadata <- read_delim("IRB_Human_Metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
metadata
# Load KEGG pathway abundance
kegg_abundance <- ko2kegg_abundance("picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv") 
kegg_abundance
# "ALDEx2": ANOVA-Like Differential Expression tool for high throughput sequencing data 
# "DESeq2": Differential expression analysis based on the negative binomial distribution using DESeq2 
# "edgeR": Exact test for differences between two groups of negative-binomially distributed counts using edgeR 
# "limma voom": Limma-voom framework for the analysis of RNA-seq data 
# "metagenomeSeq": Fit logistic regression models to test for differential abundance between groups using metagenomeSeq 
# "LinDA": Linear models for differential abundance analysis of microbiome compositional data 
# "Maaslin2": Multivariate Association with Linear Models (MaAsLin2) for differential abundance analysis 
# "Lefse": Linear discriminant analysis (LDA) effect size algorithm for high-dimensional microbiome data
# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR","limma voom","metagenomeSeq","LinDa","Maaslin2","Lefse")
################################################################################
# aldex2_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
#                                    metadata = metadata, 
#                                    group = "Diabetes_Status",
#                                    p.adjust ="BH",
#                                    daa_method = "ALDEx2",
#                                    select = NULL,
#                                    reference = "HC") 
# 
# aldex2_daa_results_df%>%filter(p_values<0.05) #28
# aldex2_daa_results_df%>%filter(p_adjust<0.05) #0
################################################################################
deseq2_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
                                   metadata = metadata, 
                                   group = "Diabetes_Status",
                                   p.adjust ="BH",
                                   daa_method = "DESeq2",
                                   select = NULL,
                                   reference = "HC") 

deseq2_daa_results_df<-deseq2_daa_results_df%>%filter(!is.na(p_values))) #14
deseq2_daa_results_df%>%filter(p_adjust<0.05) #4
deseq2_daa_results_df
################################################################################
edger_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
                                  metadata = metadata, 
                                  group = "Diabetes_Status",
                                  p.adjust ="BH",
                                  daa_method = "edgeR",
                                  select = NULL,
                                  reference = "HC") 

edger_daa_results_df%>%filter(p_values<0.05) #14
edger_daa_results_df%>%filter(p_adjust<0.05) #11
################################################################################
# limma_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
#                                   metadata = metadata, 
#                                   group = "Diabetes_Status",
#                                   p.adjust ="BH",
#                                   daa_method = "limma voom",
#                                   select = NULL,
#                                   reference = "HC") 
# 
# dim(limma_daa_results_df%>%filter(p_values<0.05)) #21
# limma_daa_results_df%>%filter(p_adjust<0.05) #0
################################################################################
# metagenomeSeq_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
#                                           metadata = metadata, 
#                                           group = "Diabetes_Status",
#                                           p.adjust ="BH",
#                                           daa_method = "metagenomeSeq",
#                                           select = NULL,
#                                           reference = "HC") 
# 
# metagenomeSeq_daa_results_df%>%filter(p_values<0.05) #3
# metagenomeSeq_daa_results_df%>%filter(p_adjust<0.05) #0
################################################################################
# LinDa_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
#                                   metadata = metadata, 
#                                   group = "Diabetes_Status",
#                                   p.adjust ="BH",
#                                   daa_method = "LinDA",
#                                   select = NULL,
#                                   reference = "HC") 
# 
# LinDa_daa_results_df%>%filter(p_values<0.05) #14
# LinDa_daa_results_df%>%filter(p_adjust<0.05) #0
################################################################################
# Maaslin2_daa_results_df<- pathway_daa(abundance = kegg_abundance, 
#                                      metadata = metadata, 
#                                      group = "Diabetes_Status",
#                                      p.adjust ="BH",
#                                      daa_method = "Maaslin2",
#                                      select = NULL,
#                                      reference = "HC") 
# 
# Maaslin2_daa_results_df%>%filter(p_values<0.05) #16
# Maaslin2_daa_results_df%>%filter(p_adjust<0.05) #0
################################################################################
# Run ggpicrust2 with input file path

# Annotate pathway results using KO to KEGG conversion
deseq2_daa_results_df
deseq2_daa_results_df_slim<-deseq2_daa_results_df%>%filter(p_values<0.05)

# a<-deseq2_daa_results_df[1:50,]
deseq2_daa_results_df%>%filter(p_adjust<0.25)
# a2 <- pathway_annotation(pathway = "KO", daa_results_df = a, ko_to_kegg = TRUE)
edger_daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = edger_daa_results_df, ko_to_kegg = TRUE,)
deseq2_daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = deseq2_daa_results_df, ko_to_kegg = F,)
deseq2_daa_annotated_sub_method_results_df

kegg_abundance
deseq2_daa_results_df




edger_daa_annotated_sub_method_results_df

deseq2_daa_annotated_sub_method_results_df$p_values
deseq2_daa_annotated_sub_method_results_df
metadata<-metadata%>%mutate(Diabetes_Status=factor(Diabetes_Status,levels=c("T2D","HC")))

metadata$Diabetes_Status
# Generate pathway error bar plot
# Please change Group to metadata$your_group_column if you are not using example dataset
p <- pathway_errorbar(abundance = kegg_abundance, 
                     daa_results_df = deseq2_daa_annotated_sub_method_results_df, 
                     Group = metadata$Diabetes_Status, 
                     p_values_threshold = 0.25, 
                     order = "pathway_class", 
                     select = NULL, 
                     ko_to_kegg = T, 
                     p_value_bar = TRUE, 
                     colors = NULL, 
                     x_lab = "pathway_name")
p$data


glimpse(deseq2_daa_annotated_sub_method_results_df)


write.table(p$data,"picrust")
################################################################################

# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
ko_abundance
p <- pathway_errorbar(abundance = ko_abundance %>% 
                       column_to_rownames("#NAME"), 
                     daa_results_df = deseq2_daa_annotated_sub_method_results_df, 
                     Group = metadata$Diabetes_Status, 
                     p_values_threshold = 0.05, 
                     order = "Diabetes_Status",
                     select = deseq2_daa_annotated_sub_method_results_df %>% 
                       arrange(p_adjust) %>% 
                       slice(1:20) %>% dplyr::select(feature) %>% pull(), 
                     ko_to_kegg = T, 
                     p_value_bar = TRUE, 
                     colors = NULL)

# Workflow for MetaCyc Pathway and EC

# Load MetaCyc pathway abundance and metadata
data("metacyc_abundance")
metacyc_abundance2<-read_table("picrust2/pathways_out/path_abun_unstrat.tsv")


data("metadata")
ggpicrust2::metacyc_abundance
metacyc_abundance<-read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv",sep="\t",header=T)%>%select(-description)
library(mosaic)
dim(metacyc_abundance)
df<-metacyc_abundance%>%rowwise()%>%mutate(tot=sum(Buffington_201_001:Buffington_201_014))
df<-df%>%arrange(tot)
df_keep<-df%>%filter(tot>10)
df_keep
metacyc_abundance<-metacyc_abundance%>%filter(pathway%in%df_keep$pathway)
dim(metacyc_abundance)

# Perform pathway DAA using LinDA method
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance %>% 
                                      column_to_rownames("pathway"), 
                                      metadata = metadata, 
                                      group = "Diabetes_Status",
                                      p.adjust = "BH", 
                                      daa_method = "DESeq2")
metacyc_daa_results_df%>%filter(p_adjust<0.05)
# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc",
                                                       daa_results_df = metacyc_daa_results_df,
                                                       ko_to_kegg = FALSE)
write.table(metacyc_daa_annotated_results_df,"pathway_daa_deseq2_results_df.tsv",sep = "\t",row.names = F)
metacyc_daa_annotated_results_df$method<-"DESeq2"
unique( metacyc_daa_annotated_results_df$method)
unique(metacyc_daa_annotated_results_df$group2)
m<-metacyc_abundance %>% column_to_rownames("pathway")
# Generate pathway error bar plot
# Please change column_to_rownames() to the feature column
# Please change Group to metadata$your_group_column if you are not using example dataset
pathway_errorbar(abundance = m,
                daa_results_df = metacyc_daa_annotated_results_df, 
                Group = "group2", 
                ko_to_kegg = FALSE, 
                p_values_threshold = 0.05,
                order = "group2", 
                select = NULL, 
                p_value_bar = TRUE, 
                colors = NULL) 
                # x_lab = "description")
metacyc_daa_results_df
# Generate pathway heatmapmetacyc_daa_results_df
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust <0.05)
dim(feature_with_p_0.05)
metacyc_daa_annotated_results_df$feature
metacyc_abundance$pathway
h<- metacyc_abundance %>% filter(pathway %in% feature_with_p_0.05$feature)

h

h<-left_join(h,metacyc_daa_annotated_results_df)
h2<-h[,2:12]
h2
rownames(h2)<-h$pathway




pathway_heatmap(abundance =h2, 
                metadata = metadata,
                group = "Diabetes_Status")


library(pheatmap)
mat<-as.matrix(h2)
mat
metacyc_daa_annotated_results_df
#######################
annot_col<-metacyc_daa_annotated_results_df%>%filter(feature%in%rownames(mat))
rownames(annot_col)<-annot_col$feature
annot_col<-annot_col %>% select(-c(feature,method,group1,group2,adj_method))
annot_col
metacyc_daa_results_df
metacyc_abundance
library(DESeq2)
metacyc_abundance<-read.delim("picrust2/pathways_out/path_abun_unstrat_descrip.tsv",sep="\t",header=T)%>%select(-description)
m<-metacyc_abundance %>% column_to_rownames("pathway")
m<-m%>%mutate(across(.cols=everything(),as.integer))
m
dds <- DESeqDataSetFromMatrix(m, colData = metadata, design = ~Diabetes_Status) 
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$Diabetes_Status <- relevel(dds$Diabetes_Status, ref = "HC")
# dds$Diabetes_Status <- factor(dds$Diabetes_Status, levels = c("T2D","HC"))
dds <- DESeq(dds,parallel = T)
res <- results(dds,)

res2<-data.frame(res)%>%filter(padj<0.05)%>%filter(!is.na(padj))
res2

MA_plot<-plotMA(res, ylim=c(-7,7))
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

class(MA_plot)

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pathway_pca(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment")

# Run pathway DAA for multiple methods
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
methods <- c("ALDEx2", "DESeq2", "edgeR")
daa_results_list <- lapply(methods, function(method) {
 pathway_daa(abundance = metacyc_abundance %>% column_to_rownames("pathway"), metadata = metadata, group = "Environment", daa_method = method)
})

# Compare results across different methods
comparison_results <- compare_daa_results(daa_results_list = daa_results_list, method_names = c("ALDEx2_Welch's t test", "ALDEx2_Wilcoxon rank test", "DESeq2", "edgeR"))