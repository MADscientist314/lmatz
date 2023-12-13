library(DESeq2)
library(phyloseq)
# first we need to make a DESeq2 object

sam<-data.frame(sample_data(ps))
sam$Diabetes_Status
counts<-data.frame(otu_table(ps))
deseq_counts <- DESeqDataSetFromMatrix(counts, colData = sam, design = ~Diabetes_Status) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
deseq_counts_vst
vst_trans_count_tab <- assay(deseq_counts_vst)
vst_trans_count_tab
# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sam)
tax
ps
#############################################
tax_phy<-tax_table(ps)
sam_phy<-sample_data(ps)
tree_phy<-phy_tree(ps)
vst_physeq <- phyloseq(vst_count_phy,tax_phy,sam_phy,tree_phy)
meta(vst_physeq)
# generating and visualizing the PCoA with phyloseq
unweighted_vst_pcoa <- ordinate(vst_physeq, 
                                method="MDS", 
                                distance="unifrac")
weighted_vst_pcoa <- ordinate(vst_physeq, 
                              method="MDS", 
                              distance="wunifrac")

vst_pcoa <- ordinate(ps, method="MDS",distance="bray")
vst_pcoa$values
sample_names(ps)
unweighted_vst_pcoa$vectors
# UniFrac(vst_physeq, 
#         weighted=TRUE, 
#         normalized=TRUE, 
#         parallel=FALSE, 
#         fast=TRUE)

eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
color_pal<-list(Diabetes_Status=c("sT2D"="#FF7100","HC"="#240A94"))
my_pal<-c("#240A94","#FF7100")

vst_physeq
plot_ordination(vst_physeq, vst_pcoa, color="Diabetes_Status") +
  geom_point(size=1) +
  labs(col="Diabetes_Status") +
  # geom_text(aes(label=rownames(sam), hjust=0.3, vjust=-0.4)) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  scale_color_manual(values=my_pal)+
  theme_bw() +
  theme(legend.position="none")
vst_pcoa

################################################################################
res<-vst_pcoa$vectors
res
res<-data.frame(res)%>%mutate(SampleID=rownames(res))
sam
res
library(tidyverse)
res2<-full_join(res,sam)

res2<-res2%>%mutate(BMI=as.numeric(BMI))
res2<-as_tibble(res2)
res2<-res2%>%mutate(Sex=factor(Sex,levels=c("Male","Female")))
res2$BMI
res2
library(ggpubr)
ggscatter(data = res2,
          x = "Axis.1",
          y = "Axis.2",
          color = "BMI",
          rug = F) +scale_color_viridis_c()
          # shape = "Sex",
          # xlab = "Axis 1 [23.6%]",
          # ylab = "Axis 2 [13.9%]",
          # title = "Unweighted Unifrac",
          # palette = my_pal,
          ellipse = T,
          ellipse.border.remove = F,
          ellipse.alpha = 0,ellipse.level = 0.95,
          star.plot = T
          # add = "reg.line"
          ) +
  # coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  theme_pubr(legend = "right") +
  scale_color_viridis_c()
#############################WEIGHTED ###################################################
eigen_vals <- unweighted_vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
color_pal<-list(Diabetes_Status=c("sT2D"="#FF7100","HC"="#240A94"))
my_pal<-c("#240A94","#FF7100")


plot_ordination(vst_physeq, weighted_vst_pcoa, color="Diabetes_Status") +
  geom_point(size=1) +
  labs(col="Diabetes_Status") +
  # geom_text(aes(label=rownames(sam), hjust=0.3, vjust=-0.4)) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  scale_color_manual(values=my_pal)+
  theme_bw() +
  theme(legend.position="none")
vst_pcoa
################################################################################
res<-weighted_vst_pcoa$vectors
res
res<-data.frame(res)%>%mutate(SampleID=rownames(res))
res2<-full_join(res,sam)
res2<-res2%>%mutate(BMI=as.numeric(BMI))
res2<-as_tibble(res2)
res2$SampleID
res2<-res2%>%mutate(Sex=factor(Sex,levels=c("Male","Female")))
ggscatter(data = res2,
          x = "Axis.1",
          y = "Axis.2",
          color = "Diabetes_Status",
          rug = F,
          shape = "Sex",
          xlab = "Axis 1 [85.8%]",
          ylab = "Axis 2 [13%]",
          title = "Weighted Unifrac",
          palette = my_pal,
          # label = "SampleID",
          ellipse = T,
          ellipse.border.remove = F,
          ellipse.alpha = 0,ellipse.level = 0.95,
          star.plot = T
          # add = "reg.line"
) + coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  theme_pubr(legend = "right") #+
# scale_color_viridis_c()
################################################################################
res<-weighted_vst_pcoa$vectors
res
res<-data.frame(res)%>%mutate(SampleID=rownames(res))
res2<-full_join(res,sam)
res2<-res2%>%mutate(BMI=as.numeric(BMI))
res2<-as_tibble(res2)
res2$SampleID
res2<-res2%>%mutate(Sex=factor(Sex,levels=c("Male","Female")))
ggscatter(data = res2,
          x = "Axis.1",
          y = "Axis.2",
          color = "Diabetes_Status",
          rug = F,
          shape = "Sex",
          xlab = "Axis 1 [85.8%]",
          ylab = "Axis 2 [13%]",
          title = "Weighted Unifrac",
          palette = my_pal,
          # label = "SampleID",
          ellipse = T,
          ellipse.border.remove = F,
          ellipse.alpha = 0,ellipse.level = 0.95,
          star.plot = T
          # add = "reg.line"
) + coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  theme_pubr(legend = "right") #+
# scale_color_viridis_c()


ggscatter(data = res2,
          x = "Axis.1",
          y = "Axis.2",
          color = "Diabetes_Status",
          rug = F,
          shape = "Sex",
          xlab = "Axis 1 [85.8%]",
          ylab = "Axis 2 [13%]",
          title = "Weighted Unifrac",
          palette = my_pal,
          # label = "SampleID",
          ellipse = T,
          ellipse.border.remove = F,
          ellipse.alpha = 0,ellipse.level = 0.95,
          star.plot = T
          # add = "reg.line"
) + coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  theme_pubr(legend = "right") #+
# scale_color_viridis_c()

colSums(counts)



library(vegan)
vst_trans_count_tab

adonis2(vst_trans_count_tab~Diabetes_Status,data = sam,permutations = 999)
euc_dist

uni<-UniFrac(vst_physeq,weighted = FALSE,normalized = FALSE,parallel = T)
wuni<-UniFrac(vst_physeq,weighted = TRUE,normalized = FALSE,parallel = T)

anova(betadisper(uni, sam$Diabetes_Status)) # 0.002
u<-adonis2(uni~Diabetes_Status,data = sam,permutations = 999) # 0.003
summary(u)
u

adonis2(wuni~Diabetes_Status,data = sam,permutations = 9999) # 0.003


betadisper(wuni, sam$Diabetes_Status)
anova(betadisper(wuni, sam$Diabetes_Status)) # 0.002
anova(betadisper(uni, sam$Diabetes_Status)) # 0.002

gghistogram(data = sam,x = "BMI",fill = "Diabetes_Status")

adonis2(wuni~Diabetes_Status*BMI,data = sam,permutations = 9999) # 0.003

sam
