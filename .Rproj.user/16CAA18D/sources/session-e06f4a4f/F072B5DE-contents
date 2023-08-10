library(DESeq2)
# first we need to make a DESeq2 object
sam$Diabetes_Status
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

vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
meta(vst_physeq)
# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
color_pal<-list(Diabetes_Status=c("sT2D"="#FF7100","HC"="#240A94"))
my_pal<-c("#FF7100","#240A94")


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

res<-data.frame(res)%>%mutate(SampleID=rownames(res))
res2<-full_join(res,sam)
res2<-res2%>%mutat(BMI=as.numeric(BMI))
res2<-as_tibble(res2)
res2
ggscatter(data = res2,
          x = "Axis.1",
          y = "Axis.2",
          color = "BMI",
          rug = F,
          shape = "Diabetes_Status",
          xlab = "Axis 1 [23.6%]",ylab = "Axis 2 [13.9%]",
          # palette = my_pal,
          ellipse = F,
          ellipse.border.remove = F,
          ellipse.alpha = 0,ellipse.level = 0.95,
          star.plot = T
          # add = "reg.line"
          ) +coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") +
  theme_pubr(legend = "right") +
  scale_color_viridis_c()
################################################################################


library(vegan)
vst_trans_count_tab
adonis2(vst_trans_count_tab~Diabetes_Status,data = sam,permutations = 999)
euc_dist
betadisper(euc_dist, sam$Diabetes_Status) # 0.002
a<-adonis(euc_dist~Diabetes_Status,data = sam,permutations = 999) # 0.003
summary(a)
a$aov.tab

adonis2(euc_dist~Diabetes_Status+BMI,data = sam,permutations = 999) # 0.003
sam
