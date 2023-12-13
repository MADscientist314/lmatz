library(phyloseq)
library(vegan)
library(microbiome)
library(ggpubr)
library(ggsci)
library(mosaic)

counts<-read.table("ASVs_counts.tsv",sep = "\t",header = T,row.names = 1)
tax<-read.table("ASVs_taxonomy.tsv",sep = "\t",header = T,row.names = 1,na.strings = "NA")
sam<-data.frame(metadata)
rownames(sam)<-sam$SampleID
colnames(counts)
tax
ps<-readRDS("ps.RDS")
# ps<-phyloseq(otu_table(counts,taxa_are_rows = T),tax_table(as.matrix(tax)),sample_data(sam))
ps

#################################################################################

#Data filtering
ps_genus<-tax_glom(physeq = ps,taxrank = "genus")
ps_genus_core<-core(ps_genus,detection = 1/100,prevalence = 1/100)
ps_genus_core_rel<-transform(x = ps_genus_core,"compositional")
################################################################################
top_tax<-top_taxa(n = 20,x = ps_genus_core_rel)

melt<-psmelt(ps_genus_core_rel)
keep<-melt%>%filter(OTU%in%top_tax)
other<-melt%>%filter(!OTU%in%top_tax)%>%mutate(genus="Other")
ps2<-full_join(keep,other)
sam$SampleID
ps2$SampleID
ps2<-ps2$>$arrange()
ggbarplot(data = ps2,
          x = "SampleID",
          y = "Abundance",
          color = "genus",
          fill = "genus")+
  scale_fill_d3("category20")+
  rotate_x_text(45)
################################################################################
sample_sums(ps)
sam$Diabetes_Status

write.table(x = estimate_richness(physeq = ps),file = "alpha_diversity_metrics.tsv",sep = "\t",row.names = F)

plot_richness(ps, 
              x="Diabetes_Status", 
              color="Diabetes_Status", 
              measures=c("Chao1", "Shannon","InvSimpson")) +
  theme_bw() + 
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


