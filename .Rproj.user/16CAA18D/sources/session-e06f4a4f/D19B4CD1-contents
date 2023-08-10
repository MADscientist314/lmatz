library(pheatmap)
library(tidyverse)

m<-metacyc_abundance%>%
  mutate(pathway=rownames(metacyc_abundance))%>%
           filter(pathway%in%rownames(res3))%>%
  select(-pathway)
metadata
m

key2<-key%>%filter(pathway%in%rownames(m))
key2
m
annot_col<-data.frame(metadata)
rownames(annot_col)<-metadata$SampleID  
annot_col<-annot_col%>%select(-c(SampleID,Age,BMI,Sex))

dim(res2)
key2
res2
res3<-res2%>%filter(log2FoldChange>2.5)
# res3<-res2%>%filter(log2FoldChange>2.5)

dim(res3)

annot_col
 # annot_col<-annot_col %>% mutate(Diabetes_Status=gsub("T2D","sT2D",Diabetes_Status))
color_pal<-list(Diabetes_Status=c("sT2D"="#FF7100","HC"="#240A94"))

m



pheatmap(m,
         #clustering_distance_rows = "manhattan",
        #clustering_distance_cols = "correlation",
         
        cutree_rows = 3,
        annotation_colors = color_pal,
         cutree_cols = 3,
          scale = "row",
         border_color = "black",cellheight = 10,cellwidth = 10,
         annotation_col =   annot_col,
         labels_row = key2$description)
dev.off()

