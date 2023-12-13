library(dada2)
library(tidyverse)
setwd("/home/jochum00/lmatz/LM_RawSequences/paired")

######################################
### DATA WRANGLING
#import the dataframe with the manifest,tsv and make vectors of the fwd and rev
#######################################

df<-read.table("../manifest.tsv",header = F,sep = "\t",row.names = NULL)
colnames(df)<-"samples"
df<-df%>%mutate(fwd=paste0(samples,".1.paired.fq"),
                rev=paste0(samples,".2.paired.fq"),
                filt_fwd=paste0(samples,".1.paired.filt.fq"),
                filt_rev=paste0(samples,".2.paired.filt.fq"),
                )

# one holding the file names of all the forward reads
forward_reads <- df$fwd
# and one with the reverse
reverse_reads <- df$rev

# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <-df$filt_fwd
filtered_reverse_reads <- df$filt_rev
################################################################################
# DATA FILTERING
################################################################################

# and just plotting the last 4 samples of the reverse reads
plotQualityProfile(reverse_reads[1:5])
plotQualityProfile(forward_reads[1:5])

################################################################################
filtered_out <- filterAndTrim(forward_reads, 
                              filtered_forward_reads,
                              reverse_reads, 
                              filtered_reverse_reads, 
                              maxEE=c(2,2),
                              rm.phix=TRUE, 
                              multithread = T,
                              minLen=50,
                              truncLen=c(250,200))
################################################################################
class(filtered_out) # matrix
dim(filtered_out) # 20 2

filtered_out
#           
# reads.in reads.out
# Buffington_201_001.1.paired.fq    17878     15287
# Buffington_201_002.1.paired.fq    19088     17182
# Buffington_201_003.1.paired.fq    18775     16559
# Buffington_201_005.1.paired.fq    22986     19944
# Buffington_201_006.1.paired.fq    19611     17239
# Buffington_201_007.1.paired.fq    23630     20788
# Buffington_201_008.1.paired.fq    19337     17473
# Buffington_201_009.1.paired.fq    13076     11231
# Buffington_201_010.1.paired.fq    22531     19794
# Buffington_201_011.1.paired.fq    19323     15999
# Buffington_201_013.1.paired.fq    17213     15100
# Buffington_201_014.1.paired.fq    21178     17680

################################################################################
# learn the error rates
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE) 
# 51069000 total bases in 204276 reads from 12 samples will be used for learning the error rates.
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE) 
# 40855200 total bases in 204276 reads from 12 samples will be used for learning the error rates.
# plot the learned error rates

plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
################################################################################
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <-df$samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- df$samples
##################################################################################
#run DADA2 
dada_forward <- dada(derep_forward, err=err_forward_reads, pool=T, multithread=TRUE)
# 12 samples were pooled: 204276 reads in 46370 unique sequences.
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool=T, multithread=TRUE)
# 12 samples were pooled: 204276 reads in 44362 unique sequences.
##################################################################################
#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE, minOverlap=150)

# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
##################################################################################
# Generate a count table
##################################################################################
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) #  12 729
##################################################################################
#Chimera identification
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) 
# Identified 195 bimeras out of 729 input sequences.
# though we only lost 17 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) 
#  0.9485947
################################################################################
#Overview of counts throughout
# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=df$samples, 
                          dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], 
                          dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN),
                          merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab

write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)
################################################################################
################################################################################
################################################################################
#Assigning taxonomy
################################################################################
################################################################################
################################################################################
## loading DECIPHER
library(DECIPHER)
# download the reference taxonomy data
## downloading DECIPHER-formatted SILVA v138 reference
# download.file(url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData", destfile="SILVA_SSU_r138_2019.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2019.RData")
## creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
## and classifying
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=48)
################################################################################
#Extracting the standard goods from DADA2
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

