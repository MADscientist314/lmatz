#############################################################################################
# Trimmomatic.sh
#############################################################################################

for n in $(cat manifest.tsv)
do 
echo "running trimmomatic on sample $n";
java -jar /media/jochum00/Aagaard_Raid1/reference_datasets/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 48 \
/home/oliviastuyck/lmatz/LM_RawSequences/fastq/$n.1.fq \
/home/oliviastuyck/lmatz/LM_RawSequences/fastq/$n.2.fq \
/home/oliviastuyck/lmatz/LM_RawSequences/fastq/$n.1.paired.fq \
/home/oliviastuyck/lmatz/LM_RawSequences/fastq/$n.1.unpaired.fq \
/home/oliviastuyck/lmatz/LM_RawSequences/fastq/$n.2.paired.fq \
/home/oliviastuyck/lmatz/LM_RawSequences/fastq/$n.2.unpaired.fq \
ILLUMINACLIP:/media/jochum00/Aagaard_Raid1/reference_datasets/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36;
done;
#############################################################################################

#############################################################################################
# picrust2.sh
#############################################################################################

###############################################################
#Place amplicon sequence variants (or OTUs) into reference phylogeny
place_seqs.py \
-s ASVs.fa \
-o placed_seqs.tre \
-p 48 \
--intermediate placement_working

###############################################################
#Run hidden-state prediction to get 16S copy numbers, E.C. number, and KO abundances per predicted genome
hsp.py -i 16S -t placed_seqs.tre -o marker_nsti_predicted.tsv.gz -p 48 -n
# EC numbers
hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 48
# KEGG orthologs
hsp.py -i KO -t placed_seqs.tre -o KO_predicted.tsv.gz -p 48

###############################################################
# Predict E.C. and KO abundances in sequencing samples 
#(adjusts gene family abundances by 16S sequence abundance) 
# be sure to include stratified contribution output

metagenome_pipeline.py -i ASVs_counts.tsv \
                       -m marker_nsti_predicted.tsv.gz \
					   --strat_out \
                       -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out

metagenome_pipeline.py -i ASVs_counts.tsv \
						--strat_out \
                       -m marker_nsti_predicted.tsv.gz \
                       -f KO_predicted.tsv.gz \
                       -o KO_metagenome_out

###########################################################################					
#Infer MetaCyc pathway abundances and coverages based on predicted E.C. number abundances 
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -o pathways_out \
                    --intermediate pathways_working \
                    -p 48

###########################################################################					
#Add descriptions as new column in gene family and pathway abundance tables 

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
					-m EC \
                    -o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
					-m KO \
                    -o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o pathways_out/path_abun_unstrat_descrip.tsv.gz
					