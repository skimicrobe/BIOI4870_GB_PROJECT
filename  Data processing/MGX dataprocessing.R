######################################################################
# Data processing for BIOI 4870 Project                              #
# Goal: Create a file to import Neo4j for Metagenomics(MGX) NODE     #                                    
# Author:  Suyeon Kim                                                #
# Current version: April 26, 2020                                    #
######################################################################


############################################################
# Set your working directory to the current file directory #                                              
############################################################
setwd("/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/MGX")

library(stringr)
library(tidyr)
library(reshape2)


##############################
# Input Data (Metagenomics)  #
##############################
### Read multiple input files and create a list. 
inputdat<- list.files("./", pattern = "taxonomic_profile")
mgx<-lapply(inputdat, read.table, sep="\t")
names(mgx)<- inputdat


################################################
# Subsetting 'species level' from each samples #
################################################
sampleID<- list()

for( i in names(mgx)){
  index<- mgx[[i]][,1]
  print(index)
  sampleID[[i]]<- mgx[[i]][which(str_count(index, "\\|") == 6),]
}

#############################################
# Find union of taxon across all samples    #
#############################################
taxon_ref<-c()

for ( r in 1:length(sampleID)){
  
  get_taxon<- as.character(sampleID[[r]][,1])
  taxon_ref<- c(taxon_ref, get_taxon)
  uni_taxon<- as.matrix(unique(taxon_ref))
}

###############################################################
# Curate and Create an input matrix for MGS node in Graph DB  #
###############################################################
  
microbiome_dat<-c()

for( m in 1:length(sampleID)){
  matched_sam<- sampleID[[m]][match(uni_taxon[,1], sampleID[[m]][,1]),2]
  microbiome_dat<- as.matrix(cbind(microbiome_dat, matched_sam))
}

for ( j in 1: length(names(sampleID))) {
  colnames(microbiome_dat)[j]<- names(sampleID)[[j]]
  colnames(microbiome_dat)[j]<- gsub("_taxonomic_profile.tsv","", colnames(microbiome_dat)[j])
}

bac_id<- sprintf("bacteria%s", seq(1:nrow(microbiome_dat)))
rownames(microbiome_dat)<- bac_id
microbiome_dat<- cbind(uni_taxon, microbiome_dat)
colnames(microbiome_dat)[1]<- "bacteria"

## replace NA values with zero 
microbiome_dat[is.na(microbiome_dat)]<- 0

## Reshape the data wide to long format 
mgx_mat <-melt(microbiome_dat[,2:6], id.vars = c("bacteria"))
bac_column <- as.character(microbiome_dat[,1])
mgx_mat <- cbind(mgx_mat, otu = rep(bac_column,5))

colnames(mgx_mat)[1:3]<-c("index", "sample_id", "ab_value")
add_ch<- gsub(".*\\|g__","g__", mgx_mat[,4])

final_mgx<- cbind(mgx_mat[,1:3], add_ch)
colnames(final_mgx)[4]<- "otu"

#################
# Save output   #
#################
write.csv(final_mgx, file="mgx_intput.csv",row.names=F,quote = F)



