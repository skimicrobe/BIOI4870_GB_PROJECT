######################################################################
# Data processing for BIOI 4870 Project                              #
# Goal: Create a file to import Neo4j for Metagenomics(MTX) NODE     #                                    
# Author:  Suyeon Kim                                                #
# Current version: April 26, 2020                                    #
######################################################################

######################################################################
# Set your working directory to the current file directory                                            #
######################################################################
setwd("/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/MTX")

library(stringr)
library(dplyr)

##############################################################################################
# Data processing : Reformat Data in which include 'geneID with species' alone with sampleID #
##############################################################################################
geneFam_processing<- function(mtx) {
  
  filt_gene<- list()
  
  for( i in names(mtx)) {
    index<- mtx[[i]][,1]
    print(index)
    filt_gene[[i]]<- mtx[[i]][grep("g__", index),]
  }
  
  
  ### Filter information 
  for( i in names(filt_gene)) {
    index_info<- filt_gene[[i]][,1]
    #print(index_info)
    filt_gene[[i]][,1]<- gsub("UniRef90_", "", index_info)
  }
  
  ### Separate geneID and species name
  new_mtx<- list()
  for ( j in names(filt_gene)) {
    split_col<- str_split_fixed(as.character(filt_gene[[j]][,1]),"\\|", n=2 )
    new_mtx[[j]]<- cbind(split_col, filt_gene[[j]][,2])  
  }
  
  ### Add column for smapleIDs
  filt_mtx<- list()
  
  for (r in 1:length(names(new_mtx))){ 
    sample_name<-as.character(names(new_mtx)[r])
    filt_mtx[[r]]<- cbind(new_mtx[[r]], sample_id= rep(sample_name, nrow(new_mtx[[r]])))
    colnames(filt_mtx[[r]])<- c("gene", "bacteria", "ab_value","sample_id") 
  }
  return(filt_mtx)
}

pathWay_processing<- function(path_data) {
  filt_path<- list()
  
  for( i in names(path_data)) {
    index<- path_data[[i]][,1]
    print(index)
    filt_path[[i]]<- path_data[[i]][grep("g__", index),]
  }
  
  
  test<-### Filter information 
  for( i in names(filt_path)) {
    index_info<- filt_path[[i]][,1]
    #print(index_info)
    filt_path[[i]][,1]<- gsub(".*\\:", "", index_info)
  }
  
  ### Separate geneID and species name
  new_path<- list()
  for ( j in names(filt_path)) {
    split_col<- str_split_fixed(as.character(filt_path[[j]][,1]),"\\|", n=2 )
    new_path[[j]]<- cbind(split_col, filt_path[[j]][,2])  
  }
  
  ### Add column for smapleIDs
  filt_path_final<- list()
 
  for (r in 1:length(names(new_path))){ 
    sample_name<-as.character(names(new_path)[r])
    filt_path_final[[r]]<- cbind(new_path[[r]], sample_id= rep(sample_name, nrow(new_path[[r]])))
    colnames(filt_path_final[[r]])<- c("pathway", "bacteria", "ab_value","sample_id") 
  
  }
  return(filt_path_final)
}

######################################
# List of Input Data (Metagenomics)  #
######################################

### Gene Family 
inputdat<- list.files("./", pattern = "genefamilies")
mtx<-lapply(inputdat, read.table, sep="\t")
names(mtx)<- inputdat

### Pathway 
inputdat<- list.files("./", pattern = "pathabundance")
path_data<-lapply(inputdat, read.csv, sep="\t")      ### read.table 
names(path_data)<- inputdat

###################
#  List of Output #
###################

###### Mtx - GENE ##############
mtx_gene<- geneFam_processing(mtx)
#test_mtx<- unlist(lapply(mtx_gene, function(x) rbind.data.frame(c(t(x)))))
mtx_mat<- as.data.frame(do.call(rbind, mtx_gene))

## Convert List to Data Frame Rows
mtx_mat[,2]<-gsub("\\.","|", mtx_mat[,2])
mtx_mat[,4]<- gsub("_genefamilies.tsv", "", mtx_mat[,4])

#########################
# Save output for gene  #
#########################
write.csv(mtx_mat, file="/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/MTX/mtx_gene_input.csv",row.names = F, quote = F)

##### MTX - Pathway #############
mtx_pathway<- pathWay_processing(path_data)

mtx_path_mat<- as.data.frame(do.call(rbind, mtx_pathway))
mtx_path_mat[,2]<- gsub("\\.","|", mtx_path_mat[,2])
mtx_path_mat[,4]<- gsub("_pathabundance.tsv", "", mtx_path_mat[,4])

###########################
# Save output for pathway #
###########################
write.csv(mtx_path_mat, file="/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/MTX/mtx_pathway_input.csv",row.names = F, quote = F)






