###########################################################################
# Data processing for BIOI 4870 Project                                   #
# Goal: Create a file to import Neo4j for Proteomics' and 'Metabolomics'  #                                    
# Author:  Suyeon Kim                                                     #
# Current version: April 27, 2020                                         #
###########################################################################

setwd("/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/")

library(reshape2)

################
# List of Data #
################
proteomics<-read.csv("/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/HMP2_proteomics_ecs.tsv", header=T, sep="\t")
metabolomics<-read.csv("/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/HMP2_metabolomics.csv", header=T, sep=",")

##################################################### 
# Reshape the data wide to long format (proteomics) #
#####################################################
#px_mat <-melt(proteomics[,2:ncol(proteomics)], id.vars = c("Gene"))      
px_mat<- melt(proteomics, id.vars= c("Gene"))
filt_px<- gsub(".*\\:", "", px_mat$Gene)
filt_mat<- cbind(filt_px, px_mat[,2:ncol(px_mat)])
colnames(filt_mat)[2:3]<-c("sample_id", "ab_value")


####################################################### 
# Reshape the data wide to long format (metabolomics) #
#######################################################
metabolites<- metabolomics[,6]
mb_mat<- cbind(metabolites, metabolomics[,8:ncol(metabolomics)])
mb_reform <-melt(mb_mat, id.vars = c("metabolites"))   
colnames(mb_reform)[2:3]<- c("sample_id","ab_value")

#################
# Save outputs  #
#################
write.csv(filt_mat, file="proteomics_intput.csv",row.names=F, quote=F)
write.csv(mb_reform, file="metabolite_input.csv",row.names=F, quote=F)


######################################################################################################################
# Due to limtied space for ODIN server, One sample (CSM5MCVN) from Protein and metabolites have used for this project.#
#######################################################################################################################
matched_px<- px[grep("CSM5MCVN", px$sample_id),]
colnames(matched_px)[1]<- "protein"
filt_p<- as.matrix(matched_px[which(matched_px$ab_value !=0),])

matched_mb<- mb[grep("CSM5MCVN", mb$sample_id),]
matched_mb$ab_value<- as.numeric(as.character(matched_mb$ab_value))
matched_mb[is.na(matched_mb)]<-0
filt_mb<- as.matrix(matched_mb[which(matched_mb$ab_value !=0),])

rm_mb_old<-filt_mb[!(!is.na(filt_mb$metabolites) & filt_mb$metabolites == ""),]
rm_mb<- filt_mb[!(is.na(filt_mb$metabolites) | filt_mb$metabolites== ""),]

##################################################
# Save outputs (filtered protein and metablites) #
##################################################
write.csv(filt_p, file="/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/matched_px_input.csv",row.names=F,quote = F)
write.csv(rm_mb, file="/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/matched_mb_input.csv", row.names=F, quote=F)






