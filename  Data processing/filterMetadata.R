######################################################################
# Goal: Filtering the metadata                                       #
# - Subset useful variables related to patients info                 #
# Author:  Suyeon Kim                                                #
######################################################################
setwd("/Users/suyeonkim/Box Sync/BIOI4870 Project - Data/")
######################################################################
# SHMP2 Metadata                                                     #
######################################################################
metadat<- read.csv("/Users/suyeonkim/Google Drive/YEAR_OF_2020/Class_2020_spring/BIOI4870_db_search_pattern_dis_bioinfo/Project/data/hmp2_metadata.csv", sep=",")

######################################################################
# Subset csv file with useful clinical info.                         #
######################################################################

patient_info<- metadat[,c("External.ID", "data_type", "diagnosis", "Rectum.Flora." , "Ileum.flora.", "sex", "Antibiotics","biopsy_location")]

#uniq_sample<- patient_info[unique(patient_info$External.ID),]
filt_mat<- patient_info[,c("External.ID", "diagnosis", "sex", "Antibiotics")]
colnames(filt_mat)<-c("External.ID", "Diagnosis_info", "Gender", "Antibiotics")

unique_samples <- filt_mat[-which(duplicated(filt_mat)),]

########################
# Save output file     #
########################

write.csv(unique_samples, file="filt_metadata.csv",row.names=F,quote = F)
