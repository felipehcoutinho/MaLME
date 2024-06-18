library(dplyr)
library(tidyr)
library(tibble)


#Prepare metadata
predictor_file<-"/mnt/lustre/scratch/fcoutinho/StG/Edited_TARA_GOV2_Viromes_Metadata.tsv"

raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

colnames(raw_predictor_df)<-c("Sample","Station","Latitude","Longitude","Depth","Biome","Region","Province","Temperature","Nitrogen","Ammonium_5m","Nitrate","Nitrate_5m","Nitrite","Nitrite_5m","Nitrate_Nitrite","Phosphate","Silicate","Iron_5m","Oxygen","Salinity","Sigma_Theta","Chlorophyll_A","POC","PIC","Residence_Time","Conductivity","PAR_1_1","PAR_1_8","PAR_1_30","Season","Season_Stage","Distance_From_Coast","Angular_Scattering_Coefficient","Optical_Backscattering_Coefficient","Backscattering_Coefficient_of_Particles","Fluorescence","Optical_Beam_Attenuation_Coefficient","Beam_Attenuation_Coefficient_of_Particles","Sunshine_Duration","Sea_Ice_Concentration","Sea_Ice_Free_Period.","Sample_Accession_Number","MG_UID","Ecological_Zone")
	
raw_predictor_df<-raw_predictor_df[,c("MG_UID","Sample","Station","Latitude","Longitude","Depth","Biome","Region","Province","Temperature","Nitrogen","Ammonium_5m","Nitrate","Nitrate_5m","Nitrite","Nitrite_5m","Nitrate_Nitrite","Phosphate","Silicate","Iron_5m","Oxygen","Salinity","Sigma_Theta","Chlorophyll_A","POC","PIC","Residence_Time","Conductivity","PAR_1_1","PAR_1_8","PAR_1_30","Season","Season_Stage","Distance_From_Coast","Angular_Scattering_Coefficient","Optical_Backscattering_Coefficient","Backscattering_Coefficient_of_Particles","Fluorescence","Optical_Beam_Attenuation_Coefficient","Beam_Attenuation_Coefficient_of_Particles","Sunshine_Duration","Sea_Ice_Concentration","Sea_Ice_Free_Period.","Sample_Accession_Number","Ecological_Zone")]


dim(raw_predictor_df)
summary(raw_predictor_df)


write.table(raw_predictor_df,file="/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Viral_Communities_Sample_Metadata.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


###Process OTU RPKM abundance data
response_file<-"/mnt/lustre/scratch/fcoutinho/StG/Corrected_Abundance/Min_Count_1k_RPKM_Abundance_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq.tsv"
#response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/short.tsv"
raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

#t_raw_response_df<-as.data.frame(t(raw_response_df))
#colnames(t_raw_response_df)<-rownames(raw_response_df)

MG_UIDs<-colnames(raw_response_df)
raw_response_df$Sequence<-rownames(raw_response_df)

dim(raw_response_df)

colnames(raw_response_df)

response_info_file<-"/mnt/lustre/scratch/fcoutinho/StG/Min_Count_1k_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_Full_Seq_Info.tsv"
response_info_df<-read.table(file=response_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	
full_response_info_df<-merge(raw_response_df,subset(response_info_df,select=c("Sequence","VPF_family","VPF_baltimore","VPF_genus","Predicted_Host_1_Domain","Predicted_Host_2_Phylum","Predicted_Host_3_Class","Predicted_Host_4_Order","Predicted_Host_5_Family","Predicted_Host_6_Genus")),all.x=TRUE,by="Sequence")

full_response_info_df<-full_response_info_df[,c("Sequence","VPF_family","VPF_baltimore","VPF_genus","Predicted_Host_1_Domain","Predicted_Host_2_Phylum","Predicted_Host_3_Class","Predicted_Host_4_Order","Predicted_Host_5_Family","Predicted_Host_6_Genus",MG_UIDs)]

colnames(full_response_info_df)

summary(full_response_info_df)

#colnames(full_response_info_df)[1]<-"OTU_UID"

write.table(full_response_info_df,file="/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Virus_Communities_RPKM_Abundance_Data.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

