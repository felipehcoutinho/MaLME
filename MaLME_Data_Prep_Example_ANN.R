#Process metadata
predictor_file<-"/mnt/lustre/scratch/fcoutinho/StG/Edited_TARA_GOV2_Viromes_Metadata.tsv"

raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(raw_predictor_df)<-c("Sample","Station","Latitude","Longitude","Depth","Biome","Region","Province","Temperature","Nitrogen","Ammonium_5m","Nitrate","Nitrate_5m","Nitrite","Nitrite_5m","Nitrate_Nitrite","Phosphate","Silicate","Iron_5m","Oxygen","Salinity","Sigma_Theta","Chlorophyll_A","POC","PIC","Residence_Time","Conductivity","PAR_1_1","PAR_1_8","PAR_1_30","Season","Season_Stage","Distance_From_Coast","Angular_Scattering_Coefficient","Optical_Backscattering_Coefficient","Backscattering_Coefficient_of_Particles","Fluorescence","Optical_Beam_Attenuation_Coefficient","Beam_Attenuation_Coefficient_of_Particles","Sunshine_Duration","Sea_Ice_Concentration","Sea_Ice_Free_Period.","Sample_Accession_Number","Run_Accession_Number","Ecological_Zone")
	

predictor_variables<-c("Temperature","Depth","Oxygen","Chlorophyll_A","Salinity","Iron_5m","Ammonium_5m","Run_Accession_Number")
	
raw_predictor_df<-subset(raw_predictor_df,select=predictor_variables)

dim(raw_predictor_df)
summary(raw_predictor_df)
###Process virome RPKM data
response_file<-"/mnt/lustre/scratch/fcoutinho/StG/Corrected_Abundance/RPKM_Abundance_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_rep_seq.tsv"

raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1)

dim(raw_response_df)
colnames(raw_response_df)

library(vegan)

shannon_diversity<-diversity(t(raw_response_df[,-1]), index = "shannon")

length(shannon_diversity)
summary(shannon_diversity)

shannon_diversity_df<-cbind(names(shannon_diversity),shannon_diversity)
colnames(shannon_diversity_df)<-c("Run_Accession_Number","Shannon_Diversity")
shannon_diversity_df<-as.data.frame(shannon_diversity_df)
shannon_diversity_df$Shannon_Diversity<-as.numeric(shannon_diversity_df$Shannon_Diversity)
dim(shannon_diversity_df)
summary(shannon_diversity_df)

#Merge and print
full_df<-merge(raw_predictor_df,shannon_diversity_df,by="Run_Accession_Number",all.x=TRUE)

dim(full_df)
summary(full_df)

write.table(full_df,file="/mnt/lustre/scratch/fcoutinho/MaLME/Example_Data_1.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)