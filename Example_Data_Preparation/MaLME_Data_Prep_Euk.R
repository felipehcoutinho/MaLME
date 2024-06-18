library(dplyr)
library(tidyr)
library(tibble)

###Process MAG abundance data
response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Delmont_et_al_2022_TARA_Euk_MAGs_Abundance.tsv"
raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

dim(raw_response_df)
colnames(raw_response_df)[1:10]

euk_response_df<-raw_response_df[,-1]
rownames(euk_response_df)<-raw_response_df[,1]
dim(euk_response_df)
rownames(euk_response_df)[1:15]

euk_response_abd_df<-euk_response_df

perc_euk_response_abd_df<-as.data.frame(t((t(euk_response_abd_df)/colSums(euk_response_abd_df))*100))
summary(colSums(perc_euk_response_abd_df))
#colnames(perc_euk_response_abd_df)[1]<-"OTU_UID"

tax_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/MaLME/Delmont_et_al_2022_TARA_Euk_MAGs_Taxonomy.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE,row.names=1)
#colnames(tax_data)[1]<-"OTU_UID"
summary(tax_data)
rownames(tax_data)[1:15]

perc_euk_response_abd_tax_df<-merge(perc_euk_response_abd_df,tax_data,by="row.names",all.x=TRUE)

perc_euk_response_abd_tax_df$Row.names[1:15]

table(perc_euk_response_abd_tax_df$Best_taxonomy_PHYLUM)

dim(perc_euk_response_abd_tax_df)
colnames(perc_euk_response_abd_tax_df[940:960])
tail(perc_euk_response_abd_tax_df[,c(940:960)])

clean_perc_euk_response_abd_tax_df<-perc_euk_response_abd_tax_df[!is.na(perc_euk_response_abd_tax_df$Best_taxonomy_PHYLUM),]
dim(clean_perc_euk_response_abd_tax_df)
colnames(clean_perc_euk_response_abd_tax_df[930:960])
tail(clean_perc_euk_response_abd_tax_df[,c(930:960)])
clean_perc_euk_response_abd_tax_df$Row.names[1:15]
colnames(clean_perc_euk_response_abd_tax_df)[1]<-"OTU_UID"

mg_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/MaLME/Delmont_et_al_2022_TARA_Metagenomes_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
summary(mg_data)

unq_fil_sizes<-as.vector(unique(mg_data$Filter_Size))

for (fsize in unq_fil_sizes) {
	#t_clean_perc_euk_response_abd_tax_df<-as.data.frame(t(clean_perc_euk_response_abd_tax_df))
	sample_subset_ids<-mg_data$Metagenome_Id[which(mg_data$Filter_Size == fsize)]
	sub_clean_perc_euk_response_abd_tax_df<-subset(clean_perc_euk_response_abd_tax_df,select=c(OTU_UID,Best_taxonomy_KINGDON,Best_taxonomy_PHYLUM,Best_taxonomy_CLASS,Best_taxonomy_ORDER,Best_taxonomy_FAMILY,Best_taxonomy_GENRE,sample_subset_ids))
	dim(sub_clean_perc_euk_response_abd_tax_df)
	outname<-paste("/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Eukaryote_Communities_Fraction_",fsize,"_Percentage_Abundance_Data.tsv",sep="")
	write.table(sub_clean_perc_euk_response_abd_tax_df,file=outname,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
}
 

#Prepare metadata
predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Delmont_et_al_2022_TARA_Sample_Metadata.tsv"

raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

colnames(raw_predictor_df)[4]<-"Temperature"
colnames(raw_predictor_df)[5]<-"Salinity"
dim(raw_predictor_df)
summary(raw_predictor_df)

mg_data<-read.table(file="/mnt/lustre/scratch/fcoutinho/MaLME/Delmont_et_al_2022_TARA_Metagenomes_Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

mg_data$Station<-as.factor(gsub("(\\d)+(\\w){4}(\\d)+$","",mg_data$Metagenome_Id,perl=TRUE))
table(mg_data$Station)

clean_predictor_df<-merge(mg_data,raw_predictor_df,by="Station",all.x=TRUE)
clean_predictor_df<-subset(clean_predictor_df,select=c(Station,Metagenome_Id,Filter_Code,Filter_Size,Latitude,Longitude,Temperature,Salinity,Si,NO3,Phos,Fe,SI_NO3,SI_T))

clean_predictor_df<-clean_predictor_df[,c("Metagenome_Id","Station","Filter_Code","Filter_Size","Latitude","Longitude","Temperature","Salinity","Si","NO3","Phos","Fe","SI_NO3","SI_T")]


colnames(clean_predictor_df)[1]<-"MG_UID"
dim(clean_predictor_df)
summary(clean_predictor_df)

write.table(clean_predictor_df,file="/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Eukaryote_Communities_Sample_Metadata.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

