library(dplyr)
library(tidyr)
library(tibble)

###Process OTU raw abundance data
response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/mitags_tab_otu.tsv"
#response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/short.tsv"
raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

dim(raw_response_df)
#rownames(raw_response_df)<-raw_response_df[,1]
t_raw_response_df<-as.data.frame(t(raw_response_df))
#colnames(t_raw_response_df)[1]<-"Full_Taxonomy"
#summary(t_raw_response_df$Full_Taxonomy)
rownames(t_raw_response_df)[26:100]

t_raw_response_df$Full_Taxonomy<-rownames(t_raw_response_df)
t_raw_response_df<-t_raw_response_df %>% separate(Full_Taxonomy, c("Domain","Phylum","Class","Order","Family","Genus","Species"),sep=";")
t_raw_response_df<-t_raw_response_df %>% separate(Domain, c("OTU_ID","Domain"),sep=" ")
dim(t_raw_response_df)
colnames(t_raw_response_df)[181:188]
#summary(t_raw_response_df[,181:187])
head(t_raw_response_df[,181:188])
#"OTU_ID",
table(t_raw_response_df$Domain)
table(t_raw_response_df$Phylum)
table(t_raw_response_df$Class)
table(t_raw_response_df$Order)

#Prokaryotes
prok_response_df<-t_raw_response_df[which((t_raw_response_df$Domain == "Bacteria") | t_raw_response_df$Domain == "Archaea"),]

dim(prok_response_df)
table(prok_response_df$Domain)
table(prok_response_df$Phylum)
table(prok_response_df$Class)

prok_response_abd_df<-subset(prok_response_df,select=-c(OTU_ID,Domain,Phylum,Class,Order,Family,Genus,Species))
dim(prok_response_abd_df)

prok_response_tax_df<-subset(prok_response_df,select=c(OTU_ID,Domain,Phylum,Class,Order,Family,Genus,Species))
dim(prok_response_tax_df)
head(prok_response_tax_df)

perc_prok_response_abd_df<-as.data.frame(t((t(prok_response_abd_df)/colSums(prok_response_abd_df))*100))
summary(colSums(perc_prok_response_abd_df))

perc_prok_response_full_df<-merge(prok_response_tax_df,perc_prok_response_abd_df,by="row.names",all.x=TRUE)
dim(perc_prok_response_full_df)

colnames(perc_prok_response_full_df)[1]<-"Full_OTU_ID"

colnames(perc_prok_response_full_df)

write.table(perc_prok_response_full_df[,-1],file="/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Prokaryote_Communities_Percentage_Abundance_Data.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#Prepare metadata
predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Salazar_et_al_2019_Table_S4_TARA_Metadata.tsv"

raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)

dim(raw_predictor_df)
summary(raw_predictor_df)

colnames(raw_predictor_df)[1]<-"MG_UID"

write.table(raw_predictor_df,file="/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Prokaryote_Communities_Sample_Metadata.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

