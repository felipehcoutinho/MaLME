#Read abundance table  (samples are rows and genomes and environmental parameters are columns)
abd_df<-read.table(file="/home/felipe.coutinho/MaLME/Example_Datasets/Marine_Prokaryote_Communities+MetaData.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1)

summary(abd_df[,c(1:15)])

#Read taxonomy table (genomes are rows and taxonomic classification are columns)
tax_df<-read.table(file="/home/felipe.coutinho/MaLME/Example_Datasets/Marine_Prokaryote_Communities_Taxonomy_Data.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1)

summary(tax_df)

#Check if order of genomes match in the tables
summary(rownames(tax_df) == colnames(abd_df)[12:ncol(abd_df)])

#Select a taxon subset
sub_abd_df<-abd_df[,rownames(tax_df)[which(tax_df$Phylum=="Cyanobacteria")]]

#Select a taxon subset and environmental parameters
sub_abd_df<-abd_df[,c("Temperature","Salinity",rownames(tax_df)[which(tax_df$Phylum=="Cyanobacteria")])]