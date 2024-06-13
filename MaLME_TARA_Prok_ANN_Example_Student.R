full_df<-read.table(file="TARA_Marine_Prokaryote_Communities_Filtered_OTU_Abundance+Metadata.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

selected_OTU_ID<-"AACY020552977.734.2233"
full_df$Response_Var<-full_df[[selected_OTU_ID]]

library(nnet)
library(Metrics)

set.seed(666)

net_full<-nnet(Response_Var ~ Temperature+Salinity+Ammonium.5m+Iron.5m+ChlorophyllA,full_df,size=5,linout=TRUE,maxit=1000)

#Calculate Pearson R between measured and predicted values of the model using the training set
cor.test(full_df$Response_Var,net_full$fitted.values, method="pearson")

#Evaluate the performance of the model on the training set
rmse(full_df$Response_Var,net_full$fitted.values)

#Extrapolate for a +2 Degrees warming scenario
full_df_ext<-full_df
full_df_ext$Temperature<-full_df_ext$Temperature+2

ext_preds<-predict(net_full,full_df_ext,type=c("raw"))

summary(net_full$fitted.values)
summary(ext_preds)

wilcox.test(net_full$fitted.values,ext_preds,paired=FALSE,exact=FALSE)

