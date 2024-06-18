#Read in predictor/response table
predictor_file<-"/mnt/smart/users/fcoutinho/Repos/MaLME/Example_Datasets/Example_RF_Input_Data_2_Prok_Subset.tsv"

raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

#Rename response variable
colnames(raw_predictor_df)[which(colnames(raw_predictor_df) == "Is_Polar")]<-"Response_Variable"

summary(raw_predictor_df[,1:10])

#Subset the table to only keep complete cases (i.e. Samples with non  NA values for all predictor and response variables)
raw_predictor_df<-raw_predictor_df[complete.cases(raw_predictor_df), ]

dim(raw_predictor_df)

#Z-transform the data
library(scales)

raw_predictor_df[,c(2:ncol(raw_predictor_df))]<-as.data.frame(scale(raw_predictor_df[,-1],center=TRUE,scale=TRUE))

#Load the caret library
library(caret)

seed_num<-666
set.seed(seed_num)

#Assign samples to the training and validation sets
train_set_row_nums<-createDataPartition(raw_predictor_df$Response_Variable,p=0.5,list=FALSE)

raw_predictor_df_train<-raw_predictor_df[train_set_row_nums,]

raw_predictor_df_valid<-raw_predictor_df[-train_set_row_nums,]

table(raw_predictor_df_train$Response_Variable)
table(raw_predictor_df_valid$Response_Variable)

#Train a model using the training set only
library(ranger)

set.seed(seed_num)
rf_model<-ranger(data=raw_predictor_df_train, formula = Response_Variable ~ . , num.trees=100, write.forest=TRUE, num.threads=1, save.memory =TRUE)

#Evaluate the performance of the model on the training set
confusionMatrix(as.factor((raw_predictor_df_train$Response_Variable+0)),as.factor(rf_model$predictions),mode ="prec_recall",positive="1")

#Get predictions for the validation set
valid_preds<-predict(rf_model,raw_predictor_df_valid)

#Evaluate the performance of the model on the validation set
confusionMatrix(as.factor((raw_predictor_df_valid$Response_Variable+0)),as.factor(valid_preds$predictions),mode ="prec_recall",positive="1")

#Fit single model with all samples
set.seed(seed_num)
rf_model<-ranger(data=raw_predictor_df, formula = Response_Variable ~ . , num.trees=100, write.forest=TRUE,importance="impurity", num.threads=1, save.memory =TRUE)

#Evaluate the performance of the model on the full set model
confusionMatrix(as.factor((raw_predictor_df$Response_Variable+0)),as.factor(rf_model$predictions),mode ="prec_recall",positive="1")

#Evaluate predictor importance  on the full set model
imp_df<-as.data.frame(cbind(names(rf_model$variable.importance),rf_model$variable.importance))
colnames(imp_df)<-c("Predictor","Importance")
imp_df$Importance<-as.numeric(imp_df$Importance)

#Merge the importance data with the available predictor information (e.g. Taxonomy, Metabolic traits, ecological niche, etc)
predictor_info_df<-read.table(file="/mnt/smart/users/fcoutinho/Repos/MaLME/Example_Datasets/Marine_Prokaryote_Communities_Taxonomy_Data.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

imp_df<-merge(imp_df,predictor_info_df,by.x="Predictor",by.y="OTU_ID",all.x=TRUE)

imp_df$Predictor<-as.factor(imp_df$Predictor)

#Look at the top 20 most important predictors
imp_df[order(imp_df$Importance,decreasing=TRUE)[1:20],]

#Plot predictor importance according to taxonomy
tax_imp_box<-ggplot(imp_df[which(imp_df$Importance > 0),])+
geom_boxplot(aes(fill=Phylum, y=log2(Importance),x=Phylum))+
theme_bw()+
theme(legend.position="top",axis.text.x = element_text(angle = 45,hjust = 1))

ggsave("RF_Example_ImportancexPhylum_Boxplot.pdf",plot=tax_imp_box,width=7,height=5,pointsize=8)
