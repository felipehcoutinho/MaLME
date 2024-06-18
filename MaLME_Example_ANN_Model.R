#Read in the example dataset table as a Data Frame
example_df<-read.table(file="/mnt/smart/users/fcoutinho/Repos/MaLME/Example_Datasets/Example_ANN_Input_Data_2_Prok.tsv",sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1)

#Print a summary of the dataframe
summary(example_df)

#Subset the table to only keep complete cases (i.e. Samples with non  NA values for all predictor and response variables)
sub_example_df<-example_df[complete.cases(example_df), ]

summary(sub_example_df)

#Load the ANN training library
library(nnet)

#Load the caret library
library(caret)

#Rename response variable
colnames(sub_example_df)[which(colnames(sub_example_df) == "Shannon_Diversity_Index")]<-"Response_Variable"

#Z-transform the data
sub_example_df<-as.data.frame(scale(sub_example_df,center=TRUE,scale=TRUE))

#Set a seed number to make your results reproducible
#seed_num<-12223
seed_num<-666
set.seed(seed_num)

#Assign samples to the test and validation sets
train_set_row_nums<-createDataPartition(sub_example_df$Response_Variable,p=0.5,list=FALSE)

sub_example_df_train<-sub_example_df[train_set_row_nums,]

sub_example_df_valid<-sub_example_df[-train_set_row_nums,]

dim(sub_example_df_train)
dim(sub_example_df_valid)

summary(sub_example_df_train)
summary(sub_example_df_valid)

#Test if the differences in the values of the response variable are significantly different between train and test sets
wilcox.test(sub_example_df_train$Response_Variable,sub_example_df_valid$Response_Variable,paired=FALSE,exact=FALSE)

#Train a model using the training set only
set.seed(seed_num)

trained_net_ts<-nnet(Response_Variable ~ Temperature+Oxygen+Depth.nominal+ChlorophyllA+Salinity+Ammonium.5m+Iron.5m,sub_example_df_train,size=3,linout=TRUE,maxit=10)

#Check the contents of the output ANN object
ls(trained_net_ts)

#Check output values
summary(trained_net_ts$fitted.values)

#Load the Metrics library to evaluate model performance
library(Metrics)

#Evaluate the performance of the model on the training set
postResample(sub_example_df_train$Response_Variable,trained_net_ts$fitted.values)

#Get predictions for the validation set
valid_preds<-predict(trained_net_ts,sub_example_df_valid,type=c("raw"))

summary(valid_preds)

#Evaluate the performance of the model on the validation set
postResample(valid_preds,sub_example_df_valid$Response_Variable)

set.seed(seed_num)

#Train the network on full set
trained_net<-nnet(Response_Variable ~ Temperature+Oxygen+Depth.nominal+ChlorophyllA+Salinity+Ammonium.5m+Iron.5m,sub_example_df,size=3,linout=TRUE,maxit=10)

#Evaluate the performance of the new model on the full set
postResample(trained_net$fitted.values, sub_example_df$Response_Variable)

#Load the ANN interpretation library
library(NeuralNetTools)

#Calculate the importance of the predictors using the Olden method
importance_olden<-olden(trained_net,bar_plot=FALSE)

#Look at the outpt and identify the most important predictor
importance_olden
