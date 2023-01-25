#This code builds random forest models for the predictions of hosts of prokaryote viruses based on a given set of predictors.
#In this case the predictors are computed hmmsearch scores of genomic proteins versus a custom database of HMM models 
#It should also work with any other set of numerical predictors (e.g. number of occurrence of pFAM od pVOG domains in genomes)

library(ranger)
library(PRROC)
options(expressions = 5e5)

tree_num<-1000 #Number of trees to grow
var_tries<-5000 #Variables to try at each split when growing trees
valid_filter_votes<-0.6 #Minimum percentage to consider a prediction of the RF as valid
min_length<-1000 #Mininum genome length to be included when training and validating models
min_peg_count<-1 #Minimum number of protein encoding genes a genome must have to be used in training and validation steps
train_frac<-0.8 #Fraction of data to be used for training when building Model1 (Trained on Train_frac fraction of RefSeq)
min_dom_prevalence<-1 #Minimum Number of genomes in which predictor has to be have value above 0 to be included in the models
filter_correl_cutoff<-0.9 #Correlation cutoff to be passed to findCorrelation function to identify which predictors to exclude
threads<-120 #Number of threads to use when growing trees and performing predictions
kfolds<-10 #Number of k-fold cross valiation to be done using model 1
node_size<-1 #Node size for the probabilistic forest
prefix<-"New" #Prefix to add to output files 


#in_file is the input matrix that is used to train and validate models
#in_file MUST be a .tsv matrix of Genomes (Rows) x Predictors (Columns)
#The structure of in_file columns , MUST be, in this order: Host_Taxonomy_6_Genus,Cluster_Representative_RefSeq,Length,PEG_Count,Original_File,Recruitment_Criteria_Updated
#Where Host_Taxonomy_6_Genus is the name of the host of the viral genome
#Cluster_Representative_RefSeq is a TRUE/FALSE/NA column indicating is the genome is the representative of its group of non-redundant genomes (Only those with a value of TRUE will be kept for building models 1 and 2)
#Length is the length of the genomes in base pairs. (Only genomes with length >= min_length will be considered during training and validation steps)
#PEG_Count is the number of protein encoding genes detected in the genome (Only genomes with PEG_Count >= min_peg_count will be considered during training and validation steps)
#Original_File is the name of the original fasta file in which the genomic sequence can be retrieved
#Recruitment_Criteria_Updated is the method used to identify the genome as a true porkatyote virus. Unless the genome has value of Cluster_Representative_RefSeq == TRUE all genomes for which the value of 

#Read in in_file to data0, Rename the first column to Host, set the rownames value of the data0$Sequence column and print a summary of the first seven colunmns
#in_file<-"/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/MMSeqs_Clusters/Filtered_Min_OG_1_Filtered_Scaffold_to_Domain_Score_Min_Score_50-Max_evalue_1e-05.tsv+Host_Taxonomy_6_Genus,Cluster_Representative_RefSeq,Length,PEG_Count,Original_File,Recruitment_Criteria_Updated.tsv"
in_file<-"/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/MMSeqs_Clusters/Filtered_0.9_OGs/Filtered_Correl_0.9_Min_OG_1_Filtered_Scaffold_to_Domain_Score_Min_Score_50-Max_evalue_1e-05.tsv+Host_Taxonomy_6_Genus,Cluster_Representative_RefSeq,Length,PEG_Count,Original_File,Recruitment_Criteria_Updated.tsv"

data0<-read.table(file=in_file,sep="\t",header=TRUE,quote="",comment="")

colnames(data0)[1]<-"Host"
rownames(data0)<-data0$Sequence
summary(data0[,1:7])

#################################################Benchmarking stage####################################################################
#####################Model 1: train on 80% of refseq and validate on the remaining 20%
#Filter to keep only viral RefSeq representative sequences with known host and faling within the esblished cutoffs for length and PEG count. Pass the filtered table to data1
data1<-data0[which(data0$Cluster_Representative_RefSeq == TRUE & data0$Host != "Unknown" & data0$Length >= min_length & data0$PEG_Count >= min_peg_count),]
summary(data1[,1:7])

#remove from data1 metadata columns that should not be used for training
data1<-subset(data1, select=-c(Cluster_Representative_RefSeq,Length,PEG_Count,Sequence,Original_File,Recruitment_Criteria_Updated))

#Assign to ColsToKeepM1 the vector of predictors in at least the number of min_dom_prevalence genomes
ColsToKeepM1<-as.vector(names(which(colSums(data1[,-1] != 0) >= min_dom_prevalence)))
#Keep in data1 only the predictors in ColsToKeepM1
data1<-subset(data1, select=c("Host",ColsToKeepM1))


for (k in 1:kfolds) {
	#Assign to train randomly selected row indexes. The rows present in train will be the ones in the training dataset and the ones absent will be the ones in the validation dataset
	train <- sample(nrow(data1), train_frac*nrow(data1), replace = FALSE)
	#Split the data into training and validation sets based on the values of train
	TrainSet <- data1[train,]
	ValidSet <- data1[-train,]

	#Convert the Host column of TrainSet and ValidSet into factors because this seems to be the only way the random forest will run
	TrainSet$Host<-factor(as.character(TrainSet$Host))
	ValidSet$Host<-factor(as.character(ValidSet$Host))

	#Set seed to always have the same output with the same input file
	#Build the model using number of trees specified by tree_num. Variables to try from var_tries. Set write.forest to TRUE so that the model can later be used for predictions. 
	start_time <- Sys.time()
	set.seed(100)
	#ranger_model_1<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE,importance="permutation",probability=TRUE,num.threads=threads)
	ranger_model_1<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE,importance="none",min.node.size=node_size,probability=TRUE,num.threads=threads)
	end_time <- Sys.time()
	end_time - start_time

	#Save image
	#save.image(file="/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/MMSeqs_Clusters/MMSeqs_Clusters_Ranger_Model_1.RData")

	#Checking classification accuracy for the training set
	predTrain <- predict(ranger_model_1, TrainSet[, colnames(TrainSet) != "Host"], type = "response", num.threads=threads)
	#Since probabilistic forest was built, predTrain$predictions hold the probabbility of each genome to be assigned to each class
	predVotes <- predTrain$predictions
	ColsMax<-max.col(predVotes)

	#Identify the winner prediction, i.e. the class with highest probabbility
	WinnerPreds<-c()
	for (i in ColsMax) {	
		WinnerPreds<-c(WinnerPreds,colnames(predVotes)[i])
	}

	#Calculate the accuracy by comparing the winner predictions to the value of TrainSet$Host
	accuracy_train<-mean(as.vector(WinnerPreds) == TrainSet$Host)

	#Build a matrix with the probabilities, the winner prediction and the true host
	pred_results_complement<-cbind(rownames(TrainSet),as.vector(WinnerPreds),as.vector(TrainSet$Host))
	colnames(pred_results_complement)<-c("Sequence","Predicted_Host","True_Host")
	pred_results_complete<-as.data.frame(cbind(predVotes,pred_results_complement))

	pred_results_complete$Veredict<-FALSE
	pred_results_complete$Winner_Score<-0

	for (i in 1:nrow(pred_results_complete)) {
		if (as.vector(pred_results_complete$Predicted_Host[i]) == as.vector(pred_results_complete$True_Host[i])) { pred_results_complete$Veredict[i]<-TRUE }
		pred_results_complete$Winner_Score[i]<-max(predVotes[i,])
	}

	accuracy_train_filtered<-mean(pred_results_complete$Veredict[which(pred_results_complete$Winner_Score >= valid_filter_votes)])  

	recall_train_filtered<-length(which(pred_results_complete$Winner_Score >= valid_filter_votes))/nrow(pred_results_complete)

	file_name<-paste(prefix,k,"Fold","Training_Ranger_Votes_Matrix_MMSeqs_Clusters_Model_1_Filtered_RefSeq_0.8.tsv",sep="_")
	write.table(pred_results_complete,file=file_name,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


	#Print table of confusion matrix
	file_name<-paste(prefix,k,"Fold","Training_Ranger_Model_1_Train_MMSeqs_Clusters_Filtered_RefSeq_0.8_Confusion_Matrix.tsv",sep="_")
	write.table(ranger_model_1$confusion.matrix,file=file_name,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

	#Print table of variable importance
	#write.table(ranger_model_1$variable.importance,file="Importance_Matrix_Training_Ranger_Model_1_Train_MMSeqs_Clusters_Filtered_RefSeq_0.8.tsv",sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

	#Generate ROC Curve
	file_name<-paste(prefix,k,"Fold","Training_Rranger_ROC_Curve_MMSeqs_Clusters_Model_1_Train_MMSeqs_Clusters_Filtered_RefSeq_0.8.pdf",sep="_")
	pdf(file_name,height=7,width=7,pointsize=8)
	roc<-roc.curve(pred_results_complete$Winner_Score,weights.class0=pred_results_complete$Veredict,curve=T)       
	plot(roc)
	dev.off()


	####################Estimations for the validation dataset
	predValid <- predict(ranger_model_1, ValidSet[, colnames(TrainSet) != "Host"], type = "response", num.threads=threads)
	predVotes <- predValid$predictions
	ColsMax<-max.col(predVotes)
	WinnerPreds<-c()
	for (i in ColsMax) {	
		WinnerPreds<-c(WinnerPreds,colnames(predVotes)[i])
	}
	accuracy_valid<-mean(as.vector(WinnerPreds) == ValidSet$Host)

		pred_results_complement<-cbind(rownames(ValidSet),as.vector(WinnerPreds),as.vector(ValidSet$Host))
		colnames(pred_results_complement)<-c("Sequence","Predicted_Host","True_Host")
		pred_results_complete<-as.data.frame(cbind(predVotes,pred_results_complement))

		pred_results_complete$Veredict<-FALSE
		pred_results_complete$Winner_Score<-0

		for (i in 1:nrow(pred_results_complete)) {
			if (as.vector(pred_results_complete$Predicted_Host[i]) == as.vector(pred_results_complete$True_Host[i])) { pred_results_complete$Veredict[i]<-TRUE }
			pred_results_complete$Winner_Score[i]<-max(predVotes[i,])
		}

		accuracy_valid_filtered<-mean(pred_results_complete$Veredict[which(pred_results_complete$Winner_Score >= valid_filter_votes)])  

		recall_valid_filtered<-length(which(pred_results_complete$Winner_Score >= valid_filter_votes))/nrow(pred_results_complete)
		
		file_name<-paste(prefix,k,"Fold","Validation_Ranger_Votes_Matrix_MMSeqs_Clusters_Model_1_Filtered_RefSeq_0.8.tsv",sep="_")
		write.table(pred_results_complete,file=file_name,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


	################Print summary of run
	params_summary<-rbind(c(in_file,k,node_size,min_dom_prevalence,train_frac,ranger_model_1$num.trees,ranger_model_1$mtry,accuracy_train,accuracy_valid,valid_filter_votes,accuracy_valid_filtered,recall_valid_filtered))
	colnames(params_summary)<-c("Input_File","K-fold","Node_Size","Min_Domain_Prevalence","Training_Fraction","Number_of_Trees","Variable_Tries","Accuracy_Training_Set","Accuracy_Validation_Set_Full","Votes_Winner_Filtering_Cutoff","Accuracy_Validation_Set_Filtered","Validation_Set_Filtered_Recall")

	file_name<-paste(prefix,"HP_Ranger_Parameter_Summary.tsv",sep="_")
	write.table(params_summary,file=file_name,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
}

################Model 2: train on 100% of refseq and validate on GLUVAB
train_frac<-1
#Filter to keep only viral RefSeq sequences with known host and faaling within the esblished cutoffs for length and PEG count
data1<-data0[which(data0$Cluster_Representative_RefSeq == TRUE & data0$Host != "Unknown" & data0$Length >= min_length & data0$PEG_Count >= min_peg_count),]

#remove columns that should not be used for training
data1<-subset(data1, select=-c(Cluster_Representative_RefSeq,Length,PEG_Count,Sequence,Original_File,Recruitment_Criteria_Updated))
ColsToKeepM2<-as.vector(names(which(colSums(data1[,-1] != 0) >= min_dom_prevalence)))
data1<-subset(data1, select=c("Host",ColsToKeepM2))
#Remove domains with present in less than 10 genomes
TrainSet <- data1

data2<-data0[which(data0$Original_File != "NA" & data0$Host != "Unknown" & data0$Length >= min_length & data0$PEG_Count >= min_peg_count),]
data2<-subset(data2, select=c("Host",ColsToKeepM2))

ValidSet <- data2

TrainSet$Host<-factor(as.character(TrainSet$Host))
ValidSet$Host<-factor(as.character(ValidSet$Host))


start_time <- Sys.time()
set.seed(100)
ranger_model_2<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE, importance="none",probability=TRUE,min.node.size=node_size,num.threads=threads)
end_time <- Sys.time()

end_time - start_time

# Checking classification accuracy for the training set
predTrain <- predict(ranger_model_2, TrainSet[, colnames(TrainSet) != "Host"], type = "response", num.threads=threads)
predVotes <- predTrain$predictions
ColsMax<-max.col(predVotes)

WinnerPreds<-c()
for (i in ColsMax) {	
	WinnerPreds<-c(WinnerPreds,colnames(predVotes)[i])
}

#WinnerPreds <- levels(TrainSet$Host)[max.col(predVotes)]
accuracy_train<-mean(as.vector(WinnerPreds) == TrainSet$Host)

pred_results_complement<-cbind(rownames(TrainSet),as.vector(WinnerPreds),as.vector(TrainSet$Host))
colnames(pred_results_complement)<-c("Sequence","Predicted_Host","True_Host")
pred_results_complete<-as.data.frame(cbind(predVotes,pred_results_complement))

pred_results_complete$Veredict<-FALSE
pred_results_complete$Winner_Score<-0

for (i in 1:nrow(pred_results_complete)) {
	if (as.vector(pred_results_complete$Predicted_Host[i]) == as.vector(pred_results_complete$True_Host[i])) { pred_results_complete$Veredict[i]<-TRUE }
	pred_results_complete$Winner_Score[i]<-max(predVotes[i,])
}

accuracy_train_filtered<-mean(pred_results_complete$Veredict[which(pred_results_complete$Winner_Score >= valid_filter_votes)])  

recall_train_filtered<-length(which(pred_results_complete$Winner_Score >= valid_filter_votes))/nrow(pred_results_complete)

file_name<-paste(prefix,"Training_Ranger_Votes_Matrix_MMSeqs_Clusters_Model_2_Filtered_RefSeq_1.tsv",sep="_")
write.table(pred_results_complete,file=file_name,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)


#Print table of confusion matrix
file_name<-paste(prefix,"Training_Ranger_model_2_Train_MMSeqs_Clusters_Filtered_RefSeq_1_Confusion_Matrix.tsv",sep="_")
write.table(ranger_model_2$confusion.matrix,file=file_name,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

#Print table of variable importance
#write.table(ranger_model_2$variable.importance,file="Importance_Matrix_Training_Ranger_model_2_Train_MMSeqs_Clusters_Filtered_RefSeq_1.tsv",sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

#Generate ROC Curve
file_name<-paste(prefix,"Training_Rranger_ROC_Curve_MMSeqs_Clusters_Model_2_Train_MMSeqs_Clusters_Filtered_RefSeq_1.pdf",sep="_")
pdf(file_name,height=7,width=7,pointsize=8)
roc<-roc.curve(pred_results_complete$Winner_Score,weights.class0=pred_results_complete$Veredict,curve=T)       
plot(roc)
dev.off()


####################Estimations for the validation dataset
predValid <- predict(ranger_model_2, ValidSet[, colnames(TrainSet) != "Host"], type = "response", num.threads=threads)
predVotes <- predValid$predictions
ColsMax<-max.col(predVotes)
WinnerPreds<-c()
for (i in ColsMax) {	
	WinnerPreds<-c(WinnerPreds,colnames(predVotes)[i])
}
accuracy_valid<-mean(as.vector(WinnerPreds) == ValidSet$Host)

	pred_results_complement<-cbind(rownames(ValidSet),as.vector(WinnerPreds),as.vector(ValidSet$Host))
	colnames(pred_results_complement)<-c("Sequence","Predicted_Host","True_Host")
	pred_results_complete<-as.data.frame(cbind(predVotes,pred_results_complement))

	pred_results_complete$Veredict<-FALSE
	pred_results_complete$Winner_Score<-0

	for (i in 1:nrow(pred_results_complete)) {
		if (as.vector(pred_results_complete$Predicted_Host[i]) == as.vector(pred_results_complete$True_Host[i])) { pred_results_complete$Veredict[i]<-TRUE }
		pred_results_complete$Winner_Score[i]<-max(predVotes[i,])
	}

	accuracy_valid_filtered<-mean(pred_results_complete$Veredict[which(pred_results_complete$Winner_Score >= valid_filter_votes)])  

	recall_valid_filtered<-length(which(pred_results_complete$Winner_Score >= valid_filter_votes))/nrow(pred_results_complete)

	file_name<-paste(prefix,"Validation_Ranger_Votes_Matrix_MMSeqs_Clusters_Model_2_Filtered_RefSeq_1.tsv",sep="_")
	write.table(pred_results_complete,file=file_name,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
	
	
################Print summary of run
params_summary<-rbind(c(in_file,NA,node_size,min_dom_prevalence,train_frac,ranger_model_2$num.trees,ranger_model_2$mtry,accuracy_train,accuracy_valid,valid_filter_votes,accuracy_valid_filtered,recall_valid_filtered))
colnames(params_summary)<-c("Input_File","K-fold","Node_Size","Min_Domain_Prevalence","Training_Fraction","Number_of_Trees","Variable_Tries","Accuracy_Training_Set","Accuracy_Validation_Set_Full","Votes_Winner_Filtering_Cutoff","Accuracy_Validation_Set_Filtered","Validation_Set_Filtered_Recall")

file_name<-paste(prefix,"HP_Ranger_Parameter_Summary.tsv",sep="_")
write.table(params_summary,file=file_name,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)


#######################################Build definitive model (Full RefSeq + PhaPhyCS)###############################################
data3<-data0[which(data0$Host != "Unknown" & data0$Length >= min_length & data0$PEG_Count >= 1),]
summary(data3[,1:7])
data3<-subset(data3, select=-c(Cluster_Representative_RefSeq,Length,PEG_Count,Sequence,Original_File,Recruitment_Criteria_Updated))

ColsToKeepM3<-as.vector(names(which(colSums(data3[,-1] != 0) >= min_dom_prevalence)))
write.table(ColsToKeepM3,file="HP_Ranger_Model_3_Valid_Cols.txt",sep="\t",append=FALSE,row.names=FALSE,col.names=FALSE,quote=FALSE)
data3<-subset(data3, select=c("Host",ColsToKeepM3))
TrainSet <- data3
TrainSet$Host<-factor(as.character(TrainSet$Host))


start_time <- Sys.time()
set.seed(100)
#ranger_model_3<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=TRUE, importance="permutation",num.threads=threads)
#ranger_model_3<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=TRUE,importance="permutation",local.importance=TRUE,min.node.size=node_size,num.threads=threads)
ranger_model_3<-ranger(y=TrainSet$Host, x=TrainSet[, colnames(TrainSet) != "Host"], num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=TRUE,importance="impurity",min.node.size=node_size,num.threads=threads)
end_time <- Sys.time()

end_time - start_time

save.image(file="/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/MMSeqs_Clusters/MMSeqs_Clusters_Ranger_Model_1+2+3.RData")

#Print table of confusion matrix
file_name<-paste(prefix,"Final_Ranger_Model_3_Train_MMSeqs_Clusters_Filtered_RefSeq_1+GLUVAB_Confusion_Matrix.tsv",sep="_")
write.table(ranger_model_3$confusion.matrix,file=file_name,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

#Print Global importance matrix
file_name<-paste(prefix,"Importance_Matrix_Final_Ranger_Model_3_Train_MMSeqs_Clusters_Filtered_RefSeq_1+GLUVAB.tsv",sep="_")
write.table(ranger_model_3$variable.importance,file=file_name,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

#Print Local importance matrix
#file_name<-paste(prefix,"Importance_Local_Matrix_Final_Ranger_Model_3_Train_MMSeqs_Clusters_Filtered_RefSeq_1+GLUVAB.tsv",sep="_")
#write.table(ranger_model_3$variable.importance.local,file=file_name,sep="\t",append=FALSE,row.names=TRUE,col.names=NA,quote=FALSE)

# Checking classification accuracy for the training set
predTrain <- predict(ranger_model_3, TrainSet[, colnames(TrainSet) != "Host"], type = "response", num.threads=threads)
predVotes <- predTrain$predictions
ColsMax<-max.col(predVotes)

WinnerPreds<-c()
for (i in ColsMax) {	
	WinnerPreds<-c(WinnerPreds,colnames(predVotes)[i])
}

accuracy_train<-mean(as.vector(WinnerPreds) == TrainSet$Host)

pred_results_complement<-cbind(rownames(TrainSet),as.vector(WinnerPreds),as.vector(TrainSet$Host))
colnames(pred_results_complement)<-c("Sequence","Predicted_Host","True_Host")
pred_results_complete<-as.data.frame(cbind(predVotes,pred_results_complement))

pred_results_complete$Veredict<-FALSE
pred_results_complete$Winner_Score<-0

for (i in 1:nrow(pred_results_complete)) {
	if (as.vector(pred_results_complete$Predicted_Host[i]) == as.vector(pred_results_complete$True_Host[i])) { pred_results_complete$Veredict[i]<-TRUE }
	pred_results_complete$Winner_Score[i]<-max(predVotes[i,])
}

accuracy_train_filtered<-mean(pred_results_complete$Veredict[which(pred_results_complete$Winner_Score >= valid_filter_votes)])  

recall_train_filtered<-length(which(pred_results_complete$Winner_Score >= valid_filter_votes))/nrow(pred_results_complete)

file_name<-paste(prefix,"Training_Ranger_Votes_Matrix_MMSeqs_Clusters_Model_3_Filtered_RefSeq_1+GLUVAB.tsv",sep="_")
write.table(pred_results_complete,file=file_name,sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

params_summary<-rbind(c(in_file,NA,node_size,min_dom_prevalence,train_frac,ranger_model_3$num.trees,ranger_model_3$mtry,accuracy_train,NA,valid_filter_votes,NA,NA))
colnames(params_summary)<-c("Input_File","K-fold","Node_Size","Min_Domain_Prevalence","Training_Fraction","Number_of_Trees","Variable_Tries","Accuracy_Training_Set","Accuracy_Validation_Set_Full","Votes_Winner_Filtering_Cutoff","Accuracy_Validation_Set_Filtered","Validation_Set_Filtered_Recall")

file_name<-paste(prefix,"HP_Ranger_Parameter_Summary.tsv",sep="_")
write.table(params_summary,file=file_name,sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)

#############Remove unnecessary variables and save the three models
rm(accuracy_train,accuracy_train_filtered,accuracy_valid,accuracy_valid_filtered,ColsMax,data0,data1,data2,data3,data4,end_time, filter_correl_cutoff,i,in_file,min_dom_prevalence,min_length,min_peg_count,params_summary,pred_results_complement,pred_results_complete,predTrain,predValid,predValidResults,predVotes,recall_train_filtered,recall_valid_filtered,roc,start_time,threads,train,train_frac,TrainSet,tree_num,valid_file1,valid_filter_votes,ValidSet,var_tries,WinnerPreds)

save.image(file="MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData")

#Read in the first validation file
#valid_file1<-"/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/Dzunkova_2019_Validation/Dzunkova_MMSeqs_M3_Scaffold_to_Domain_Score_Min_Score_50-Max_evalue_1e-05.tsv"
#valid_file1<-"/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/Arkhipova_2018_Validation/Arkhipova_MMSeqs_M3_Scaffold_to_Domain_Score_Min_Score_50-Max_evalue_1e-05.tsv"
valid_file1<-"/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/NCBI_NR_NonRefSeq_Genomes_Validation/NCBI_NR_Non_RefSeq_Genome_to_Domain_Score_Min_Score_50-Max_evalue_1e-05.tsv"


data4<-read.table(file=valid_file1,sep="\t",header=TRUE,quote="",comment="")
colnames(data4)[1]<-"Host"
rownames(data4)<-data4$Sequence
summary(data4[,1:7])
ValidSet<-subset(data4, select=c(ColsToKeepM3))

predValid <- predict(ranger_model_3, ValidSet, type = "response", num.threads=threads)

ColsMax<-max.col(predValid$predictions)
WinnerPreds<-c()
for (i in ColsMax) {	
	WinnerPreds<-c(WinnerPreds,colnames(predVotes)[i])
}

predValidResults<-cbind(data4[,1:7],WinnerPreds)
predValidResults$Veredict<-FALSE
predValidResults$Winner_Score<-0
predValidResults<-cbind(predValidResults,predValid$predictions)

for (i in 1:nrow(predValidResults)) {
	if (as.vector(predValidResults$WinnerPreds[i]) == as.vector(predValidResults$Host_6_Genus[i])) { predValidResults$Veredict[i]<-TRUE }
	predValidResults$Winner_Score[i]<-max(predValid$predictions[i,])
}


write.table(predValidResults,file="Final_Ranger_Model_3_Validation_MMSeqs_Clusters_Filtered_RefSeq_1+GLUVAB_Predictions_NCBI_NR_Non_RefSeq.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
#write.table(predValidResults,file="Final_Ranger_Model_3_Validation_MMSeqs_Clusters_Filtered_RefSeq_1+GLUVAB_Predictions_Arkhipova.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

################Print summary of run (needs editing)
params_summary<-rbind(c(in_file,NA,NA,min_dom_prevalence,train_frac,ranger_model_3$num.trees,ranger_model_3$mtry,accuracy_train,NA,valid_filter_votes,NA,NA))
colnames(params_summary)<-c("Input_File","Min_Col_Sum","Min_Row_Sum","Min_Domain_Prevalence","Training_Fraction","Number_of_Trees","Variable_Tries","Accuracy_Training_Set","Accuracy_Validation_Set_Full","Votes_Winner_Filtering_Cutoff","Accuracy_Validation_Set_Filtered","Validation_Set_Filtered_Recall")


write.table(params_summary,file="/home/rohit/felipe/Databases/RefSeqVir_Oct_19/RF_Host_Pred/External_Validation_Set/MMSeqs_Clusters/HP_Ranger_Parameter_Summary.tsv",sep="\t",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)
