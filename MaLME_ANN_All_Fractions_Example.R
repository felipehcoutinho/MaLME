library(scales)
library(foreach)
library(doParallel)
library(NeuralNetTools)
library(nnet)
library(Metrics)
library(dplyr)
	
max_opt_tries<-10
max_rmse<-0.2
n_neurons<-3

response_type<-"Prokaryotes"

if (response_type == "Eukaryotes") {
	#Eukaryotes
	predictor_variables<-c("Temperature","Salinity","Si","NO3","Phos","Fe")
	#Eukaryotes (939 x 14)
	predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Eukaryote_Communities_Sample_Metadata.tsv"
	#Eukaryotes (676 x194)
	response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Eukaryote_Communities_Fraction_20-180µm_Percentage_Abundance_Data.tsv"
} else if (response_type == "Prokaryotes") {
	#Prokaryotes
	predictor_variables<-c("Temperature","Salinity","PO4","Ammonium.5m","NO2NO3","Iron.5m")
	#Prokaryotes (180 x 38)
	predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Prokaryote_Communities_Sample_Metadata.tsv"
	#Prokaryotes (18703 x 187)
	response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Prokaryote_Communities_Percentage_Abundance_Data.tsv"
} else if (response_type == "Viruses") {
	#Viruses
	predictor_variables<-c("Temperature","Oxygen","Chlorophyll_A","Salinity","Iron_5m","Ammonium_5m")
	#Viruses (76753 x 139)
	predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Viral_Communities_Sample_Metadata.tsv"
	#Viruses (131 x 45)
	response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Virus_Communities_RPKM_Abundance_Data.tsv"
}


raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

raw_predictor_df<-subset(raw_predictor_df,select=c("MG_UID",predictor_variables))

raw_predictor_df<-raw_predictor_df[complete.cases(raw_predictor_df),]

#Build the scenarios of climate change for the raw predictors
num_cols<-nums<-unlist(lapply(raw_predictor_df, is.numeric), use.names = FALSE)

z_predictor_df<-raw_predictor_df

#summary(z_predictor_df)

#Load and filter response data
raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

num_col_idx<-unlist(lapply(raw_response_df, is.numeric), use.names = FALSE)

MG_UIDs<-colnames(raw_response_df)[colnames(raw_response_df) %in% z_predictor_df$MG_UID]

response_variables<-rownames(raw_response_df)[which(rowSums(raw_response_df[,num_col_idx]) > 0)]

#response_variables<-c("AACY020065672.95.1597")
response_variables<-response_variables[5:15]

z_response_df<-raw_response_df

t_z_response_df<-as.data.frame(t(subset(z_response_df,select=MG_UIDs)))
t_z_response_df$MG_UID<-as.factor(MG_UIDs)

dim(t_z_response_df)
summary(t_z_response_df[,1:11])

#my.cluster <- parallel::makeCluster(47, type = "PSOCK")
#doParallel::registerDoParallel(cl = my.cluster)

set.seed(666)
#all_results_list<-foreach(rvar = response_variables) %dopar% {
all_results_list<-foreach(rvar = response_variables) %do% {
	print(paste("Building ANNs for: ",rvar),sep="")
	z_train_df<-merge(subset(t_z_response_df,select=c(rvar,"MG_UID")),subset(z_predictor_df,select=c("MG_UID",predictor_variables)), by="MG_UID",all.x=TRUE)
	#dim(z_train_df)
	#summary(z_train_df)
	rownames(z_train_df)<-z_train_df$MG_UID
	z_train_df$MG_UID<-NULL
	z_train_df$Scenario<-NULL
	colnames(z_train_df)[which(colnames(z_train_df) == rvar)]<-"Z_Response_Variable"
	#dim(z_train_df)
	summary(z_train_df)
	
	passed<-FALSE
	best_ANN<-NA
	optimization_num<-0
	best_rmse<-+Inf
	best_pearson_cor_t<-NA
	while ((best_rmse > max_rmse) & (optimization_num < max_opt_tries)) {
		optimization_num<-optimization_num+1
		shouter<-paste("Response",rvar,"Iteration",optimization_num,sep=" ")
		#print(shouter)
		trained_net<-nnet(Z_Response_Variable ~ .,z_train_df,size=n_neurons,linout=TRUE,maxit=1000,decay=0.001,reltol=0.001)
		
		if ((sd(trained_net$fitted.values) == 0) | (sd(z_train_df$Z_Response_Variable) == 0)) {next}
		
		pearson_cor_t<-cor.test(z_train_df$Z_Response_Variable,trained_net$fitted.values,method="pearson")
		
		rmse_t<-rmse(z_train_df$Z_Response_Variable,trained_net$fitted.values)
		
		if (rmse_t < best_rmse) {
			best_ANN<-trained_net
			best_rmse<-rmse_t
			best_pearson_cor_t<-pearson_cor_t$estimate
			passed<-TRUE
		}
	}
	
	if (passed == TRUE) {
		importance_olden<-olden(best_ANN,bar_plot=FALSE)
		importance_olden$Normalized_Importance<-importance_olden$importance/max((abs(importance_olden$importance)))
		importance_olden$Predictor_Variable<-rownames(importance_olden)

		scenario_predictions<-as.data.frame(predict(best_ANN,z_predictor_df,type=c("raw")))
		colnames(scenario_predictions)[1]<-"Predicted_Z_Response_Variable"
		scenario_predictions$Scenario<-z_predictor_df$Scenario#scenario_names

		scenario_predictions$Response_Variable_Name<-rvar
		scenario_predictions$MG_UID<-z_predictor_df$MG_UID
	} else {
		importance_olden<-NA
		scenario_predictions<-NA
	}
	
	print(paste("Successfully built ANNs for: ",rvar),sep="")
	
	results_list<-list("Respose_Variable_Name" = rvar, "ANN" = best_ANN, "RMSE" = best_rmse, "PCC" = best_pearson_cor_t, "Olden_Importance" =  importance_olden, "Scenario_Predictions" = scenario_predictions)
}

#parallel::stopCluster(cl = my.cluster)

print("Finished building Networks")

quit(status=1)