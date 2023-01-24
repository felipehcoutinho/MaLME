library(ranger)
library(scales)
library(Metrics)

options(expressions = 5e5)

###Step 6 train networks with full data:
threads<-48
tree_num<-10000 #Number of trees to grow
node_size<-1 #Node size for the probabilistic forest
#Cutoff of correlation for considering an ANN valid
cor_cutoff<-0.7
#Max RMSE cutoff for considering an ANN valid
max_rmse<-+Inf
#Variables used for prediction
response_variables<-c("Temperature","Salinity","Si","NO3","Phos","Fe")
#passed_importance_vals_garson<-c()
passed_importance_vals<-c()
#Fraction of data to be used for training
#tfrac<-0.5
#type if data trasnformation to perfom
transformation<-"none"#"none"#"zscore"#
best_rfs<-list()


#Load, filter, and scale response data
predictor_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Eukaryote_Communities_Fraction_20-180µm_Percentage_Abundance_Data.tsv"

raw_predictor_df<-read.table(file=predictor_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,row.names=1,check.names=FALSE)

num_col_idx<-unlist(lapply(raw_predictor_df, is.numeric), use.names = FALSE)

MG_UIDs<-colnames(raw_predictor_df)[num_col_idx]

predictor_variables<-rownames(raw_predictor_df)

t_raw_predictor_df<-as.data.frame(t(raw_predictor_df))
t_raw_predictor_df<-as.data.frame(apply(t_raw_predictor_df[num_col_idx,],2,as.numeric))
t_raw_predictor_df$MG_UID<-as.factor(MG_UIDs)

summary(t_raw_predictor_df)#[,1:10])

predictor_variables<-colnames(t_raw_predictor_df)[num_col_idx]

#Load and prepare response data
response_file<-"/mnt/lustre/scratch/fcoutinho/MaLME/Marine_Eukaryote_Communities_Sample_Metadata.tsv"

raw_response_df<-read.table(file=response_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE,check.names=FALSE)

raw_response_df<-raw_response_df[which(raw_response_df$MG_UID %in% MG_UIDs),]

dim(raw_response_df)
summary(raw_response_df)

#Merge response and predictor tables

full_lsubtdata<-merge(raw_response_df,subset(t_raw_predictor_df,select=c("MG_UID",predictor_variables)), by="MG_UID",all.x=TRUE)

var_tries<-round(length(predictor_variables)/3 )#Variables to try at each split when growing trees

all_best_stats<-c()
all_rf_stats<-c()
best_rf_stats<-c()
all_preds<-rbind()

set.seed(3435)
for (resp_var in response_variables) {
	#best_prod<-0
	best_rmse<-+Inf
	#best_pred<-NA
	best_RF<-NA
	#best_tmatrix<-NA
	#best_vmatrix<-NA
	all_cort<-c()
	#all_corv<-c()
	all_rmset<-c()
	#all_rmsev<-c()
	pass_cort<-c()
	#pass_corv<-c()
	pass_rmset<-c()
	#pass_rmsev<-c()
	best_stats<-c()
	
	lsubtdata<-full_lsubtdata
	lsubtdata<-lsubtdata[complete.cases(lsubtdata),]
	Response<-lsubtdata[[resp_var]]
	MG_UIDs<-as.vector(lsubtdata$MG_UID)
	lsubtdata<-subset(lsubtdata,select=predictor_variables)
	#summary(lsubtdata)
	#colnames(lsubtdata)[1:25]
	#summary(Response)
	
		print(paste("Building RF for response variable:",resp_var,sep=" "))
		#Calculate RF with random starting weights
		rf_model<-ranger(y=Response, x=lsubtdata, num.trees=tree_num, mtry=var_tries, write.forest=TRUE,probability=FALSE,importance="impurity",min.node.size=node_size,num.threads=threads)
		
		#Calculate pearson correlation for the measured and predicted values
		cor_t<-cor.test(rf_model$predictions,Response,method="pearson")
		#cor_v<-cor_t
		#Calculate Root Mean Square Error
		rmse_t<-rmse(rf_model$predictions,Response)
		#rmse_v<-rmse_t
		#Keep track of cort_t cor_v and prod values
		all_cort<-c(all_cort,cor_t$estimate)
		#all_corv<-c(all_corv,cor_v$estimate)

		all_rmset<-c(all_rmset,rmse_t)
		#all_rmsev<-c(all_rmsev,rmse_v)
		
		if ((cor_t$estimate < cor_cutoff) || (rmse_t > max_rmse)) {passed<-FALSE } else {passed<-TRUE}
		
		all_rf_stats<-rbind(all_rf_stats,c(resp_var,1,cor_t$estimate,rmse_t,passed))
		
		#Check if correlation and RMSE values are within allowed cutoff
		if (passed == FALSE) {
		pass_cort<-c(pass_cort,NA)
		#pass_corv<-c(pass_corv,NA)
		pass_rmset<-c(pass_rmset,NA)
		#pass_rmsev<-c(pass_rmsev,NA)
		} else {
		pass_cort<-c(pass_cort,cor_t$estimate)
		#pass_corv<-c(pass_corv,cor_v$estimate)
		pass_rmset<-c(pass_rmset,rmse_t)
		#pass_rmsev<-c(pass_rmsev,rmse_v)
		}
		#Keep track of the impotance values assigned to predictors for all the passed networks
		passed_importance<-as.data.frame(cbind(names(rf_model$variable.importance),rf_model$variable.importance))
		colnames(passed_importance)<-c("Predictor","importance")
		passed_importance$importance<-as.numeric(passed_importance$importance)
		passed_importance$Response<-resp_var
		passed_importance$RMSEV<-rmse_t
		passed_importance$CorV<-cor_t$estimate
		passed_importance$Method<-"Impurity"
		passed_importance$Normalized_Importance<-passed_importance$importance/max((abs(passed_importance$importance)))
		passed_importance_vals<-rbind(passed_importance_vals,passed_importance)


		#Store best matrix and best ANN based on RMSE
		if ((rmse_t <= best_rmse) & (cor_t$estimate >= cor_cutoff)) {
			best_rmse<-rmse_t
			best_RF<-rf_model
			best_stats<-passed_importance
		}
		
	all_best_stats<-rbind(all_best_stats,best_stats)
	if (is.na(best_RF)[1] == FALSE) {	
		best_rfs[[resp_var]]<-best_RF
		pred_data<-as.data.frame(cbind(rf_model$predictions,Response))
		colnames(pred_data)<-c("Predicted","Measured")
		pred_data$Response_Var<-resp_var
		pred_data$MG_UID<-as.factor(MG_UIDs)
		all_preds<-rbind(all_preds,subset(pred_data,select=c("Measured","Predicted","Response_Var","MG_UID")))

	}
}

colnames(all_rf_stats)<-c("Response","Iteration","CorV","RMSEv","Passed")

head(all_rf_stats)

write.table(all_rf_stats,file="Reverse_Ecology_All_RF_Stats.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)

#load("Reverse_Ecology_RF.RData")

make_plots<-FALSE
if (make_plots == TRUE) {
	library(ggplot2)
	library(gplots)
	library(reshape2)
	library(RColorBrewer)
	#library("GGally")
	library(ggpubr, lib.loc="/mnt/lustre/bio/users/fcoutinho/Rlibs/")
	#library(vegan)


	#Plot measured x predicted values for each response variable in scatterplot
	fig_1A<-ggplot(all_preds,aes(y=Predicted,x=Measured))+geom_point(alpha=0.5)+geom_smooth(method = "lm",formula= y~x, se = TRUE)+facet_wrap(Response_Var ~ ., scales="free",nrow=1)+theme_bw()
		
	#Read in scaffold info
	response_info_file<-"/mnt/lustre/scratch/fcoutinho/StG/Min_Count_1k_Cluster_95ID_80Cov_GLUVAB3_VIBRANT_Genomes_Min_5kbp_Full_Seq_Info.tsv"
	scaffold_data<-read.table(file=response_info_file,sep="\t",header=TRUE,quote="",comment="",stringsAsFactors=TRUE)
	scaffold_data$Predictor<-scaffold_data$Sequence

	all_best_stats<-merge(all_best_stats,scaffold_data, by="Predictor",all.x=TRUE)
	
	f_best_stats<-all_best_stats[which(all_best_stats$Normalized_Importance > 0.01),]
	
	write.table(f_best_stats,file="Reverse_Ecology_Filtered_Predictior_Stats.tsv",sep="\t",append=FALSE,row.names=FALSE,col.names=TRUE,quote=FALSE)
	
	host_f_all_best_stats<-aggregate(all_best_stats$Normalized_Importance, by=list(Host_Phylum = all_best_stats$Predicted_Host_2_Phylum, Response=all_best_stats$Response), FUN=sum)

	colnames(host_f_all_best_stats)[3]<-"Sum_of_Normalized_Importances"
		
	host_f_all_best_stats<-host_f_all_best_stats[which(host_f_all_best_stats$Sum_of_Normalized_Importances > 0.1),]
	
	fig_1C<-ggplot(host_f_all_best_stats)+geom_bar(aes(x=Host_Phylum, y=Sum_of_Normalized_Importances, fill=Response),position="dodge",stat="identity",colour="black",alpha=0.9)+theme_bw()+theme(legend.position="top",axis.text.x = element_text(angle = 45,hjust = 1))
	#+facet_grid(. ~ Response, scales="free_x")#+scale_fill_manual(name="Host_Phylum",values=host_coloring)
	#ggplot(mdata, aes(x=Sample, y=Abundance, fill=Pathway))+geom_bar(position="stack",stat="identity",colour="black",alpha=0.9)+ylab("RPKM")+xlab("Station")+theme_bw()+theme(legend.position="right",text = element_text(size = 18),legend.key.size = unit(0.5, "cm"),axis.text.x = element_text(angle = 45,hjust = 1,size=9))
	
	#f_all_best_stats<-all_best_stats[which(all_best_stats$Normalized_Importance > 0.05),]

	#Plot importance data summarized by host phylum (Best ANNs)
	#fig_1C<-ggplot(f_all_best_stats, aes(x=Normalized_Importance, fill=Response))+geom_histogram(colour="black",bins=10)+ylab("Frequency")+theme_bw()+theme(axis.text.x = element_text(size =12),axis.text.y = element_text(size=14),axis.title=element_text(size=16),legend.text=element_text(size=16),legend.position="none")+facet_grid(Predicted_Host_2_Phylum ~ Response, scales="free_x")

	#host_f_all_best_stats<-f_all_best_stats[which(!is.na(f_all_best_stats$Predicted_Host_2_Phylum)),]

	#fig_1C<-ggplot(host_f_all_best_stats)+geom_boxplot(aes(fill=Response, y=Normalized_Importance,x=Predicted_Host_2_Phylum))+theme_bw()+theme(legend.position="top",axis.text.x = element_text(angle = 45,hjust = 1))+facet_grid(. ~ Response, scales="free_x")#+scale_fill_manual(name="Host_Phylum",values=host_coloring)

	#Plot importance data summarized by viral family(Best ANNs)

	#fig_1B<-ggplot(f_all_best_stats, aes(x=Normalized_Importance, fill=Response))+geom_histogram(colour="black",bins=10)+ylab("Frequency")+theme_bw()+theme(axis.text.x = element_text(size =12),axis.text.y = element_text(size=14),axis.title=element_text(size=16),legend.text=element_text(size=16),legend.position="none")+facet_grid(VPF_family ~ Response, scales="free_x")
	
	#fam_f_all_best_stats<-f_all_best_stats[which(!is.na(f_all_best_stats$VPF_family)),]

	#fig_1B<-ggplot(fam_f_all_best_stats)+geom_boxplot(aes(fill=Response, y=Normalized_Importance,x=VPF_family))+theme_bw()+theme(legend.position="top",axis.text.x = element_text(angle = 45,hjust = 1))+facet_grid(. ~ Response, scales="free_x")#+scale_fill_manual(name="Family",values=vir_fam_coloring)
	
	fam_f_all_best_stats<-aggregate(all_best_stats$Normalized_Importance, by=list(Family = all_best_stats$VPF_family, Response=all_best_stats$Response), FUN=sum)

	colnames(fam_f_all_best_stats)[3]<-"Sum_of_Normalized_Importances"
	
	fam_f_all_best_stats<-fam_f_all_best_stats[which(fam_f_all_best_stats$Sum_of_Normalized_Importances > 0.1),]
	
	fig_1B<-ggplot(fam_f_all_best_stats)+geom_bar(aes(x=Family, y=Sum_of_Normalized_Importances, fill=Response),position="dodge",stat="identity",colour="black",alpha=0.9)+theme_bw()+theme(legend.position="top",axis.text.x = element_text(angle = 45,hjust = 1))

	fig_1<-ggarrange(fig_1A, fig_1B, fig_1C, nrow = 3, labels = "AUTO", heights=c(1,2,2)) 
	ggsave("StG_RF_Reverse_Ecology_Figure_1.pdf",plot=fig_1,width=12,height=15,pointsize=8)
	
}

#save.image("Reverse_Ecology_RF.RData")

quit(status=1)