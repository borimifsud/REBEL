################################################################
################################################################
# XGBoost code example used in Villiers et al 2021
# The code below runs xgboost models to predict a fragments regulatory category 
# (GainedUp, GainedDown, GainedNoChange, LostUp, LostDown, LostNoChange)
# based on its motif composition.
# The code takes the ATAC-seq regions, performs a motif 
# search to create a 436 vector representing the motif compostion of the fragment
# then assigns each fragment to a regulatory outcome.
# For each category, a dataset is created combining all the fragments for that feature
# with an equal number of random other fragments from the other regulatory categories
# this is split into a train and test set (80:20)
# The second part of this code takes each category, splits the data into 5 equal groupings
# uses the model parameters from the first AUC score to predict each of the 5 groupings.
# Each equal grouping and their model is input into SHAP to generate SHAP scores.
# These SHAP scores are used for analysis in the paper section: 
# Machine learning algorithms identify distinct TF motif compositions that define responses to PML-RARA binding.
################################################################
################################################################

		# load libraries and processed data
		# load in all the pre-processed datasets in R object form 
		# from ./code/U937_PR9_Processed_Data_Descriptions.R 
		source("./code/U937_PR9_Processed_Data_Descriptions.R")
		graphs.dir <- "./graphs/"

		suppressPackageStartupMessages({
			library(stringr)
			library(GenomicRanges)
			library(limma)
			library(GenomicInteractions)
			library(foreach)
			library(doParallel)

			library("xgboost")  # the main algorithm
			library("archdata") # for the sample dataset
			library("caret")    # for the confusionmatrix() function (also needs e1071 package)
			library("dplyr")    # for some data preperation
			library("Ckmeans.1d.dp") # for xgb.ggplot.importance
			require(methods)

			library("SHAPforxgboost")
			library("ggplot2")
			library("data.table")
			library("here")
	
		})


################################################################
# Start the motif calling of ATAC-seq peaks and peak annotations
################################################################

		#setup parallel backend to use x processors
		cores=detectCores()
		cl <- makeCluster(cores[1]-1) #not to overload your computer
		registerDoParallel(cl)

		#function to get ATAC-seq peaks into format for HOMER motif calling
		GetHomerReady <- function(GR,name){
				GR <- as.data.frame(GR)[,c(1,2,3,5,4)]
				GR[,4] <- 1:length(GR[,4])
				write.table(GR,file=paste0(name,".txt"), sep = "\t", quote = F, row.names = F, col.names = F)
			}

		GetHomerReady(Final_Combined_ATAC,"Final_Combined_ATAC")

		#command line (only part) requires the homer suite to be downloaded: http://homer.ucsd.edu/homer/
		# All_HOMER_Motifs.txt (included) = a file with all motif PWM matrices in a format Homer requires (436 motifs)
		# -200,200 specifies we want to look only at 200bp right and left of the peak centre
		# PATH=$PATH:/Users/william/Tools/Homer/bin/
		#
		# annotatePeaks.pl Final_Combined_ATAC.txt hg19 -m ./code/All_HOMER_Motifs.txt > Final_Combined_ATAC_440_200_200.txt -size -200,200 -noann
		# The output from homer (Final_Combined_ATAC_440_200_200.txt) requires reformatting
		# so requires a bit of fiddling to get it into an interpretable format
		# see http://homer.ucsd.edu/homer/ngs/annotation.html 
		# for homer output and info
		# Final_Combined_ATAC_440_200_200.txt is supplied


		# next step is to create the matrix used for prediction using the motifs called on each ATAC-seq fragment - file is a subset of the first 8000 ATAC-seq fragments
		ATAC_occurence1 <- read.csv('Final_Combined_ATAC_440_200_200_subset1.txt.gz', sep = "\t", head = T)
        ATAC_occurence2 <- read.csv('Final_Combined_ATAC_440_200_200_subset2.txt.gz', sep = "\t", head = T)
        ATAC_occurence3 <- read.csv('Final_Combined_ATAC_440_200_200_subset3.txt.gz', sep = "\t", head = T)

            ATAC_occurence <- rbind(ATAC_occurence1, ATAC_occurence2)
            ATAC_occurence <- rbind(ATAC_occurence, ATAC_occurence3)
			ATAC_occurence <- ATAC_occurence[order(ATAC_occurence[,1]),]
			colnames(ATAC_occurence) <- toupper(colnames(ATAC_occurence))
			ATAC_occurence[is.na(ATAC_occurence)] <- 0
			ATAC_occurence[ATAC_occurence == ""] <- 0
			haha <- as.data.frame(ATAC_occurence[,22:457])
			colnames(haha) <- substr(colnames(haha),start=1,stop=5)
			NAMES <- colnames(haha)

		#extract from the homer table the number of times the motif comes up on each fragment
		DAT_where <- list()
		DAT_list <- list()
		
				for(i in 1:length(row.names(haha))){
	 
					 DD <- haha[i,]
					DAT_list[[i]] <- as.numeric(apply(DD,1,function(x){str_count(x,"[(]")}))
					DAT_where[[i]] <- as.numeric(apply(DD,1,function(x){str_extract(x, "[^(]+")}))

					cat("done",i,"of",length(row.names(haha)),"\n")
				}
		
		
		#add the numer of occurences to the Final_Combined_ATAC GR object
		Final_Combined_ATAC <- Final_Combined_ATAC[1:8000,]
        Final_Combined_ATAC$DAT_where <- DAT_where
		Final_Combined_ATAC$NUM <- DAT_list


		#Assign an ATAC-seq peak (with the motif info) to
		#problem= 1 HindIII fragment has more than 1 ATAC-seq peak

		Differential <- c(CONSENSUS_GI, CONSENSUSDEC_GI)

		#roughly assign any of the ATAC-seq peaks to a differential interaction end (so we know the ones to consider downstream)		
		Differential$ATAC1 <- 'NA'
		Differential$ATAC1[queryHits(findOverlaps(anchorOne(Differential),Final_Combined_ATAC))] <- subjectHits(findOverlaps(anchorOne(Differential),Final_Combined_ATAC))
		Differential$ATAC2 <- 'NA'
		Differential$ATAC2[queryHits(findOverlaps(anchorTwo(Differential),Final_Combined_ATAC))] <- subjectHits(findOverlaps(anchorTwo(Differential),Final_Combined_ATAC))
		Differential$PR1 <- 'NA'
		Differential$PR1[queryHits(findOverlaps(anchorOne(Differential),SEACR_Peaks))] <- subjectHits(findOverlaps(anchorOne(Differential),SEACR_Peaks))
		Differential$PR2 <- 'NA'
		Differential$PR2[queryHits(findOverlaps(anchorTwo(Differential),SEACR_Peaks))] <- subjectHits(findOverlaps(anchorTwo(Differential),SEACR_Peaks))

		#which Differential interactions have an ATAC-seq peak on BOTH ends (these are the ones were interested in)
		Differential <- Differential[which(Differential$ATAC1 != 'NA')]
		Differential <- Differential[which(Differential$ATAC2 != 'NA')]

		#subset Final_Combined_ATAC to those fragments just overlapping PML-RARA bound regions
		Differential <- Differential[c(which(Differential$PR1 != 'NA'),which(Differential$PR2 != 'NA'))]

		#Want to identify each ATAC-seq peak that overlaps with each interaction end
		#starting with the interaction end 1 (bait1)
	
			#bait1
			Doh1 <- Differential
			ALL_O1 <- findOverlaps(anchorOne(Doh1),Final_Combined_ATAC)
			ALL_Q1 <- queryHits(ALL_O1)
			ALL_S1 <- subjectHits(ALL_O1)
			Doh1 <- Doh1[ALL_Q1]
			Doh1$ATAC <- ALL_S1
			Doh1$Differential <- ALL_Q1
			Doh1$End <- 1
	
			#bait2
			Doh2 <- Differential
			ALL_O2 <- findOverlaps(anchorTwo(Doh2),Final_Combined_ATAC)
			ALL_Q2 <- queryHits(ALL_O2)
			ALL_S2 <- subjectHits(ALL_O2)
			Doh2 <- Doh2[ALL_Q2]
			Doh2$ATAC <- ALL_S2
			Doh2$Differential <- ALL_Q2
			Doh2$End <- 2
	
		#annotate Features onto DohALL
		up_regulated <- UP$Symbol
		up_regulated <- up_regulated[!is.na(up_regulated)]
		down_regulated <- DOWN$Symbol
		down_regulated <- down_regulated[!is.na(down_regulated)]

		#Now we have all the ATAC-seq peaks (some duplicated) that are involved in any type of differential Interaction	
		Doh_ALL <- c(Doh1,Doh2)

		#add the motif counts to these peaks
		Doh_ALL$NUM <- Final_Combined_ATAC$NUM[match(Doh_ALL$ATAC, 1:length(Final_Combined_ATAC))]
		Doh_ALL$DAT_where <- Final_Combined_ATAC$DAT_where[match(Doh_ALL$ATAC, 1:length(Final_Combined_ATAC))]
		Doh_ALL$DAT_list <- Final_Combined_ATAC$DAT_list[match(Doh_ALL$ATAC, 1:length(Final_Combined_ATAC))]

		#annotate with the regulatory outcome
		Doh_ALL$UP_DOWN_NoCHANGE <- 'NoChange'
		Doh_ALL[which(Doh_ALL$bait1 %in% up_regulated),]$UP_DOWN_NoCHANGE <- "Up"
		Doh_ALL[which(Doh_ALL$bait1 %in% down_regulated),]$UP_DOWN_NoCHANGE <- "Down"
		Doh_ALL[which(Doh_ALL$bait2 %in% up_regulated),]$UP_DOWN_NoCHANGE <- "Up"
		Doh_ALL[which(Doh_ALL$bait2 %in% down_regulated),]$UP_DOWN_NoCHANGE <- "Down"
		Doh_ALL$Feature <- paste0(Doh_ALL$dir,Doh_ALL$UP_DOWN_NoCHANGE)

		save(Doh_ALL,file="Doh_ALL_Metadata_for_Motifs.Rdata")
		save(Differential,file="Differential_Metadata_for_Motifs.Rdata")
		#these contain the barcoded info used for diwnstream interpretation 

		#convert the data into a data.frame containing the motif counts and feature for each ATAC-seq fragment
		dd  <-  t(as.data.frame(matrix(unlist(Doh_ALL$NUM), nrow=length(unlist(Doh_ALL$NUM[1])))))
		dd <- as.data.frame(dd)
		colnames(dd)[1:length(NAMES)] <- NAMES
		row.names(dd) <- 1:length(Doh_ALL)

		#Add annotations using RNA-seq data
		dd$GAIN_LOSS <- Doh_ALL$dir
		dd$UP_DOWN_NoCHANGE <- 'NoChange'
		up_regulated <- UP$Symbol
		up_regulated <- up_regulated[!is.na(up_regulated)]
		down_regulated <- DOWN$Symbol
		down_regulated <- down_regulated[!is.na(down_regulated)]

		dd[which(Doh_ALL$bait1 %in% up_regulated),]$UP_DOWN_NoCHANGE <- "Up"
		dd[which(Doh_ALL$bait1 %in% down_regulated),]$UP_DOWN_NoCHANGE <- "Down"
		dd[which(Doh_ALL$bait2 %in% up_regulated),]$UP_DOWN_NoCHANGE <- "Up"
		dd[which(Doh_ALL$bait2 %in% down_regulated),]$UP_DOWN_NoCHANGE <- "Down"

		#bait2bait,baittono, notobait?
		dd$B2B <- "NA"
		dd$B2B[which(is.na(Doh_ALL$bait1) == 'FALSE' & is.na(Doh_ALL$bait2) == 'FALSE')] <- 'B2B'
		dd$B2B[which(is.na(Doh_ALL$bait1) == 'TRUE' & is.na(Doh_ALL$bait2) == 'FALSE')] <- 'N2B'
		dd$B2B[which(is.na(Doh_ALL$bait1) == 'FALSE' & is.na(Doh_ALL$bait2) == 'TRUE')] <- 'B2N'
		#Pr2PR,baittono, notobait?
		dd$PR2PR <- "NA"
		dd$PR2PR[which(Doh_ALL$PR1 != 'NA' & Doh_ALL$PR2 != 'NA')] <- 'PR2PR'
		dd$PR2PR[which(Doh_ALL$PR1 == 'NA' & Doh_ALL$PR2 != 'NA')] <- 'N2PR'
		dd$PR2PR[which(Doh_ALL$PR1 != 'NA' & Doh_ALL$PR2 == 'NA')] <- 'PR2N'

		#create feature names
		dd_feature <- paste0(dd$GAIN_LOSS,dd$UP_DOWN_NoCHANGE)
		dd$Feature <- dd_feature
		#Ensure Interaction number is added as a column to help retain 
		#interaction pairs (and remember which differential interaction corresponds o which peak)
		dd$Interaction_Number <- Doh_ALL$Differential
		dd$End <- Doh_ALL$End

		#write out the final 436(+n)xn dataset the motif counts for each ATAC-seq fragment
		write.csv(dd,"New_Non_Biased_ATAC_df_full_strength_B2B_Motif_occurence_differential_Interactions_Double_with_Feature.csv")


################################################################
# Start XGBoost
# Xcategory vs Equal mix of the others (Balanced)
################################################################
# parameter tuning principles from: https://www.hackerearth.com/practice/machine-learning/machine-learning-algorithms/beginners-tutorial-on-xgboost-parameter-tuning-r/tutorial/
# or can just use default parameters (default used below)
# 1:6 refers to each of the 6 categories we want to test:
# DECDown, DECNoChange, DECUp, INCDown, INCNoChange, INCUp


		confusionTEST1 <- list()
		confusionTEST2 <- list()
					library(parallel)
					library(parallelMap) 
					parallelStartSocket(cpus = 4)
					
		cc <- list()
		date="01.21"

		#read back in the processed motif count data
		my_dat <- read.csv("New_Non_Biased_ATAC_df_full_strength_B2B_Motif_occurence_differential_Interactions_Double_with_Feature.csv",head = T,row.names = 1)

		iWant <- names(table(my_dat$Feature))
		cat(iWant,"\n")

		AUC_OUT_Default <- list()	
		confusion_Default <- list()			

		for(k in c(1:6)){	

					set.seed(717)
					iWant <- names(table(my_dat$Feature))
					POSITIVE <- iWant[k]
					NEGATIVE <- iWant[which(iWant != POSITIVE)]
					toTest <- paste0(POSITIVE,"_vs_Mix")

					#ensure no duplicate points
					my_dat <- as.data.frame(my_dat)
					dat <- unique(my_dat[,c(1:436,441)])
					dat <- as.data.frame(dat)
					#dat$Feature <- my_dat$Feature[match(row.names(dat),row.names(my_dat))]
					dat$Interaction_Number <- my_dat$Interaction_Number[match(row.names(dat),row.names(my_dat))]

					# proportion of mix
					# 2 and 5 represent the nochange categories which are much larger than the others so 
					# require a larger mix
		
		
				if(k == 2 || k == 5){
			
						hmm <- dat[which(dat$Feature != POSITIVE),]
						yes_hmm <- dat[which(dat$Feature == POSITIVE),]
						hmm <- hmm[sample(length(hmm[,1]),dim(yes_hmm)[1]),]
						dat <- rbind(yes_hmm,hmm)
						dat$Feature <- as.character(dat$Feature)

						}else{

							hmm <- dat[which(dat$Feature %in% NEGATIVE),]
							yes_hmm <- dat[which(dat$Feature == POSITIVE),]
							PROP <- table(hmm$Feature)/sum(table(hmm$Feature))
							mix1c <- list()		
							for(i in 1:5){
							mix1a <- iWant[which(iWant != POSITIVE)][i]
							mix1b <- hmm[which(hmm$Feature ==mix1a),]
							mix1c[[i]] <- mix1b[sample(length(mix1b[,1]),(dim(yes_hmm)[1]*0.2)),]
						}
		
					hmm <- do.call("rbind",mix1c)
					dat <- rbind(yes_hmm,hmm)

					dat$Feature <- as.character(dat$Feature)
				}

					dat[,1:436] <- apply(dat[,1:436],2,function(x){log10(as.numeric(x) + 0.5)})
					dat$Feature[which(dat$Feature == POSITIVE)] <- 1
					dat$Feature[which(dat$Feature != 1)] <- 0


					dat$Feature <- as.numeric(dat$Feature)
					data_variables <- dat[,-438]		

					#split data into 75%
					train_index <- as.numeric(sample(1:nrow(data_variables), nrow(data_variables)*0.75))

					dat <- as.matrix(dat[,-c(438)])		
					data_variables <- as.matrix(data_variables[,-437])
					#data_label <- dat[,"Feature"][sample(length(dat[,"Feature"]))]
					data_label <- dat[,437]
					data_matrix <- xgb.DMatrix(data = data_variables, label = data_label,missing=0)
					# split train data and make xgb.DMatrix
					train_data   <- data_variables[train_index,]
					train_label  <- data_label[train_index]
					train_matrix <- xgb.DMatrix(data = train_data, label = train_label,missing=0)
					# split test data and make xgb.DMatrix
					test_data  <- data_variables[-train_index,]
					test_label <- data_label[-train_index]
					test_matrix <- xgb.DMatrix(data = test_data, label = test_label,missing=0)

					#create tasks
					# split train data and make xgb.DMatrix
					train_data   <- data_variables[train_index,]
					# split test data and make xgb.DMatrix
					test_data  <- data_variables[-train_index,]

					train_data <- as.data.frame(train_data)
					train_data$Feature <- factor(train_label)
					test_data <- as.data.frame(test_data)
					test_data$Feature <- factor(test_label)

					library(mlr)
					traintask <- makeClassifTask (data = train_data,target = "Feature",positive = "1")
					testtask <- makeClassifTask (data = test_data,target = "Feature",positive = "1")

					#create learner
					#needs to produce probablities
					library(mlr)
					lrn <- makeLearner("classif.xgboost",predict.type = "prob")
					lrn$par.vals <- list(booster= "gbtree", objective="binary:logistic", eval_metric="auc", nrounds=100L, eta=0.1)


					#can run using mlr for best parameter selection, or just use default parameters
					params1 <- makeParamSet( 
					
						# number of splits in each tree
						makeDiscreteParam("max_depth",values = c(1:15)), 
						makeDiscreteParam("min_child_weight",values = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10)), 
						makeDiscreteParam("subsample",values = 0.8), 
						makeDiscreteParam("colsample_bytree", values = 0.8),
						# The number of trees in the model (each one built sequentially)
						makeDiscreteParam("nrounds", values = 100),
						# "shrinkage" - prevents overfitting
						makeDiscreteParam("eta",values = 0.1),
						# L2 regularization - prevents overfitting
						makeDiscreteParam("gamma", values=0)		
						#makeNumericParam("scale_pos_weight",lower = 1,upper = 100)
						#from https://www.kaggle.com/dkyleward/xgboost-using-mlr-package-r

					)
			
					paramsDefault <- makeParamSet( 
						# number of splits in each tree
						makeDiscreteParam("max_depth",values = 5), 
						makeDiscreteParam("min_child_weight",values = 1), 
						makeDiscreteParam("subsample",values = 0.8), 
						makeDiscreteParam("colsample_bytree", values = 0.8),
						# The number of trees in the model (each one built sequentially)
						makeDiscreteParam("nrounds", values = 100),
						# "shrinkage" - prevents overfitting
						makeDiscreteParam("eta",values = 0.1),
						# L2 regularization - prevents overfitting
						makeDiscreteParam("gamma", values=0)		
						#from https://www.kaggle.com/dkyleward/xgboost-using-mlr-package-r

					)

					#set resampling strategy
					#how the data is split and each model is tested? Prformance is averaged across all 10
					rdesc <- makeResampleDesc("CV",stratify = T,iters=5L)

					#search strategy
					ctrl = makeTuneControlGrid()
					#set parallel backend

					#parameter tuning
					#change par.set to default or parameters you want to tune
					ptm <- proc.time()
					mytune <- tuneParams(learner = lrn, task = traintask, resampling = rdesc, measures = auc, par.set = paramsDefault, control = ctrl, show.info = T)
					mytune$y 
					#0.9526276
					cat(proc.time() - ptm , "\n")

					#set hyperparameters
					lrn_tune <- setHyperPars(lrn,par.vals = mytune$x)

					#train model
					xgmodel <- train(learner = lrn_tune,task = traintask)

					#predict model
					xgpred <- predict(xgmodel,testtask)
					confusion_Default[[k]] <- confusionMatrix(xgpred$data$response,xgpred$data$truth,positive = "1")
					AUC_OUT_Default[[k]] <- mytune$y	
					names(AUC_OUT_Default) <- POSITIVE	


					write.csv(dat,paste0("Data_Used_For_",toTest,"_NEW_params_date_",date,".csv"))
					CD <- confusion_Default
					save(CD,file=paste0(toTest,"_Confusion_date_",date,".Rdata"))
					save(xgmodel,file=paste0(toTest,"_Final_Model_date_",date,".Rdata"))
					AUCO <- AUC_OUT_Default
					save(AUCO,file=paste0(toTest,"_AUC_Final_Model_date_",date,".Rdata"))
	
		}

	# The above produces 4 files:
	# 1 = Data_Used_For_x.csv has the original testing dataset used to generate the AUC scores reported
	# 2 = _Confusion_date_.Rdata is an object containing the final statistics for each of the models predictive scores
	# 3 = _Final_Model_date_.Rdata contains the parameters used for the initial model
	# 4 = _AUC_Final_Model_date_.Rdata contains the AUC score for the initial model

		#One vs All Scores
		#plot the AUC scores for each model as a barplot
		SCORES <- list.files(pattern="_AUC_Final_Model_date_",full.names=T)
		SCORES_in <- lapply(SCORES,function(x){load(x)
			SCORES_in <- AUCO})

		png(paste0("./",graphs.dir,"ML_Scores_AUC.png"), width = 2.5, height = 10, units = 'in', res = 600)
		barplot(c(as.numeric(SCORES_in[[1]][1]),
				as.numeric(SCORES_in[[2]][2]),
				as.numeric(SCORES_in[[3]][3]),
				as.numeric(SCORES_in[[4]][4]),
				as.numeric(SCORES_in[[5]][5]),
				as.numeric(SCORES_in[[6]][6])),col=c("pink","indianred","purple","blue","lightseagreen","green"),ylim=c(0,1),cex.axis=2)
		dev.off()
		
		#the above plot can be compared to the "ML_Scores_AUC_Expected.png" plot in the graphs directory


################################################################
# Next step is to use the above models and split the dataset into 5 equal groupings,
# perform a prediction and generate SHAP scores encompassing 100% of the data for downstream analysis
################################################################

		Model <- list()
		brr <- list.files(pattern="_AUC_Final_Model_date_",full.names=T)
		for(i in 1:6){
		load(brr[[i]])
			Model[[i]] <- xgmodel
		}

		set.seed(717)
		
		CC <- list()
		shap_data_sub_OUT <- list()
		dd <- list()		
		for(k in c(1:6)){
			NUMBER=k
			
			###		
			###Define the feature you want to test (need to assign this as 1 being the positve feature)		
			###		
				set.seed(717)
				iWant <- names(table(my_dat$Feature))
				POSITIVE <- iWant[k]
				NEGATIVE <- iWant[which(iWant != POSITIVE)]
				toTest <- paste0(POSITIVE,"_vs_Mix")
		
					my_dat <- as.data.frame(my_dat)
					dat <- unique(my_dat[,c(1:436,441)])
					dat <- as.data.frame(dat)
					#dat$Feature <- my_dat$Feature[match(row.names(dat),row.names(my_dat))]
					dat$Interaction_Number <- my_dat$Interaction_Number[match(row.names(dat),row.names(my_dat))]
		
				if(k == 2 || k == 5){
			
					hmm <- dat[which(dat$Feature != POSITIVE),]
					yes_hmm <- dat[which(dat$Feature == POSITIVE),]
					hmm <- hmm[sample(length(hmm[,1]),dim(yes_hmm)[1]),]
					dat <- rbind(yes_hmm,hmm)
					dat$Feature <- as.character(dat$Feature)

				}else{
		
					#proportion of mix
					hmm <- dat[which(dat$Feature %in% NEGATIVE),]
					yes_hmm <- dat[which(dat$Feature == POSITIVE),]
					PROP <- table(hmm$Feature)/sum(table(hmm$Feature))
					mix1c <- list()		
					for(i in 1:5){
						mix1a <- iWant[which(iWant != POSITIVE)][i]
						mix1b <- hmm[which(hmm$Feature ==mix1a),]
						mix1c[[i]] <- mix1b[sample(length(mix1b[,1]),dim(yes_hmm)[1]*0.2),]
					}
					
					hmm <- do.call("rbind",mix1c)
					dat <- rbind(yes_hmm,hmm)
					dat$Feature <- as.character(dat$Feature)

				}

					dat[,1:436] <- apply(dat[,1:436],2,function(x){log10(as.numeric(x) + 0.5)})
					dat$Feature[which(dat$Feature == POSITIVE)] <- 1
					dat$Feature[which(dat$Feature != 1)] <- 0
			
					dat$Feature <- as.numeric(dat$Feature)
					# Full data set
					data_variables <- dat[,-c(438)]		
					original_train_index <- as.numeric(sample(1:nrow(data_variables), nrow(data_variables)*0.75))
								
					#split data into 100% but randomly shuffled and keeping pairs
					ALL_train_index <- as.numeric(sample(1:nrow(data_variables), nrow(data_variables)*1))
		n = 5
			RANDY <- split(ALL_train_index, sort(ALL_train_index%%n))
			RANDY_unlist <- unlist(RANDY)
			data_variables <- data_variables[RANDY_unlist,]
					#all RAN in order
				traintask <- list()
				testtask <- list()
				CC <- list()
				xgpred <- list()
				xgmodel <- list()
				final_df_SHAP_F <- list()
		

		for(i in 1:5){

					lrn_tune <- setHyperPars(Model[[grep(toTest,brr)]]$learner)

						NUMBER <- k
		
						train_index <- RANDY_unlist[!(RANDY_unlist %in% RANDY[[i]])]
						data_label <- as.numeric(data_variables[,"Feature"])

						train_data   <- data_variables[train_index,]
						train_label  <- data_label[train_index]

						test_data  <- data_variables[-train_index,]
						test_label <- data_label[-train_index]
												
								train_data <- as.data.frame(train_data)
								train_data$Feature <- factor(train_label)
								test_data <- as.data.frame(test_data)
								test_data$Feature <- factor(test_label)

						traintask[[i]] <- makeClassifTask (data = train_data,target = "Feature",positive = "1")
						testtask[[i]] <- makeClassifTask (data = test_data,target = "Feature",positive = "1")
				

					#train model

					xgmodel[[i]] <- train(learner = lrn_tune,task = traintask[[i]])
					xgpred[[i]] <- predict(xgmodel[[i]] ,testtask[[i]])
			
					PRED <- xgpred[[i]] 
					CC[[i]] <- confusionMatrix(PRED$data$response,PRED$data$truth,positive = "1")
		#max(xgmodel$evaluation_log$train_auc)
				dd[[k]] <- CC



					###
					#SHAP me UP
					###

						PREDICTION <- data.frame(truth = as.numeric(PRED$data$truth), response = as.numeric(PRED$data$response))
						PREDICTION$truth[which(PREDICTION$truth == 2)] <- toTest
						PREDICTION$truth[which(PREDICTION$truth == 1)] <- "Other"
						PREDICTION$response[which(PREDICTION$response == 2)] <- toTest
						PREDICTION$response[which(PREDICTION$response == 1)] <- "Other"

						final_model <- xgmodel[[i]]$learner.model
						testing_data <- as.matrix(test_data[,1:436])
						shap_values <- shap.values(xgb_model = final_model, X_train = testing_data)

						shap_data <- shap_values$shap_score
						shap_data_sub <- shap_data
						shap_data_sub$predicted <- PREDICTION$response
						shap_data_sub$actual <- PREDICTION$truth
						shap_data_sub$InteractionRow <- row.names(PRED$data)
                   
                    final_df_SHAP <- test_data
					final_df_SHAP$iteration <- i
					final_df_SHAP_F[[i]] <- final_df_SHAP
					shap_data_sub$iteration <- i
					shap_data_sub_OUT[[i]] <- shap_data_sub
				}

				final_df_SHAP_ALL <- do.call("rbind",final_df_SHAP_F)
				shap_data_sub_OUT_ALL <- do.call("rbind",shap_data_sub_OUT)
				#AUC_OUT_Final <- mytune$y	

				write.csv(final_df_SHAP_ALL,paste0("my_data_raw_for_",toTest,".csv"))
				write.csv(shap_data_sub_OUT_ALL,paste0("my_data_Fixed_SHAP_for_",toTest,".csv"))
				dd[[k]] <- CC
				save(dd,file=paste0("final_confusion_matrix_for_",toTest,".Rdata"))
				save(xgmodel,file=paste0("final_model_parameters_for_",toTest,".Rdata"))

}

# Only correctly predicted datapoints are taken forward for analysis
# The above should produce 4 main outputs for each of the regulatory categories
# 1 = my_data_raw_for_x.csv is the final combined dataset used for predictions
# 2 = my_data_Fixed_SHAP_for_.csv is the final SHAP scores for the combined dataset used for predictions
# The first 1:436 columns should be paired across the 2 files, illustrating the shap score for that particular datapoint
# 3 = final_confusion_matrix_for_x.Rdata is an object containing the final statistics for each of the models predictive scores
# this is split into 5 showing each of the iterations to combine the whole dataset
# 4 = the final xgboost model parameters used for each iteration





