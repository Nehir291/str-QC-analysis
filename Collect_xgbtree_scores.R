library(caret)
library(reshape2)
library("purrr")
library(party)
library(grid)
library("patchwork")
library(dplyr)
library("tidyr")

#Upload calls
#Calls <- read.delim("/Users/enk15/Downloads/scripts_previous/HG00514.1_json_files.variants_reduced.tsv",header=TRUE)
#metrices <- read.table("/Users/enk15/Downloads/scripts_previous/HG00514_metrics.csv",header=TRUE)
#name <- "HG00514"

args = commandArgs(TRUE)
Calls <- read.delim(args[1])
metrices <- read.table(args[2],header=TRUE)
name <- args[3]


index<- which(names(Calls) == "Repeat.Size..bp...Allele.1")
colnames(Calls)[index[1]] <- "Genotype_size_short"

index<- which(names(Calls) == "Repeat.Size..bp...Allele.2")
colnames(Calls)[index[1]] <- "Genotype_size_long"

#Upload image metrices
merged_data <- merge(metrices, Calls, by.x = "VariantId",by.y = "VariantId")
merged_data <-  merged_data %>% separate(AlleleDepth, c('Depth_short', 'Depth_long'),sep = "/")
merged_data <-  merged_data %>% separate(Genotype.x, c('Genotype_short', 'Genotype_long'),sep = "/")

merged_data$Genotype_short <- as.numeric(merged_data$Genotype_short)
merged_data$Genotype_long <- as.numeric(merged_data$Genotype_long)
merged_data$Depth_short <- round(as.numeric(merged_data$Depth_short), digits = 0)
merged_data$Depth_long <- round(as.numeric(merged_data$Depth_long), digits = 0)


#Collect Q score
merged_data$Q_score_hap1 <- lapply(1:nrow(merged_data),function(i) (1/exp(4*(as.numeric(merged_data$CI.size..Allele.1[i])/as.numeric(merged_data$Genotype_short[i])))))

merged_data$Q_score_hap2 <- lapply(1:nrow(merged_data),function(i) (1/exp(4*(as.numeric(merged_data$CI.size..Allele.2[i])/as.numeric(merged_data$Genotype_long[i])))))

#Collect Interruption score by genotype
merged_data$Interruption_score_hap1_genotype <- lapply(1:nrow(merged_data),function(i) (1/exp(4*(as.numeric(merged_data$Total_upper_Orange[i])/as.numeric(merged_data$Genotype_short[i])))))
merged_data$Interruption_score_hap2_genotype <- lapply(1:nrow(merged_data),function(i) (1/exp(4*(as.numeric(merged_data$Total_bottom_Orange[i])/as.numeric(merged_data$Genotype_long[i])))))

#Collect interruption score by depth
merged_data$Interruption_score_hap1_depth <- lapply(1:nrow(merged_data),function(i) (1/exp(4*(as.numeric(merged_data$Total_upper_Orange[i])/as.numeric(merged_data$Depth_short[i])))))
merged_data$Interruption_score_hap2_depth <- lapply(1:nrow(merged_data),function(i) (1/exp(4*(as.numeric(merged_data$Total_bottom_Orange[i])/as.numeric(merged_data$Depth_long[i])))))



# Comparison of each metrices between TRUE and FALSE
merged_data$Genotype_short <- as.numeric(merged_data$Genotype_short)
merged_data$Genotype_long <- as.numeric(merged_data$Genotype_long)
merged_data$Depth_short <- round(as.numeric(merged_data$Depth_short), digits = 0)
merged_data$Depth_long <- round(as.numeric(merged_data$Depth_long), digits = 0)
merged_data$Interruption_score_hap1_depth <- round(as.numeric(merged_data$Interruption_score_hap1_depth), digits = 2)
merged_data$Interruption_score_hap2_depth<- round(as.numeric(merged_data$Interruption_score_hap2_depth), digits = 2)
merged_data$Interruption_score_hap1_genotype <- round(as.numeric(merged_data$Interruption_score_hap1_genotype) ,digits = 2)
merged_data$Interruption_score_hap2_genotype <- round(as.numeric(merged_data$Interruption_score_hap2_genotype),digits = 2)
merged_data$Q_score_hap1<- round(as.numeric(merged_data$Q_score_hap1), digits = 2)
merged_data$Q_score_hap2<- round(as.numeric(merged_data$Q_score_hap2), digits = 2)

merged_data$Path <-  paste(merged_data$SampleId,"_",merged_data$VariantId,".svg",sep="")
merged_data$Path <-  sub("-","_",merged_data$Path)
merged_data$Path <-  sub("-","_",merged_data$Path)
merged_data$Path <-  sub("-","_",merged_data$Path)


#To test mL tool, first prepare the input data
merged_data_red <- merged_data %>% select(VariantId,CI.size..Allele.1,CI.size..Allele.2                           
                                            , Genotype_short,Genotype_long ,Depth_short,Depth_long, Times_Upper_Interruptions,Count_Upper_reads_with_interruptions       
                                            ,Times_Upper_Gaps,Times_Bottom_Interruptions,Count_Bottom_reads_with_interruptions     
                                            ,Times_Bottom_Gaps,Total_interrupting_nucleotides_upper,Median_upper,Total_interrupting_nucleotides_bottom      
                                            ,Median_bottom,Total_upper_Orange,Median_orange_upper                      
                                            ,Total_read_upper,Count_Orange_Upper_reads_with_interruptions,Total_bottom_Orange,Median_orange_bottom                       
                                            ,Total_read_bottom,Count_Orange_Bottom_reads_with_interruptions,SampleId,                                       
                                            RepeatUnitLength,Genotype_size_short,Genotype_size_long,Q_score_hap1,Q_score_hap2,                               
                                            Interruption_score_hap1_genotype,Interruption_score_hap2_genotype,Interruption_score_hap1_depth,Interruption_score_hap2_depth,Path)   


#Method 1: full set of Pacbio data + subset of BW data with new metrices Q score, interruption score

fit.xgbTree_manually_selected <- readRDS("fit.xgbTree_PB_manually_selected.rda")
predictions_xgbTree_manually_selected <- predict(fit.xgbTree_manually_selected, merged_data_red, type="prob")
merged_data_red_ID <- merged_data_red %>% select (Path)
predictions_xgbTree_manually_selected_scored <- cbind(merged_data_red_ID,predictions_xgbTree_manually_selected)
colnames(predictions_xgbTree_manually_selected_scored)<- c("Path","xgbTree_True","xgbTree_False")
predictions_v1 <- merge(predictions_xgbTree_manually_selected_scored,merged_data,by.x="Path",by.y="Path" )

predictions_v2<-  predictions_v1 %>% select (Path,xgbTree_True)

setwd("/Users/enk15/Downloads/scripts_previous")
write.table(predictions_v2, file = paste0(name,"_xgbtree_scores.tsv"), sep = "\t",
            row.names = FALSE,col.names = TRUE,quote=FALSE)
