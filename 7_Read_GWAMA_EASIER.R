# Author:  JLLP
# Title: Read ewas resutls after met analysis with gwama 
# Date: 13/01/2023



# Description: Read results of EWAS. Check its quality by ploting/computing 
# lambda and some plots. Also, get significant or top CpGs and plot them. 

## GWAMA Output Format

# 1. chromosome - Marker chromosome 
# 2. position - Marker position (bp) 
# 3. rs_number - Marker ID 
# 4. reference_allele - Effect allele 
# 5. other_allele - Non effect allele 
# 6. eaf - effect allele frequency as reported from the 
# UK10K + 1000 Genomes reference panel. 
# 7. beta - Overall beta value for meta-analysis 
# 8. se - standard error of beta.
# 9. beta_95L - Lower 95% CI for BETA 
# 10. beta_95U - Upper 95% CI for BETA 
# 11. z - Z-score 
# 12. p-value - Meta-analysis p-value 
# 13. -log10_p-value - Absolute value of logarithm of meta-analysis 
# p-value to the base of 10. 
# 14. q_statistic - Cochran's heterogeneity statistic 
# 15. q_p-value - Cochran's heterogeneity statistic's p-value 
# 16. i2 - Heterogeneity index I2 by Higgins et al 2003 
# 17. n_studies - Number of studies with marker present 
# 18. n_samples - Number of samples with marker present (will be NA 
# 	            if marker is present in any input file where N column 
# 				is not present) 
#     19. effects - Summary of effect directions ('+' - positive effect of 
# 	          reference allele, '-' - negative effect of reference 
# 			  allele, '0' - no effect (or non-significant) effect of 
# 			  reference allele, '?' - missing data)
# 			  
# ANNOTATION TO LOOK AT:  "pos" y "chr"	  
################################################################################
#                      SET WORKING ENVIRONMENT                                 #
################################################################################  

# rm(list=ls())
# gc()

# Load libraries 

library(readr)
library(qqman)
library(dplyr)
library(tidyverse)


################################################################################
#                      DEFINE PARAMETERS                                       #
################################################################################ 

# Read .out from meta-analysis 
    #Fixed/Random.out -> son els ouputs directes de GWAMa
    #Fixed_Modif.iut -> contenent els adjusments 

# Main Models
files <- c('results/GWAMA_Results/Meta_Model_1CT_Filtr/Meta_Model_1CT_Fixed_Modif.out', 
           'results/GWAMA_Results/Meta_Model_1CT_Filtr/Meta_Model_1CT_Random_Modif.out', 
           
           "results/GWAMA_Results/Meta_Model_2CT_Filtr/Meta_Model_2CT_Fixed_Modif.out",
           "results/GWAMA_Results/Meta_Model_2CT_Filtr/Meta_Model_2CT_Random_Modif.out",
           
           "results/GWAMA_Results/Meta_Model_3CT_Filtr/Meta_Model_3CT_Fixed_Modif.out", 
           "results/GWAMA_Results/Meta_Model_3CT_Filtr/Meta_Model_3CT_Random_Modif.out")

# Sensitivity Models           
files <- c("results/GWAMA_Results/Meta_ethnicity_Filtr/Meta_ethnicity_Fixed_Modif.out",
           "results/GWAMA_Results/Meta_ethnicity_Filtr/Meta_ethnicity_Random_Modif.out")
           #"results/GWAMA_Results/Meta_Model_healthy_diet_Filtr/Meta_Model_healthy_diet_Fixed_Modif.out",
           #"results/GWAMA_Results/Meta_Model_healthy_diet_Filtr/Meta_Model_healthy_diet_Random_Modif.out")

# Leave one out Models           
files <- c("results/GWAMA_Results/Meta_Out_ALSPAC_Filtr/Meta_Out_ALSPAC_Fixed_Modif.out",
           "results/GWAMA_Results/Meta_Out_ALSPAC_Filtr/Meta_Out_ALSPAC_Random_Modif.out")

# Sub HELIX Models


outputfolder <- 'results/Read_GWAMA/'

################################################################################
#                      DEFINE FUNCTION                                         #
################################################################################ 


# Function to calculate lambda & confidence intervals
inflation <- function(pvector) {
  chisq <- qchisq(1 - pvector, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq))))
  upper_ci <-  lambda + (SE.median / qchisq(0.5, 1))
  lower_ci <-  lambda - (SE.median / qchisq(0.5, 1))
  lambda_ci <- c(lambda, upper_ci, lower_ci)
  return(lambda_ci)
}

# Function to do a qqplot
do_qqplot <- function(df,var_name){
  
  p <- ggpubr::ggqqplot(df, x = var_name) 
  p <- p + ggtitle(paste0('QQ plot ', var_name))
  return(p)
}

# Function to add lambda and CI intervals to a qqplot
add_lambda_to_plot <- function(pvector, qqplot) {
  lambda <- inflation(pvector)
  # Add lamba + CI 
  lambda_CI <- sprintf("λ = %.2f [%.2f,%.2f]", lambda[1], lambda[2], lambda[3]) 
  qqplot <- qqplot + annotate(geom = "text", x = -Inf, y = Inf,
                              hjust = -0.15, vjust = 1 + 0.15 * 3, label = lambda_CI, size = 5)
  return(qqplot)
}

#------------------------------------------------------------------------------#
#                 GET INFLATION AND TRANSFORM FILES                            #
#------------------------------------------------------------------------------#


for (i in seq_along(files)) {
  
  filename <- gsub('.out', "", basename(files[i]))
  cat('\nReading... ', filename, "\n")
  gwama_result <- read.table(files[i], header = TRUE)
  
  # Save results meta analysis into rds & csv
  cat('Saving results in rds & csv format... ')
  write_rds(gwama_result, file = paste0(outputfolder, filename,'.rds') )
  data.table::fwrite(gwama_result, file = paste0(outputfolder, filename,'.csv'),sep=";")
  
  # Get lambda and CI
  cat('Calculating lambda and CI ')
  lambda_ci <- as.data.frame(t(inflation(gwama_result$p.value)))
  colnames(lambda_ci) <- c('lambda','upper_ci', 'lower_ci')
  
  new_filename <- paste0(outputfolder, 'lambda_',filename, '.csv')
  write_csv(lambda_ci, new_filename)
}

#------------------------------------------------------------------------------#
#                            GET  TOP    CPGs                                  #
#------------------------------------------------------------------------------#

# Get CpGs (if any) at thresholds: Bonferroni (yes/no), FDR (0.05) and suggestive (1e-05)

for (i in seq_along(files)) {
  
  filename <- gsub('.out', "", basename(files[i]))
  cat('\nReading... ', filename, "\n")
  gwama_result <- read.table(files[i], header = TRUE)
  
  top <- subset(gwama_result, gwama_result$Bonferroni == 'yes')
  if (nrow(top) != 0) {write.table(top, paste0(outputfolder, 'bonferroni_',filename, '.csv'), sep=";")}
  
  top <- subset(gwama_result, gwama_result$FDR < 0.05)
  if (nrow(top) != 0){write.table(top, paste0(outputfolder, 'fdr_',filename, '.csv'), sep=";")}  
    
  top <- subset(gwama_result, gwama_result$p.value < 1e-05)
  if (nrow(top) != 0){write.table(top, paste0(outputfolder, 'suggestive_',filename, '.csv'),sep=",")}
}
































# Use write_csv from reader isntead of write.csv, 2x faster
# use fwrite from data.table is faster than write_csv


# # EWAS catalog - retrive asscoaiiton for top cpgs
# library(devtools)
# library(jsonlite)
# # https://rdrr.io/github/jrs95/ewascatalog/
# # install_github("ewascatalog/ewascatalog-r") doesnt work, copied function hehe 
# # @author James R Staley <js16174@bristol.ac.uk>
# library()
# ewascatalog <- function(cpgquery=NULL,regionquery=NULL,genequery=NULL,traitquery=NULL){
#   if(is.null(cpgquery) & is.null(regionquery) & is.null(genequery) & is.null(traitquery)) stop("no query entered")
#   if((length(cpgquery[1])+length(regionquery[1])+length(genequery[1])+length(traitquery[1]))>1) stop("only one query type allowed")
#   if(!is.null(regionquery)){
#     ub <- as.numeric(sub("*.-", "", sub(".*:", "",regionquery)))
#     lb <- as.numeric(sub("-.*", "", sub(".*:", "",regionquery)))
#     dist <- ub - lb
#     if(any(dist>10000000)) stop("region query can be maximum of 10mb in size")
#   }
#   if(!is.null(cpgquery)){
#     query <- paste0("cpgquery=",cpgquery)
#   }
#   if(!is.null(regionquery)){
#     query <- paste0("regionquery=",regionquery)
#   }
#   if(!is.null(genequery)){
#     query <- paste0("genequery=",genequery)
#   }
#   if(!is.null(traitquery)){
#     query <- paste0("traitquery=",sub(" ", "+", tolower(traitquery)))
#   }
#   if(length(query)>100) stop("a maximum of 100 queries can be requested at one time")
#   results <- data.frame()
#   for(i in 1:length(query)){
#     json_file <- paste0("http://www.ewascatalog.org/api/?",query[i])
#     print(json_file)
#     json_data <- fromJSON(file=json_file)
#     if(length(json_data)==0){
#       cat("No results for",sub(".*=","",query),"\n")
#       next
#     }
#     fields <- json_data$fields
#     tables <- as.data.frame(matrix(unlist(json_data$results), ncol=length(fields), byrow=T))
#     names(tables) <- fields
#     results <- rbind(results,tables)
#   }
#   return(results)
# }
# 
# 
# CpG_res <- ewascatalog(cpgquery="cg00029284")



# Region res <- ewascatalog(regionquery="6:15000000-25000000")

# Gene res <- ewascatalog(genequery="FTO")

# Trait res <- ewascatalog(traitquery="body mass index")
