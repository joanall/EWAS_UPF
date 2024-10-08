# Title: Do Epigenome Wide Association Analysis 
# Author:  JLLP
# Date: 17/02/2023

# Description: Run a RLM model to asses association between UPF and DNAm where,
# the output are b values winsorized and the main exposure is UPF.  


#==============================================================================#
#                           DEFINE USER PARAMETERS                             #
#==============================================================================#

# Define directory to save the outputs & create it 

output_directory <- "results/Results_Models_HELIX/"

# Note: Working directory is root folder -> use absolutes paths
# Check if directory to save results exists, if not create one 

if (!dir.exists(output_directory)) { dir.create(output_directory)
} else { print("Directory already exists")}


# Define inout path for metadata and mehtylation data  
metadata_directory <- "db/final_metadata.RDS"
methylation_directory <- "db/final_bVals.RDS"


# Define ethnicity (depend on each cohort, write largest ethnicity group)
ethnicity <- "caucasian" #in HELIX: all or caucasian

var_name_ethnicity <- 'child_caucasian'
category_ethnicity <- 'yes'



# Define model -> from shell 
# Note: When running bash command we specifiy in --array which models we want 
# to use. By default uses last element of array (i.e. model.4.CT)

args <- commandArgs(trailingOnly = TRUE)

print(paste0("Input model:  ", args[1])) 

# sbatch --array=1-8 Call_Run_EWAS.sh 




#==============================================================================#
#                                                                              #
#                           SET WORKING ENVIRONMENT                            #----
#                                                                              #
#==============================================================================#
# rm(list=ls())
# gc()

# Load libraries 

library(tidyverse)
library(dplyr)
library(sfsmisc)
library(parallel)


#------------------------------------------------------------------------------#
#                           Define models                                      #----
#------------------------------------------------------------------------------#

exposure <- 'propUPF_scaled'

model.1 <- c("cohort","child_sex", "child_age") # basic "child_caucasian"

model.2 <- c(model.1, "maternal_age", "maternal_bmi", "maternal_smoking",
                            "maternal_education")

model.3 <- c(model.2,
             "child_bmi","child_sedentary_behaviour") # child related vars (+ ethnicity)

model.4 <- c(model.3,
             "total_vegetables_intake","total_fruit_intake") # diet related vars


cells <-  c("NK_6", "Bcell_6", "CD4T_6", "CD8T_6", "Gran_6", "Mono_6")

model.1.CT <- c(model.1, cells)

model.2.CT <- c(model.2, cells)

model.3.CT <- c(model.3, cells)

model.4.CT <- c(model.4, cells)

which_model <- eval(parse(text=args[1]))
print(which_model)



#------------------------------------------------------------------------------#
#                          Define Functions                                    #----
#------------------------------------------------------------------------------#

# Define function to take care of extrem methylation values 
winsorize_row <- function(data, pct = winsorize.pct) {
  
  stopifnot(is.matrix(data))
  
  # Estimates quantiles for each row (column) in a matrix
  quantiles <- matrixStats::rowQuantiles(data, probs = c(pct,1-pct), na.rm = T) 
  
  # Define limit values in numeric named vector 
  
  low <- quantiles[,1]
  upper <- quantiles[,2]
  
  # contar els samples que estan per sota/damunt del limit? -> sempre 6? 
  
  outliers.lower <- rowSums(data < low, na.rm=T)
  outliers.upper <- rowSums(data > upper, na.rm=T)
  
  # For those CpGS with an outline value, get their index/position in methylation data 
  idx <- which(data < low, arr.ind=T)
  
  # Replace outlines for last inline value
  data[idx] <- low[idx[,1]]
  
  idx <- which(data > upper, arr.ind=T)
  data[idx] <- upper[idx[,1]]
  
  # Check there aren't NA values 
  n <- rowSums(!is.na(data))
  
  log <- data.frame(outliers.lower, outliers.upper, n)
  # 
  return(list(data=data, log=log))
}



# Define fucntion to run rlm 
#  We need to specify:
#   - formula (outcome + exposure + covariates)
#   - data (meth values + variables vals)


rlm.single.cpg <- function(cpg, exposure, covariates, metadata, meth_vals) {
  
  exp_cov <- c(exposure, covariates)
  
  df_exp_cov <- metadata %>% dplyr::select(all_of(exp_cov))
  
  # Formula of robust linear regression
  f <- as.formula(paste0(cpg, "~", paste(exp_cov, collapse=" + ")))
  

  # Fit model with our data
  fitted_model <- MASS::rlm(formula = f, data=cbind(meth_vals, df_exp_cov))
  

  # P VALUES
  p.value <- as.numeric(f.robftest(fitted_model, var = exposure)["p.value"])

  result <- cbind (as.data.frame(summary(fitted_model)$coefficients)[c(exposure),], p.value)
  
  rownames(result) <- cpg
  result <- rename(result, c('estimate' = 'Value' ,
                             'std.error' = 'Std. Error',
                             't.value' = 't value'))
  
  return(result)
}

#===============================================================================


#==============================================================================#
#                                                                              #
#                               LOAD & ARRANGE DATA                            #----
#                                                                              #
#==============================================================================#


# Load metadata with co-variables 
metadata <- readRDS(file = metadata_directory)

# Load methylation data
meth_vals <- readRDS(file = methylation_directory)


# Take care of methylation outliers (do winsorization) if B values 
  meth_vals_copy <- as.data.frame(meth_vals)
replace.outliers <- winsorize_row(t(as.matrix(meth_vals_copy)), 0.005)

# From winsorization results get dataframe 
meth_vals <- as.data.frame(t(as.matrix(replace.outliers$data)))

# Check order 
all(colnames(meth_vals) == colnames(meth_vals_copy))  #TRUE
all(rownames(meth_vals) == rownames(meth_vals_copy))  #TRUE

# Save winsorized data
save(meth_vals, file= paste0(output_directory, "/3_meth_vals_winsorized.RData"))

# Save winsorization log
wins_log <- as.data.frame(replace.outliers$log)

write.table(wins_log, file = paste0(output_directory,"/winslog",Sys.Date(),".txt"),
            sep = "\t", col.names = T, row.names = T, append = F, quote=FALSE)

rm(replace.outliers)
rm(meth_vals_copy)
  

# Include only largest ethnicity
if (exists("var_name_ethnicity") && var_name_ethnicity != "") {
  
  metadata_copy <- data.frame(metadata)
  
  metadata <- metadata_copy %>%
    dplyr::filter(!!sym(var_name_ethnicity) == category_ethnicity)  # keep only specified ethnicity
  
  # Remove smaller ethnicity group samples from methylation data 
  meth_vals_copy <- data.frame(meth_vals)
  meth_vals <- meth_vals_copy[rownames(meth_vals_copy) %in% rownames(metadata), ]
  
  print(paste0('Sample size: ', nrow(metadata)))
  print(paste0('Number of CpGs: ', ncol(meth_vals))) 
  
  # Remove copies of original inputs
  rm(metadata_copy)
  rm(meth_vals_copy)
}



#==============================================================================#
#                                                                              #
#                              RUN PARALLELIZATION                             #----
#                                                                              #
#==============================================================================#

cpg <- colnames(meth_vals)

rlm.result <- parallel::mclapply(X = cpg, 
                                FUN = rlm.single.cpg,
                                exposure = exposure,
                                covariates = which_model,
                                metadata = metadata,
                                meth_vals = meth_vals,
                                mc.cores=8)

final_result <- dplyr::bind_rows(rlm.result)

complete_result <- final_result %>%
  dplyr::mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(p.value) #order ascending by p value

# Save EWAS results
if (exists("var_name_ethnicity") && var_name_ethnicity != "") {
  save_to <- paste0(output_directory,"rlm.result.ethn.", args[1],".rds")
  
} else {
  save_to <- paste0(output_directory,"rlm.result.", args[1],".rds")
  
  
}

write_rds(complete_result, file = save_to)





