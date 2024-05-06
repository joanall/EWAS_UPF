# Author:  JLLP
# Title: Built dataframe results EWAS
# Date: 03/01/2023

# Description: Read results of EWAS and :
# - 1) Compute the adjusted P values if it's not done, code example as follows
      # df <- df %>%
      #   dplyr::mutate(q.value = p.adjust(p.value, method = "BH")) %>%
      #   dplyr::arrange(p.value) #order ascending by p value

# - 2) Make it readable format to run the Qaulity Control of EASIER package. 

################################################################################
#                      SET WORKING ENVIRONMENT                                 #
################################################################################

library(readr)


# -----------------------------------------------------------------------------#
#                                                                              #
#                                   HELIX                                      #                                         
#                                                                              #        
#------------------------------------------------------------------------------#

# Define input and output directory 
input_directory <- "results/EWAS_RESULTS/HELIX"
output_directory <- "results/preQC/HELIX/" #Note: preQC have to be created before, dir.create just creates one sub-directory not 2 

# Check if directory exists, if not create it
if (!dir.exists(output_directory)) { dir.create(output_directory) } else { print("Directory already exists")}

# Get files 
results_files <- list.files(path=input_directory, pattern = "3_rlm.result.model.[1-4].CT",  full.name=TRUE)
# Add file with main model and only white European ethnicity 
results_files <- c('results/EWAS_RESULTS/HELIX/3_rlm.result.cauc.model.3.CT.B.rds')


for (i in seq_along(results_files)){
  
  cat('Reading file... ', results_files[i], "\n")
  
  rlm_result <- readRDS(results_files[i]) 
  colnames(rlm_result) <- c('estimate','std.error', 't.value', 'p.value', 'q.value')

  #Renaming of variables to match 
  rlm_result$probeID <- rownames(rlm_result)
  
  names(rlm_result)[names(rlm_result) == 'estimate'] <- 'BETA'
  names(rlm_result)[names(rlm_result) == 'std.error'] <- 'SE'
  names(rlm_result)[names(rlm_result) == 'p.value'] <-  'P_VAL'
  
  rlm_result <- rlm_result[,c('BETA','SE', 't.value', 'P_VAL', 'q.value','probeID')]
  # Save
  new_filename <- basename(gsub('.rds', '.txt',results_files[i]))
  write.table(rlm_result, paste0(output_directory,new_filename), 
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  readr::write_rds(rlm_result, file = paste0(output_directory,basename(results_files[i])))
  
}


# -----------------------------------------------------------------------------#
#                                                                              #
#                             Generation XXI (GXXI)                            #                                         
#                                                                              #        
#------------------------------------------------------------------------------#

# Define input and output directory 
input_directory <- "results/EWAS_RESULTS/GXXI"
output_directory <- "results/preQC/GXXI/" #Note: preQC have to be created before, dir.create just creates one sub-directory not 2 


# Check if directory exists, if not create it
if (!dir.exists(output_directory)) { dir.create(output_directory) } else { print("Directory already exists")}

# Get files 
results_files <- list.files(path=input_directory, pattern = "3_rlm.result.model.[1-4].CT",  full.name=TRUE)
# Add file with main model and only white European ethnicity 
results_files <- c('results/EWAS_RESULTS/GXXI/3_rlm.result.model.3.CT.B.rds')

for (i in seq_along(results_files)) {
  
  cat('Reading file... ', results_files[i], "\n")
  
  rlm_result <- readRDS(results_files[i]) 
  
  colnames(rlm_result) <- c('estimate','std.error', 't.value', 'p.value', 'q.value')
  
  #Renaming of variables to match 
  rlm_result$probeID <- rownames(rlm_result)
  
  names(rlm_result)[names(rlm_result) == 'estimate'] <- 'BETA'
  names(rlm_result)[names(rlm_result) == 'std.error'] <- 'SE'
  names(rlm_result)[names(rlm_result) == 'p.value'] <-  'P_VAL'
  
  # Save
  new_filename <- basename(gsub('.rds', '.txt',results_files[i]))
  write.table(rlm_result, paste0(output_directory,new_filename), 
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  readr::write_rds(rlm_result, file = paste0(output_directory,basename(results_files[i])))
  
}


# -----------------------------------------------------------------------------#
#                                                                              #
#                                  ALSPAC                                      #                                         
#                                                                              #        
#------------------------------------------------------------------------------#

# Define input and output directory 
input_directory <- "results/EWAS_RESULTS/ALSPAC"
output_directory <- "results/preQC/ALSPAC/" #Note: preQC have to be created before, dir.create just creates one sub-directory not 2 

# Check if directory exists, if not create it
if (!dir.exists(output_directory)) { dir.create(output_directory) } else { print("Directory already exists")}

# Get files 
results_files <- list.files(path=input_directory, pattern = "rlm.result.model.[1-4].CT",  full.name=TRUE)
# Add file with main model and only white European ethnicity -> no need, all are white EUR
results_files <- "results/EWAS_RESULTS/ALSPAC/3_rlm.result.model.3.CT.B.rds"

for (i in seq_along(results_files)) {
  
  cat('Reading file... ', results_files[i], "\n")
  
  #rlm_result <- read.csv(results_files[i], sep=",") 
  rlm_result <- readRDS(results_files[i])
  
  #rownames(rlm_result) <- rlm_result$X
  #rlm_result <- rlm_result[,-1]
  #colnames(rlm_result) <- c('estimate','std.error', 't.value', 'p.value', 'q.value')
  
  #Renaming of variables to match 
  rlm_result$probeID <- rownames(rlm_result)
  names(rlm_result)[names(rlm_result) == 'Estimate'] <- 'BETA'
  names(rlm_result)[names(rlm_result) == 'Std. Error'] <- 'SE'
  names(rlm_result)[names(rlm_result) == 't_score'] <-  't_score'
  names(rlm_result)[names(rlm_result) == 'p_value'] <-  'P_VAL'
  names(rlm_result)[names(rlm_result) == 'q.value'] <-  'fdr'
  
  
  # Save
  # Change name to: rlm.result.model.2.CT_date_date to 3_rlm.result.model.X.CT.B.txt 
  
  new_filename <- basename(gsub('.rds', '.txt',results_files[i]))
  
  #new_filename <- paste0('3_',strsplit(basename(results_files[i]), "_")[[1]][1], '.B.txt')
  write.table(rlm_result, paste0(output_directory,new_filename), 
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  #new_filename <- paste0('3_',strsplit(basename(results_files[i]), "_")[[1]][1], '.B.rds')
  readr::write_rds(rlm_result, file = paste0(output_directory,basename(results_files[i])))
  
  
}



# -----------------------------------------------------------------------------#
#                                                                              #
#                             Generation R (GenR)                              #                                         
#                                                                              #        
#------------------------------------------------------------------------------#


# Define input and output directory 
input_directory <- "results/EWAS_RESULTS/GenR"
output_directory <- "results/preQC/GenR/" #Note: preQC have to be created before, dir.create just creates one sub-directory not 2 


# Check if directory exists, if not create it
if (!dir.exists(output_directory)) { dir.create(output_directory) } else { print("Directory already exists")}

# Get files 
results_files <- list.files(path=input_directory, pattern = "3_rlm.result.cauc.model.[1-4].CT",  full.name=TRUE)
# Add file with main model and only white European ethnicity -> all of them are
results_files <- "results/EWAS_RESULTS/GenR/3_rlm.result.cauc.model.3.CT.rds"
for (i in seq_along(results_files)) {
  
  cat('Reading file... ', results_files[i], "\n")
  
  rlm_result <- readRDS(results_files[i]) 
  colnames(rlm_result) <- c('estimate','std.error', 't.value', 'p.value', 'q.value')
  
  #Renaming of variables to match 
  rlm_result$probeID <- rownames(rlm_result)
  
  names(rlm_result)[names(rlm_result) == 'estimate'] <- 'BETA'
  names(rlm_result)[names(rlm_result) == 'std.error'] <- 'SE'
  names(rlm_result)[names(rlm_result) == 'p.value'] <-  'P_VAL'
  
  # Save
  new_filename <- basename(gsub('.rds', '.txt',results_files[i]))
  write.table(rlm_result, paste0(output_directory,new_filename), 
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  readr::write_rds(rlm_result, file = paste0(output_directory,basename(results_files[i])))
}

