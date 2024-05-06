# Author:  NS & JLLP 
# Title:  Association analysis between CELL types and variables  
# Date: 29/11/2022

# Description : Check association between UPF and covariates to avoid cofounding. 
# Apply RL model and adjust for variables use in future main model in ewas.


################################################################################
#                         SET WORKING ENVIRONMENT                              #
################################################################################

# Load libraries

library(tidyverse)
library(plyr)
library(kableExtra)


################################################################################
#                             DEFINE VARIABLES                                 #
################################################################################


output_directory <- "results/Descriptive_Association_HELIX/"

if (!dir.exists(output_directory)) { dir.create(output_directory)
  
} else { print("Directory already exists")}


metadata_directory <- "db/final_metadata.RDS"

# ethnicity <- "caucasian"  
ethnicity <- "all" 


################################################################################
#                                  LOAD DATA                                   #
################################################################################


# Load metadata (exposure + methylation variables)

metadata <- readRDS(file = metadata_directory)

if (ethnicity == "caucasian") {
  
  # Exclude non-caucasian children 
  
  metadata_copy <- data.frame(metadata)
  
  metadata <- metadata_copy %>%
    dplyr::filter(child_caucasian == "yes")  # keep only caucasian children
  
  output_directory <- paste0(output_directory, "2_Cauc_")
  
} else if (ethnicity != "caucasian") {
  output_directory <- paste0(output_directory, "2_All_")}




################################################################################
#                       DEFINE MODELS TO TEST ASSOCIATION                      #
################################################################################
  
exposure <- "propUPF_scaled"

model.1 <- c("cohort", "child_sex", "child_age") # basic

model.2 <- c(model.1,
             "maternal_age", "maternal_bmi", "maternal_smoking",
             "maternal_education") # maternal related vars

model.3 <- c(model.2,
             "child_bmi","child_sedentary_behaviour") # child related vars

model.4 <- c(model.3,
             "total_vegetables_intake","total_fruit_intake") # diet related vars

list_models <- list(model.1, model.2, model.3, model.4)




################################################################################
#                  Examine associations of UPF with cell types                 #
################################################################################

# Declare cell type vars 

cell_types <- c("NK_6", "Bcell_6", "CD4T_6", "CD8T_6", "Gran_6", "Mono_6")

# All data 

test_association <- function (y, exposure,vars, data) {
  
  model <- formula(paste(y, "~", paste(c(exposure,vars), collapse = "+" )))
  fit <- lm(model, data = data[ ,c(y,exposure,vars)])
  coefficients <- cbind(y,as.data.frame(summary(fit)$coefficients)[c(exposure),])
  return(coefficients) 
  
}

# Get associacation for all models 
  
for (i in 3:length(list_models)) {
  cellfits <- lapply (cell_types, 
                      function(y) test_association(y, exposure, 
                                                   list_models[[i]], 
                                                   metadata))
  results <- ldply (cellfits, data.frame) # transform list to dataframe
  write.csv(results, 
            file = paste0(output_directory, "Association_CT_vars_model.", i, ".csv"),  quote=FALSE) 
  
  kable(results,booktabs = T) %>%
    kable_styling(bootstrap_options = c("condensed"), # rows a little more tight
                  font_size = 11.5,
                  full_width = F,
                  row_label_position = "c") %>%
    add_header_above(c("Descriptive Statistics" = 5 )) %>%
    save_kable(paste0(output_directory, "Association_CT_vars_model.", i, ".pdf"))
}


################################################################################
#          Examine associations of UPF with veggetables & fruit intake         #
################################################################################

# Declare cell type vars 

veg_fruit_intake <- c("total_fruit_intake", "total_vegetables_intake")

test_association_cat <- function (y, exposure,vars, data) {
  
  model <- formula(paste(y, "~", paste(c(exposure,vars), collapse = "+" )))
  fit <- glm(model, data = data[ ,c(y,exposure,vars)], family ='binomial')
  coefficients <- cbind(y,as.data.frame(summary(fit)$coefficients)[c(exposure),])
  return(coefficients) 
  
}

i <- 1
health_fits <- lapply (veg_fruit_intake, 
                       function(y) test_association_cat(y, exposure, 
                                                        list_models[[i]], 
                                                        metadata))

results <- ldply (health_fits, data.frame) # transform list to dataframe



for (i in 1:length(list_models)) {
  health_fits <- lapply (veg_fruit_intake, 
                         function(y) test_association_cat(y, exposure, 
                                                          list_models[[i]], 
                                                          metadata))
  results <- ldply (health_fits, data.frame) # transform list to dataframe
  write.csv(results, 
            file = paste0(output_directory, "Association_veg_fruit_upf_model.", i, ".csv"),  quote=FALSE) 
  
  kable(results,booktabs = T, digits = 3) %>%
    kable_styling(bootstrap_options = c("condensed"), # rows a little more tight
                  font_size = 11.5,
                  full_width = F,
                  row_label_position = "c") %>%
    add_header_above(c("Descriptive Statistics" = 5 )) %>%
    save_kable(paste0(output_directory, "Association_veg_fruit_upf_vars_model.", i, ".pdf"))
}




