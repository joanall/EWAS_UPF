# Title:  Descriptive & Association Analysis
# Author:  JLLP
# Date: 29/11/2022

# Description: Script used in all cohorts to perform descriptives analysis 
# in plot and tables format. Also look at the assocoation with CEll Type. 
# Previous to run he EWAS. 


#==============================================================================#
#                           DEFINE USER PARAMETERS                             #
#==============================================================================#


# Write directory to store the results

output_directory <- "new_version_14_02_2024/results/Descriptive_Association_HELIX_28082024/"
if (!dir.exists(output_directory)) {dir.create(output_directory, recursive = TRUE)
} else { print("Directory already exists")}


# Define path of input data 
metadata_directory <- "db/final_metadata.RDS"

# If there's more than one cohort set cohort to TRUE and write the names 
cohort <- TRUE
names_cohort <- c("BIB", "EDEN", "INMA","KANC","MOBA","RHEA")


# Define variables to include in descriptive statistics table

vardesc <- c("propUPF_scaled", # exposure 
             
             "child_sex", "child_age", # child related
             "child_bmi", "child_caucasian", 
             "child_sedentary_behaviour",
             "total_vegetables_intake", "total_fruit_intake",
             
             "maternal_education", "maternal_bmi",  # mother related
             "maternal_age",  "maternal_smoking")



if (isTRUE(cohort)) {vardesc <- c(vardesc,'cohort')}

#  Define categorical variables to include in descriptive statistics plots
cat_variables <- c( 'child_sex', 'child_caucasian', 
                    'maternal_education', 'maternal_smoking',
                    "total_vegetables_intake", "total_fruit_intake")

if (isTRUE(cohort)) {cat_variables <- c(cat_variables,'cohort')}


# Define continuous variables to include in descriptive statistics plots
cont_variables <- c('propUPF_scaled', 'child_age', 
                    'child_sedentary_behaviour', 'child_bmi',
                    'maternal_bmi','maternal_age')


# Define category of ethnicity (if it applies)
# Eg. In HELIX data we have more than one ethnic group, we do association
# analysis with including only the children in the largest ethnicity

var_name_ethnicity <- 'child_caucasian'
category_ethnicity <- 'yes'






#==============================================================================#
#                                                                              #
#                           SET WORKING ENVIRONMENT                            #
#                                                                              #
#==============================================================================#


# Load libraries

library(tableone)
library(kableExtra)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(gridtext)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(plyr)

#------------------------------------------------------------------------------#
#                        Define functions                                      #
#------------------------------------------------------------------------------#



# Style title of plot 
get_title <- function(var) {
  
  if (grepl("_",var) == TRUE){
    new_name <- str_to_title(str_replace_all(var, "_", " "))
    
  } else if (grepl("_",var) == FALSE) {
    new_name <- str_to_title(var)
  }
  
  return(new_name)
}


# Do a barplot 
barplot <- function (var, data, strata = NULL) {
  
  
  if (is.null(strata) ) {
    p <- ggplot(data) +  geom_bar(aes(x=!!sym(var)), stat="count")  
    
  } else {
    
    p <- ggplot(data) +  geom_bar(aes(fill=!!sym(strata), x=!!sym(var)), stat="count",
                                  position = position_dodge(preserve = 'single'))
    
  }
  # Define plot styling 
  title_name <- get_title(var)
  p <- p + theme_tufte() +
    ggtitle(title_name) + theme(plot.title = element_text(hjust = 0.5, size = 10)) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  return (p)
  
}

# Do histogram
hist <- function (var, data, strata = NULL) {
  
  
  if (is.null(strata) ) {
    
    p <- ggplot(data) +  geom_histogram(aes(x=!!sym(var)),  bins=7)  
    
  } else {
    
    p <- ggplot(data) +  
      geom_histogram(aes(fill=!!sym(strata), color = !!sym(strata), x=!!sym(var)),  bins=7,
                                        position="identity", alpha = 0.3) 
  }
  
  # Define plot styling 
  title_name <- get_title(var)
  
  p <- p + theme_tufte() +
    ggtitle(title_name) + theme(plot.title = element_text(hjust = 0.5, size = 7)) +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) 
  
  return (p)
  
}

# Do lm with continous variables 
test_association <- function (y, exposure,vars, data) {
  
  model <- formula(paste(y, "~", paste(c(exposure,vars), collapse = "+" )))
  fit <- lm(model, data = data[ ,c(y,exposure,vars)])
  coefficients <- cbind(y,as.data.frame(summary(fit)$coefficients)[c(exposure),])
  return(coefficients) 
  
}

# Do lm with categorical variables
test_association_cat <- function (y, exposure,vars, data) {
  
  model <- formula(paste(y, "~", paste(c(exposure,vars), collapse = "+" )))
  fit <- glm(model, data = data[ ,c(y,exposure,vars)], family ='binomial')
  coefficients <- cbind(y,as.data.frame(summary(fit)$coefficients)[c(exposure),])
  return(coefficients) 
  
}

#==============================================================================#
#                                                                              #
#                           DESCRIPTIVE ANALYSIS                               #----
#                                                                              #
#==============================================================================#

# Load metadata (exposure + methylation variables)
metadata <- readRDS(file = metadata_directory) 

#------------------------------------------------------------------------------#
#                           Descriptive Tables                                 #----
#------------------------------------------------------------------------------#
descriptive_object <- CreateTableOne(vars = vardesc, 
                                     data = metadata, includeNA=T) 

descriptive_table <- print(descriptive_object,  noSpaces = TRUE, varLabels=T)

# Save descriptive table 
write.table(descriptive_table , 
            file = paste0(output_directory,"/Descriptive_Table_non_stratified.csv"), 
            sep=";")


#------------------------------------------------------------------------------#
#                           Descriptive Plots                                  #----
#------------------------------------------------------------------------------#

# Get barplots  for categorical variables 
barplots <- lapply(cat_variables, function(var) barplot(var, metadata)) 

# Get histograms  for categorical variables 
histplots <- lapply(cont_variables, function(var) hist(var, metadata))

# Join plots and save them 
descriptive_plots <- c(barplots, histplots)

grid_descriptive_plots <- gridExtra::grid.arrange(grobs = descriptive_plots,
                                                  ncol = 3)

ggsave(grid_descriptive_plots,height = 8,
       file = paste0(output_directory,"/Descriptive_Plots_non_stratified.pdf"))






#==============================================================================#
#                                                                              #
#                       DESCRIPTIVE ANALYSIS by COHORT                         #-----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                           Descriptive Tables                                 #----
#------------------------------------------------------------------------------#

if (isTRUE(cohort)) {
  
  # Get summary statistics stratified By cohort
  descriptive_object <- CreateTableOne(vardesc[-2], strata="cohort", 
                                       data=metadata, includeNA =T)
  
  descriptive_table <- print(descriptive_object,  noSpaces = TRUE, varLabels=T)
  
  write.table(descriptive_table, file = paste0(output_directory,"/Descriptive_Table_by_cohort.csv"),
              sep=";")
  
}  



#------------------------------------------------------------------------------#
#                           Descriptive Plots                                  #----
#------------------------------------------------------------------------------#


if (isTRUE(cohort)) {
  
  # Get barplots  for categorical variables stratified by cohort 
  
  barplots <- lapply(cat_variables[-1], #not include h_cohort var
                     function(var) barplot(var, metadata, strata= 'cohort')) 
  
  # Get histograms  for categorical variables stratified by cohort 
  
  histplots <- lapply(cont_variables, function(var) hist(var, metadata, 'cohort')) 
  
  # Join plots and save them 
  
  descriptive_plots <- c(barplots, histplots)
  
  
  # Get common legend (1/2)
  
  leg <- get_legend(descriptive_plots[[1]]) # extract legend from one plot
  
  legend <- as_ggplot(leg)  # turn into gg object
  
  
  # Remove legend from each plot
  
  descriptive_plots_no_legend <- lapply(descriptive_plots,function(x) x + 
                                          theme( legend.position = "none")) 
  
  
  # Join plots 
  grid_descriptive_plots <- gridExtra::grid.arrange(grobs = descriptive_plots_no_legend,
                                                    legend,
                                                    ncol = 3)
  
  
  ggsave(grid_descriptive_plots,height = 8,
         file = paste0(output_directory,"/Descriptive_Plots_by_cohort.pdf"))
  
}

#===============================================================================

#==============================================================================#
#                                                                              #
#                          ASSOCIATION ANALYSIS                                #-----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                           Define models                                      #----
#------------------------------------------------------------------------------#

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

# Declare cell type vars 
veg_fruit_intake <- c("total_fruit_intake", "total_vegetables_intake")


# Declare cell type vars 
cell_types <- c("NK_6", "Bcell_6", "CD4T_6", "CD8T_6", "Gran_6", "Mono_6")



#------------------------------------------------------------------------------#
#                     Test association with Cell Types                         #----
#------------------------------------------------------------------------------#

# All data 
for (i in 1:length(list_models)) {
  cellfits <- lapply (cell_types, 
                      function(y) test_association(y, exposure, 
                                                   list_models[[i]], 
                                                   metadata))
  results <- ldply (cellfits, data.frame) # transform list to dataframe
  write.csv(results, 
            file = paste0(output_directory, "Association_CT_vars_model.", i, ".csv"),  quote=FALSE) 

}


#------------------------------------------------------------------------------#
#                     Test association with Vegetables & Fruit                 #----
#------------------------------------------------------------------------------#

for (i in 1:length(list_models)) {
  health_fits <- lapply(veg_fruit_intake, function(y) {
    test_association_cat(y, exposure, list_models[[i]], metadata)
  })
  
  results <- ldply(health_fits, data.frame) # transform list to dataframe
  
  # Save the results without the sep argument
  write.csv(results, 
            file = paste0(output_directory, "Association_veg_fruit_upf_model.", i, ".csv"),  
            quote = FALSE)
}



#==============================================================================#
#                                                                              #
#                     ASSOCIATION ANALYSIS by ETHNICTY                         #-----
#                                                                              #
#==============================================================================#


if (exists("var_name_ethnicity") && var_name_ethnicity != "") {
  # Exclude non-caucasian children 
  metadata_copy <- data.frame(metadata)
  
  metadata <- metadata_copy %>%
    dplyr::filter(!!sym(var_name_ethnicity) == category_ethnicity)  # keep only caucasian children
  
  output_directory <- paste0(output_directory, "Ethn_")
  
  # Test association for cell types
  for (i in 3:length(list_models)) {
    cellfits <- lapply(cell_types, function(y) 
      test_association(y, exposure, list_models[[i]], metadata))
    
    results <- ldply(cellfits, data.frame) # transform list to dataframe
    write.csv(results, 
              file = paste0(output_directory, "Association_CT_vars_model.", i, ".csv"),  
              quote = FALSE) 
  }
  
  # Test Association with vegetables and fruit 
  for (i in 1:length(list_models)) {
    health_fits <- lapply(veg_fruit_intake, function(y) 
      test_association_cat(y, exposure, list_models[[i]], metadata))
    
    results <- ldply(health_fits, data.frame) # transform list to dataframe
    write.csv(results, 
              file = paste0(output_directory, "Association_veg_fruit_upf_model.", i, ".csv"),  
              quote = FALSE) 
  }
}



