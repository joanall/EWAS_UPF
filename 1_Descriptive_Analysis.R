# Author:  JLLP
# Title:  Get descriptives of our data 
# Date: 29/11/2022

# Description: Get the descirptives (tables & plots) of data from HELIX cohort
# used to asses assocaition UPF and DNAm. 

################################################################################
#                           SET WORKING ENVIRONMENT                            #
################################################################################

# Load libraries

library(tableone)
library(kableExtra)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(gridtext)
library(ggpubr)
library(dplyr)

################################################################################
#                           DEFINE USER PARAMETERS                             #
################################################################################


#  Define output/results directory 

output_directory <- "results/Descriptive_Association_HELIX/"

# Check if directory exists, if not create it

if (!dir.exists(output_directory)) { dir.create(output_directory)
  
} else { print("Directory already exists")}


# Define path of metadata 

metadata_directory <- "db/final_metadata.RDS"

# If there's more than one cohort or not (TRUE or FALSE, by default is FALSE) 

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



if (isTRUE(cohort)) { vardesc <- c(vardesc,'cohort')}

#  Define categorical variables to include in descriptive statistics plots


cat_variables <- c( 'child_sex', 'child_caucasian', 
                    'maternal_education', 'maternal_smoking',
                    "total_vegetables_intake", "total_fruit_intake")

if (isTRUE(cohort)) { cat_variables <- c(cat_variables,'cohort')}


# Define continuous variables to include in descriptive statistics plots

cont_variables <- c('propUPF_scaled', 'child_age', 
                    'child_sedentary_behaviour', 'child_bmi',
                    'maternal_bmi','maternal_age')



################################################################################
#                           DEFINE FUNCTIONS                                   #
################################################################################


# Define function style title of plot 

get_title <- function(var) {
  
  if (grepl("_",var) == TRUE){
    new_name <- str_to_title(str_replace_all(var, "_", " "))
    
  } else if (grepl("_",var) == FALSE) {
    new_name <- str_to_title(var)
  }
  
  return(new_name)
}


# Define function to plot a barplot 

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

# Define function to plot a histogram  

hist <- function (var, data, strata = NULL) {
  
  
  if (is.null(strata) ) {
    
    p <- ggplot(data) +  geom_histogram(aes(x=!!sym(var)),  bins=7)  
    
  } else {
    
    p <- ggplot(data) +  geom_histogram(aes(fill=!!sym(strata), color = !!sym(strata), x=!!sym(var)),  bins=7,
                                        #position = position_dodge(preserve = 'single')
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


################################################################################
#                           DEFINE FUNCTIONS                                   #
################################################################################

# Load metadata (exposure + methylation variables)

metadata <- readRDS(file = metadata_directory ) # 1138 x 34



################################################################################
#               GET DESCRIPTIVE TABLE NON - STRATIFIED                         #
################################################################################


# Non stratified 

descriptive_object <- CreateTableOne(vars = vardesc, 
                                     data = metadata, includeNA=T) 

descriptive_table <- print(descriptive_object,  noSpaces = TRUE, varLabels=T)


# Save descriptive table 

write.table(descriptive_table , 
            file = paste0(output_directory,"/1_Descriptive_Table_non_stratified.csv"), sep=";")




# Formatting table - rename rownames 

# Note: Careful with the order, if it's different names will not match/be correct


if (isFALSE(cohort)){
  
  rownames(descriptive_table_formatted) <- c("N" ="n",
                                             "Scaled Proportion of UPF (mean (SD)) (UPF ingested / total servings)",
                                             "Child Sex = male (%)",
                                             "Child Age (mean (SD)) (years)", 
                                             "Child BMI (mean (SD)) (kg/m2)",
                                             "Child Caucasian = yes (%)",
                                             "Child Sedentary Behaviour (mean (SD)) (min/day)",
                                             "Maternal Education (%)",
                                             "Low", "Middle", "High",
                                             "Maternal BMI (mean (SD)) (kg/m2)",
                                             "Maternal Age (mean (SD))) (years)", 
                                             "Maternal Smoking = yes (%)")
}

if (isTRUE(cohort)){
  rownames(descriptive_table_formatted) <- c("N",
                                             "Scaled Proportion of UPF (mean (SD)) (UPF ingested / total servings)",
                                             "Child Sex = male (%)",
                                             "Child Age (mean (SD)) (years)", 
                                             "Child BMI (mean (SD)) (kg/m2)",
                                             "Child Caucasian = yes (%)",
                                             "Child Sedentary Behaviour (mean (SD)) (min/day)",
                                             "Total Vegetables Intake (%) (times / week)",
                                             "< 6.0", "6.0-8.5", "> 8.5",
                                             "Total Fruit Intake (%) (times / week)",
                                             "< 7", "7-14.1","> 14.1",
                                             "Maternal Education (%)",
                                             "Low", "Middle", "High",
                                             "Maternal BMI (mean (SD)) (kg/m2)",
                                             "Maternal Age (mean (SD))) (years)", 
                                             "Maternal Smoking = yes (%)",
                                             "Cohort (%)", 
                                             names_cohort) }  



# Formatting table - Styling 

kable(descriptive_table_formatted,booktabs = T) %>%
  kable_styling(bootstrap_options = c("condensed"), # rows a little more tight
                font_size = 11.5,
                full_width = F,
                row_label_position = "c") %>%
  add_header_above(c("Descriptive Statistics" = 2 )) %>%
  save_kable(paste0(output_directory,"/1_Descriptive_Table_non_stratified.pdf"))



################################################################################
#                  DESCRIPTIVE PLOTS NON-STRATIFIED                            #
################################################################################

# Get barplots  for categorical variables 

barplots <- lapply(cat_variables, function(var) barplot(var, metadata)) 

# Get hsitograms  for categorical variables 

histplots <- lapply(cont_variables, function(var) hist(var, metadata))

# Join plots and save them 

descriptive_plots <- c(barplots, histplots)

grid_descriptive_plots <- gridExtra::grid.arrange(grobs = descriptive_plots,
                                                  ncol = 3)

ggsave(grid_descriptive_plots,height = 8,
       file = paste0(output_directory,"/1_Descriptive_Plots_non_stratified.pdf"))








################################################################################
#                    GET DESCRIPTIVE TABLE STRATIFIED                          #
################################################################################


if (isTRUE(cohort)) {
  
  # Get summary statistics stratified By cohort
  
  descriptive_object <- CreateTableOne(vardesc[-2], strata="cohort", 
                                       data=metadata, includeNA =T)
  
  descriptive_table <- print(descriptive_object,  noSpaces = TRUE, varLabels=T)
  
  write.table(descriptive_table, file = paste0(output_directory,"/1_Descriptive_Table_by_cohort.csv"),sep=";")
  
  
  # Formatting table - rename rownames 
  
  # Note: Careful with the order, if it's different names will not match/be correct
  
  descriptive_table_formatted <- descriptive_table
  rownames(descriptive_table_formatted) <- c("N",
                                             "Scaled Proportion of UPF (mean (SD)) (UPF ingested / total servings)",
                                             "Child Sex = male (%)",
                                             "Child Age (mean (SD)) (years)", 
                                             "Child BMI (mean (SD)) (kg/m2)",
                                             "Child Caucasian = yes (%)",
                                             "Child Sedentary Behaviour (mean (SD)) (min/day)",
                                             "Total Vegetables Intake (%) (times / week)",
                                             "< 6.0", "6.0-8.5", "> 8.5",
                                             "Total Fruit Intake (%) (times / week)",
                                             "< 7", "7-14.1","> 14.1",
                                             "Maternal Education (%)",
                                             "Low", "Middle", "High",
                                             "Maternal BMI (mean (SD)) (kg/m2)",
                                             "Maternal Age (mean (SD))) (years)", 
                                             "Maternal Smoking = yes (%)")
  
  
  
  # Formatting table - Styling (not finished)
  
  kable(descriptive_table_formatted,booktabs = T) %>%
    kable_styling(bootstrap_options = c("condensed"), # rows a little more tight
                  font_size = 11.5,
                  full_width = F,
                  row_label_position = "c") %>%
    add_header_above(c("Descriptive Statistics Stratified" = 9 )) %>%
    save_kable(paste0(output_directory,"/1_Descriptive_Table_by_cohort.pdf"))
}




################################################################################
#                     DESCRIPTIVE PLOTS STRATIFIED                             #
################################################################################

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
         file = paste0(output_directory,"/1_Descriptive_Plots_by_cohort.pdf"))
  
}