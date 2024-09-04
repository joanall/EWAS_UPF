# Title: Analyze the results after the meta analysis
# Author:  JLLP 
# Date: 13/01/2023
# Description: Check the quality of the meta analysis by visualizing QQ plots
# and computing lambdas. Select the top CpGs and get Manhattan and forest plots.

## GWAMA Output Format----------------------------------------------------------

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
## ANNOTATION TO LOOK AT:  "pos" y "chr"	  
#-------------------------------------------------------------------------------



#==============================================================================#
#                           DEFINE USER PARAMETERS                             #
#==============================================================================#

cohort_files <- c(
  'new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
  'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
  'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
  'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Cohort names corresponding to the files
cohort_names <- c("HELIX", "GXXI", "ALSPAC", "GenR")

# Load the meta-analysis results from the .rds files
fixed_file <- "new_version_14_02_2024/results/GWAMA_Results/Meta_Model_3_Filtr/Meta_Model_3_Fixed_Modif.out"
random_file <- "new_version_14_02_2024/results/GWAMA_Results/Meta_Model_3_Filtr/Meta_Model_3_Random_Modif.out"


output_folder <- 'new_version_14_02_2024/results/ANALYZE_DMPs_28082024/MODEL_3/'

if (!dir.exists(output_folder)) {dir.create(output_folder, recursive = TRUE)
} else { print("Directory already exists")}


#==============================================================================#
#                                                                              #
#                           SET WORKING ENVIRONMENT                            #----
#                                                                              #
#==============================================================================#
# rm(list=ls())
# gc()

# Load libraries 
library(readr)
library(qqman)
library(dplyr)
library(tidyverse)


#------------------------------------------------------------------------------#
#                          Define Functions                                    #----
#------------------------------------------------------------------------------# 

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


make_outfile_readable <- function(filename, outfolder, return_df = FALSE) {
  
  gwama_result <- read.table(filename, header = TRUE)
  # Tranfrom to rds format
  save_to <- paste0(outfolder, gsub(".out", ".rds",basename(filename)))
  write_rds(gwama_result, file = save_to)
  
  # Transform to table 
  save_to <- paste0(outfolder, gsub(".out", ".csv",basename(filename)))
  data.table::fwrite(gwama_result, file = save_to, sep=";")
  
  if (return_df == TRUE) {return(gwama_result)}
  
}

#==============================================================================#
#                                                                              #
#                              ANALYZE DMPs                                    #----
#                                                                              #
#==============================================================================#


# Fixed effects 
fixed_results <- make_outfile_readable(fixed_file, output_folder, return_df = TRUE)  # Get data as dataframe
random_results <- make_outfile_readable(random_file, output_folder, return_df = TRUE)  # Get data as dataframe

# Get inflation factor 
lambda_ci <- as.data.frame(t(inflation(fixed_results$p.value)))
colnames(lambda_ci) <- c('lambda','upper_ci', 'lower_ci')

new_filename <- paste0(output_folder, 'lambda_', basename(gsub('.out', '.csv', fixed_effect, '.csv')))
write_csv(lambda_ci, new_filename)

# Get inflation factor 
lambda_ci <- as.data.frame(t(inflation(random_results$p.value)))
colnames(lambda_ci) <- c('lambda','upper_ci', 'lower_ci')

new_filename <- paste0(output_folder, 'lambda_', basename(gsub('.out', '.csv', random_effect, '.csv')))
write_csv(lambda_ci, new_filename)



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


#==============================================================================#
#                                                                              #
#                            GET MANHATTAN PLOT                                #----
#                                                                              #
#==============================================================================#


#------------------------------------------------------------------------------#
#                             GET MANHATTAN PLOT                               #
#------------------------------------------------------------------------------#

#results <- readRDS(file_results)
  
# Get numeric for chromosome (e.g. from "chr6" to 6) 

fixed_results$chromosome <- as.numeric(gsub('chr', "", fixed_results$chr))

# Add a column for -log10 p-values
fixed_results$logp <- -log10(fixed_results$p.value)

# Prepare data to make plot
don <- fixed_results %>%
  # Compute chromosome size
  group_by(chromosome) %>%
  summarise(chr_len = max(as.numeric(pos))) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(fixed_results, by = c("chromosome" = "chromosome")) %>%
  # Add a cumulative position of each SNP
  arrange(chromosome, pos) %>%
  mutate(BPcum = as.numeric(pos) + tot)

axisdf <- don %>%
  group_by(chromosome) %>%
  summarize(center = (max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE)) / 2)

# Define thresholds
bonferroni_threshold <- -log10(0.05 / nrow(fixed_results))
suggestive_threshold <- -log10(1E-05)

# Create the Manhattan plot
manhattan_plot <- ggplot(don, aes(x = BPcum, y = logp)) +
  geom_point(aes(color = as.factor(chromosome)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("grey", "black"), 22)) +
  scale_x_continuous(label = axisdf$chromosome, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(don$logp) + 2)) +
  geom_hline(yintercept = bonferroni_threshold, color = "#e7298a", linetype = "dashed") +
  geom_hline(yintercept = suggestive_threshold, color = "#78c679", linetype = "dashed") +
  geom_text(data = subset(don, logp > suggestive_threshold), aes(label = CpGId), vjust = -1, size = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(x = "Chromosome", y = "-log10(p)")

# Save the plot
save_to <- paste0(output_folder, gsub(".out", "_manhattan.jpeg", basename(fixed_file)))
ggsave(save_to, plot = manhattan_plot, width = 14, height = 7, dpi = 300)


#==============================================================================#
#                                                                              #
#                             GET FOREST PLOT                                  #----
#                                                                              #
#==============================================================================#


# Define the CpG site of interest
cpg_of_interest <- "cg03999434"
plot_title <- "g03999434"

#------------------------------------------------------------------------------#
#                        GET DATA INDIVIDUAL COHORTS                           #
#------------------------------------------------------------------------------#

# Initialize a list to store the cohort data
cohort_data_list <- list()

# Iterate over each file to load and filter the data
for (i in seq_along(cohort_files)) {
  # Read the cohort file
  cohort_data <- read.table(cohort_files[i], header = TRUE, stringsAsFactors = FALSE)
  
  # Filter the data for the specific CpG of interest
  filtered_data <- cohort_data[cohort_data$probeID == cpg_of_interest,] 
  
  if (nrow(filtered_data) == 0) {
    # If the CpG is not found in the cohort, create an empty row
    empty_data <- data.frame(
      Cohort = cohort_names[i],
      Effect_Size = NA,
      Standard_Error = NA
    )
    cohort_data_list[[i]] <- empty_data
  } else {
    # Add the cohort name to the dataframe
    filtered_data$Cohort <- cohort_names[i]
    
    # Select and rename relevant columns
    filtered_data$Effect_Size <- filtered_data$BETA
    filtered_data$Standard_Error <- filtered_data$SE
    
    # Subset the necessary columns
    cohort_data_selected <- filtered_data[, c("Cohort", "Effect_Size", "Standard_Error")]
    
    # Append the filtered and renamed data to the list
    cohort_data_list[[i]] <- cohort_data_selected
  }
}

# Combine all the cohort data into a single dataframe
combined_cohort_data <- do.call(rbind, cohort_data_list)

# Add lower and upper CIs to the cohort data, handling NAs appropriately
combined_cohort_data <- combined_cohort_data %>%
  mutate(
    lower_ci = ifelse(is.na(Effect_Size) | is.na(Standard_Error), NA, Effect_Size - 1.96 * Standard_Error),
    upper_ci = ifelse(is.na(Effect_Size) | is.na(Standard_Error), NA, Effect_Size + 1.96 * Standard_Error)
  )


#------------------------------------------------------------------------------#
#                           ARRANGE DATA FOR FOREST                            #
#------------------------------------------------------------------------------#


forest_data <- combined_cohort_data

# Assuming forest_data is already defined and contains your data
mtg <- meta::metagen(TE = forest_data$Effect_Size, 
                     seTE = forest_data$Standard_Error, 
                     sm = "MD", 
                     studlab = forest_data$Cohort, 
                     comb.random = TRUE, 
                     comb.fixed = TRUE)

# Manually format the effect size and CI with a semicolon and added space
effect_size_formatted <- format(round(mtg$TE, 4), nsmall = 4)
ci_lower_formatted <- format(round(mtg$lower, 3), nsmall = 3)
ci_upper_formatted <- format(round(mtg$upper, 3), nsmall = 3)

# Replace any string containing "NA" (including those with spaces) with empty strings
effect_size_formatted[grepl("NA", effect_size_formatted)] <- ""
ci_lower_formatted[grepl("NA", ci_lower_formatted)] <- ""
ci_upper_formatted[grepl("NA", ci_upper_formatted)] <- ""

# Format CI as [lower; upper] and add space between MD and CI
ci_formatted <- ifelse(ci_lower_formatted == "" & ci_upper_formatted == "",
                       "", 
                       paste0("[", ci_lower_formatted, "; ", ci_upper_formatted, "]"))
mtg$right_text <- paste(effect_size_formatted, "    ", ci_formatted)

#------------------------------------------------------------------------------#
#                                PLOTTING                                      #
#------------------------------------------------------------------------------#


save_to <- paste0(output_folder, 'forest_', cpg_of_interest, '.jpeg')

# Create the forest plot and save it as a JPEG with reduced white border and bigger plot size
jpeg(filename = save_to, width = 7, height = 5, units = "in", res = 300)
par(mar = c(3, 3, 2, 1))  # Adjust margins

# Generate the forest plot
meta::forest(mtg, 
             xlab = "Effect Size (MD)",  # Label for x-axis
             smlab = plot_title,         # Label for summary measure
             leftcols = c("studlab"), 
             leftlabs = c("Cohort"), 
             rightcols = c("right_text"), # Use custom column here
             rightlabs = c("MD        95%-CI"), # Set custom header
             fontsize = 7, 
             digits = 3,  # General digits setting (e.g., for heterogeneity stats)
             print.pval = FALSE, 
             addrow.overall = FALSE,
             col.fixed = "red", 
             col.random = "blue", 
             print.tau2 = FALSE, 
             col.diamond.fixed = "red", 
             col.diamond.random = "blue", 
             overall = TRUE, 
             test.overall = FALSE,
             print.I2 = TRUE, 
             print.pval.Q = TRUE, 
             fs.test.overall = 7, 
             fs.hetstat = 5, 
             fs.axis = 5, 
             pooled.totals = TRUE,
             text.random = "Random effects model", 
             text.fixed = "Fixed effect model")

dev.off()  # Close the JPEG device

