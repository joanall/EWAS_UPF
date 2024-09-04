# Title: Quality Control of EWAS result 
# Author:  JLLP
# Date: 03/01/2023

# Description: Read results of EWAS and :
# 1) Compute the adjusted P values if it's not done, code example as follows
# df <- df %>%
#   dplyr::mutate(q.value = p.adjust(p.value, method = "BH")) %>%
#   dplyr::arrange(p.value) #order ascending by p value

# 2) Make it readable format to run the Quality Control of EASIER package. 

# 3) Run Quality Control based on EASIER package

#==============================================================================#
#                           SET WORKING ENVIRONMENT                            #
#==============================================================================#

library(readr)
library(EASIER)

##  Uncomment this code to install EASIER package
#
# # Install devtools
# install.packages("devtools")
#
# # Install required packages
# devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")

# # Install EASIER package
# devtools::install_github("isglobal-brge/EASIER@HEAD")


#==============================================================================#
#                                                                              #
#                           PREPARE FILES FOR QC                               #----
#                                                                              #
#==============================================================================#


#------------------------------------------------------------------------------#
#                                   HELIX                                      #----                                      
#------------------------------------------------------------------------------#

# Define input and output directory 
input_directory <- "new_version_14_02_2024/results/EWAS_RESULTS/HELIX/"
output_directory <- "new_version_14_02_2024/results/QC_28082024/HELIX/"  

# Check if directory exists, if not create it
if (!dir.exists(output_directory)) { 
  dir.create(output_directory, recursive = TRUE) 
  } else {print("Directory already exists")}


# Define EWAS results files
results_files <- c(paste0(input_directory, '3_rlm.result.model.1.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.2.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.3.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.4.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.cauc.model.3.CT.B.rds'))

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
  
  readr::write_rds(rlm_result, file = paste0(output_directory, 
                                             basename(results_files[i])))
  
}


#------------------------------------------------------------------------------#
#                             Generation XXI (GXXI)                            #----                                         
#------------------------------------------------------------------------------#

# Define input and output directory 
input_directory <- "new_version_14_02_2024/results/EWAS_RESULTS/GXXI/"
output_directory <- "new_version_14_02_2024/results/QC_28082024/GXXI/"  


# Check if directory exists, if not create it
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE) 
  } else {print("Directory already exists")}


# Define EWAS results files
results_files <- c(paste0(input_directory, '3_rlm.result.model.1.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.2.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.3.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.4.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.cauc.model.3.CT.B.rds'))

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
  
  readr::write_rds(rlm_result, file = paste0(output_directory,
                                             basename(results_files[i])))
  
}


#------------------------------------------------------------------------------#
#                                  ALSPAC                                      #----                                         
#------------------------------------------------------------------------------#

# Define input and output directory 
input_directory <- "new_version_14_02_2024/results/EWAS_RESULTS/ALSPAC/"
output_directory <- "new_version_14_02_2024/results/QC_28082024/ALSPAC/"  

# Check if directory exists, if not create it
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE) 
  } else { print("Directory already exists")}


# Define EWAS results files
results_files <- c(paste0(input_directory, '3_rlm.result.model.1.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.2.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.3.CT.B.rds'),
                   paste0(input_directory, '3_rlm.result.model.4.CT.B.rds'))


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
  new_filename <- basename(gsub('.rds', '.txt',results_files[i]))
  
  write.table(rlm_result, paste0(output_directory,new_filename), 
              append = FALSE, sep = "\t", dec = ".",
              row.names = TRUE, col.names = TRUE)
  
  readr::write_rds(rlm_result, file = paste0(output_directory,basename(results_files[i])))
  
  
}



#------------------------------------------------------------------------------#
#                             Generation R (GenR)                              #----                                         
#------------------------------------------------------------------------------#


# Define input and output directory 
input_directory <- "new_version_14_02_2024/results/EWAS_RESULTS/GenR/"
output_directory <- "new_version_14_02_2024/results/QC_28082024/GenR/"  

# Check if directory exists, if not create it
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE) 
} else { print("Directory already exists")}


# Define EWAS results files
results_files <- c(paste0(input_directory, '3_rlm.result.cauc.model.1.CT.rds'),
                   paste0(input_directory, '3_rlm.result.cauc.model.2.CT.rds'),
                   paste0(input_directory, '3_rlm.result.cauc.model.3.CT.rds'),
                   paste0(input_directory, '3_rlm.result.cauc.model.4.CT.rds'),
                   paste0(input_directory, '3_rlm.result.cauc.model.3.CT.rds'))

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

#===============================================================================

#==============================================================================#
#                                                                              #
#                              QUALITY CONTROL                                 #----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                         DEFINE FUNCTION QC                                   #----                                      
#------------------------------------------------------------------------------#

# The function is copied from the EASIER package developed by Dolors Pelegrí

do_quality_control <- function(files, results_folder, prefixes, venn_diagrams,
                               artype, N, n, colname_NforProbe, pcMissingSamples,
                               exclude) {
  
  # Variable declaration to perform precision plot
  medianSE <- numeric(length(files))
  value_N <- numeric(length(files))
  
  if(length(n) == length(N))
    value_n <- numeric(length(files))
  
  cohort_label <- character(length(files))
  
  # Prepare output folder for results (create if not exists)
  if(!dir.exists(file.path(getwd(), results_folder )))
    suppressWarnings(dir.create(file.path(getwd(), results_folder)))
  
  ## Remove duplicates, Exclude CpGs and adjust data (BN and FDR)
  for ( i in 1:length(files) )
  {
    
    # Prepare output subfolder for cohort-model results (create if not exists)
    if(!dir.exists(file.path(getwd(), results_folder, prefixes[i] )))
      suppressWarnings(dir.create(file.path(getwd(), results_folder, prefixes[i])))
    
    # Creates an empty file to resume all data if an old file exist  removes
    # the file and creates a new one
    fResumeName <- paste0( file.path(getwd(), results_folder, prefixes[i]),"/",prefixes[i], "_descriptives.txt")
    if ( file.exists(fResumeName) ) {
      file.remove(fResumeName)
    }
    file.create(fResumeName)
    
    # Read data.
    cohort <- read.table(files[i], header = TRUE, as.is = TRUE)
    print(paste0("Cohort file : ",files[i]," - readed OK", sep = " "))
    
    # Remove rows with NA from data
    cohort <- clean_NA_from_data(cohort)
    
    # Descriptives - Before CpGs deletion #
    cohort <- descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], artype = artype[i], before = TRUE)
    
    # Remove duplicates
    # cohort <- remove_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i], '/',prefixes[i],'_descriptives_duplic.txt'), paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
    
    test_duplicate_CpGs(cohort, "probeID", paste0(results_folder,'/',prefixes[i],'_duplicates.txt') )
    
    # Remove cpGs with low representation
    # first, we test if colname_NforProbe and pcMissingSampes are defined
    if( !exists("colname_NforProbe") ) { colname_NforProbe <- NULL }
    if( !exists("pcMissingSamples") ) { pcMissingSamples <- NULL }
    
    cohort <- filterLowRepresentedCpGsinCohort(cohort, colname_NforProbe, pcMissingSamples, N[i], fileresume = fResumeName )
    
    # Exclude CpGs not meet conditions
    cat ('Excluding CpGs... ')
    if("MASK_snp5_ethnic" %in% exclude ){
      cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = ethnic[i], filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
    } else {
      #..# if( !is.null(exclude) && exclude!='') {
      cohort <- exclude_CpGs(cohort, "probeID", exclude, ethnic = "", filename = paste0(results_folder, '/',prefixes[i], '/',prefixes[i],'_excluded.txt'), fileresume = fResumeName, artype = artype[i] )
      #..# }
    }
    
    
    # Descriptives - After CpGs deletion #
    cat ('Writing descriptives... ')
    descriptives_CpGs(cohort, c("BETA", "SE", "P_VAL"), fResumeName, N[i], artype = artype[i], before = FALSE )
    
    # Adjust data by Bonferroni and FDR
    cat ('Adjusting P values ... ')
    cohort <- adjust_data(cohort, "P_VAL", bn=TRUE, fdr=TRUE, fResumeName, N[i]  )
    
    # Write QC complete data to external file
    cat ('Saving data... ')
    write_QCData(cohort, paste0(results_folder, '/',prefixes[i], '/',prefixes[i]))
    
    ## Visualization - Plots
    cat ('Making visualization plots... ')
    #. Problems in some workstations and servers.# rasterpdf::raster_pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'), res = 300)
    #..# Problems in some cases --> Get png plots : pdf(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QCplots.pdf'))
    
    # Distribution plot
    png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_SE_plot.png'), type="cairo")
    plot_distribution(cohort$SE, main = paste('Standard Errors of', prefixes[i]), xlab = 'SE')
    dev.off()
    png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_pvals_plot.png'), type="cairo")
    plot_distribution(cohort$P_VAL, main = paste('p-values of', prefixes[i]), xlab = 'p-value')
    dev.off()
    
    # QQ plot
    png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_QQ_plot.png'), type="cairo")
    qqman::qq(cohort$P_VAL, main = sprintf('QQ plot of %s (lambda = %f)', prefixes[i], lambda=get_lambda(cohort,"P_VAL")))
    dev.off()
    
    # Volcano plot.
    png(paste0(results_folder, '/',prefixes[i], '/',prefixes[i], '_QC_Volcano_plot.png'), type="cairo")
    plot_volcano(cohort, "BETA", "P_VAL", main =paste('Volcano plot of', prefixes[i]) )
    dev.off()
    
    # Add mandatory data for precisionplot
    medianSE[i] <-  median(cohort$SE)
    value_N[i] <- N[i]
    cohort_label[i] <-  prefixes[i]
    
    # if n is defined for dichotomic condition :
    if(length(n) == length(N))  value_n[i] <- n[i]
    
    # Store data for Beta Box-Plot
    if( i == 1)
      betas.data <- list()
    betas.data[[prefixes[i]]] <- cohort[,"BETA"]
    
  }
  
  # Create QC Summary -> no se'm crea pro..
  cat ('Creating QC summary ... ')
  if( file.exists(paste0(results_folder,"/tmp_pretQC.txt")) & file.exists(paste0(results_folder,"/tmp_postQC.txt")) & file.exists(paste0(results_folder,"/tmp_postQCAdj.txt") )) {
    preQC <- read.table (file = paste0(results_folder,"/tmp_pretQC.txt"), header = TRUE, sep = "\t")
    postQC <- read.table (file = paste0(results_folder,"/tmp_postQC.txt"), header = TRUE, sep = "\t")
    postQCAdj <- read.csv(file = paste0(results_folder,"/tmp_postQCAdj.txt"), header = TRUE, sep = "\t")
    
    write.table( cbind(prefixes,ethnic, postQC[,1:2], preQC, postQC[,3:length(postQC)], postQCAdj), file = paste0(results_folder,"/Summary_QCs.txt" ),
                 row.names = FALSE, col.names = TRUE, sep = "\t")
    
    do.call(file.remove, list(list.files(results_folder, full.names = TRUE, pattern = "tmp_*")))
  }
  
  
  # Data for Precision Plot
  cat ('Creating precision plot ... ')
  precplot.data <- cbind.data.frame( SE = medianSE, invSE = (1/medianSE), N = value_N, sqrt_N = sqrt(N), cohort = cohort_label )
  cols.numeric <- c("SE","invSE", "N", "sqrt_N")
  precplot.data[cols.numeric] <- sapply(precplot.data[cols.numeric],as.numeric)
  
  if(length(n) == length(N)){
    precplot.data.n <- cbind.data.frame( SE = medianSE, invSE = (1/medianSE), N = value_n, sqrt_N = sqrt(n), cohort = cohort_label )
    precplot.data.n[cols.numeric] <- sapply(precplot.data.n[cols.numeric],as.numeric)
  }
  
  # BoxPlot with Betas in all Models and cohorts
  cat ('Creating boxplot ... ')
  plot_betas_boxplot(betas.data, paste(results_folder, 'BETAS_BoxPlot.png', sep="/"))
  
  ##  Post model analysis  ## -> aqui peta algo 
  if ( length(files) > 1)
  {
    # Precision_Plot(N)
    plot_precisionp(precplot.data, paste(results_folder,  "precision_SE_N.png",sep='/'), main = "Precision Plot - 1/median(SE) vs sqrt(n)")
    
    # Precision_Plot(n)
    if(length(n) == length(N))
      plot_precisionp(precplot.data.n, paste(results_folder,  "precision_SE_N.png", sep='/'), main = "Subgroup Precision Plot -  1/median(SE) vs sqrt(n)")
    
    # Venn_Diagrams()
    for (i in 1:length(venn_diagrams))
      plot_venndiagram(venn_diagrams[[i]], qcpath = results_folder, plotpath = results_folder, bn='padj.bonf', fdr='padj.fdr')
    
  }
}


#------------------------------------------------------------------------------#
#                         Define common parameters                             #----                                      
#------------------------------------------------------------------------------#


# Minimum sample representation percentage required for CpGs
colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9 # it is not used here but in the next step that is the meta-analysis 

# Exclude  
exclude <- c('control_probes','noncpg_probes','MASK_mapping','MASK_sub30_copy',
             'MASK_extBase','MASK_typeINextBaseSwitch','Unrel_450_EPIC_blood','Sex')



#------------------------------------------------------------------------------#
#                                 MODEL 1                                      #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.1.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.1.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.1.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.1.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/MODEL_1' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI','ALSPAC' , 'GenR')
venn_diagrams <- list(c('HELIX', 'GXXI','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K','450K')
N <- c(1138,732,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)



#------------------------------------------------------------------------------#
#                                 MODEL 2                                      #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.2.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.2.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.2.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.2.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/MODEL_2' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI','ALSPAC' , 'GenR')
venn_diagrams <- list(c('HELIX', 'GXXI','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K','450K')
N <- c(1138,732,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)



#------------------------------------------------------------------------------#
#                                 MODEL 3                                      #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/MODEL_3' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI','ALSPAC' , 'GenR')
venn_diagrams <- list(c('HELIX', 'GXXI','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K','450K')
N <- c(1138,732,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)



#------------------------------------------------------------------------------#
#                   SENSITIVITY 1 (healthy diet)                               #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.4.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.4.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.4.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.4.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/HEALTHY_DIET' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI','ALSPAC' , 'GenR')
venn_diagrams <- list(c('HELIX', 'GXXI','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K','450K')
N <- c(1138,732,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)


#------------------------------------------------------------------------------#
#                       SENSITIVITY 2 (ethnicity)                              #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.cauc.model.3.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.cauc.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/ETHNICITY' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI','ALSPAC' , 'GenR')
venn_diagrams <- list(c('HELIX', 'GXXI','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K','450K')
N <- c(1014,711,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)



#------------------------------------------------------------------------------#
#                             Out HELIX                                        #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_HELIX' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('GXXI','ALSPAC' , 'GenR')
venn_diagrams <- list(c('GXXI','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('EPIC','450K','450K')
N <- c(732,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)


#------------------------------------------------------------------------------#
#                              Out GXXI                                        #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_GXXI' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX','ALSPAC' , 'GenR')
venn_diagrams <- list(c('HELIX','ALSPAC' , 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','450K','450K')
N <- c(1138,892,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)

#------------------------------------------------------------------------------#
#                             Out ALSPAC                                       #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_ALSPAC' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI', 'GenR')
venn_diagrams <- list(c('HELIX', 'GXXI', 'GenR'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K')
N <- c(1138,732,390) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)

#------------------------------------------------------------------------------#
#                             Outt Gen R                                       #----                                      
#------------------------------------------------------------------------------#

# Input files to perform QC 
files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
           'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.cauc.model.3.CT.B.txt',
           'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt')

# Result folder
results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_GenR' #if doesn't exist, it creates it

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI','ALSPAC' )
venn_diagrams <- list(c('HELIX', 'GXXI','ALSPAC'))

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K')
N <- c(1138,732,892) # sample size x cohort
n <- c(NA) 


do_quality_control(files = files, results_folder = results_folder,
                   prefixes = prefixes, venn_diagrams = venn_diagrams,
                   artype = artype, N = N, n = n, 
                   colname_NforProbe = colname_NforProbe, 
                   pcMissingSamples = pcMissingSamples,
                   exclude = exclude)


