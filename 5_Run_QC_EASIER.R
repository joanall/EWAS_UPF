# Author:  Dolors Pelegrí
# Title: Quality Control Script to use with EASIER package
# Date: 03/01/2023

################################################################################
#                      SET WORKING ENVIRONMENT                                 #
################################################################################

## -------------------------------------
##  Install EASIER Package Code
## -------------------------------------
##
##  Uncomment this code to install EASIER package
#
# # Install devtools
# install.packages("devtools")
#
# # Install required packages
# devtools::source_url("https://raw.githubusercontent.com/isglobal-brge/EASIER/HEAD/installer.R")

# # Install EASIER package
# devtools::install_github("isglobal-brge/EASIER@HEAD")

##  END -  Install EASIER Package Code
## -------------------------------------


library(EASIER)


# Note: From results, we manually create at the beginning QC_EASIER. 
# The outputs folder will be like results/QC_EASIER/model_name 
# This way we have for each model the results for each cohort

# -----------------------------------------------------------------------------#
#                                                                              #
#               INFORMATION TO DEFINE VARIABLES CORRECTLY                      #                                         
#                                                                              #        
#------------------------------------------------------------------------------#

# Input files PATH:
# results/preQC/HELIX/ ; results/preQC/GXXI/ ; results/preQC/ALSPAC/ ; results/preQC/GenR/ ; results/preQC/subHELIX/

# PARTICIPATING COHORTS DATA: 

# Prefixes for each file (always same order!)
# -> all files will be like prefix_descriptives.txt, prefix_excluded.txt etc. 
# -  HELIX, GXII, ALSPAC, GenR 
# N : 1138, 732,  892,   390 
# Models can be with or withou CT. EX:  3_rlm.result.model.1.CT.B.txt' ; 3_rlm.result.model.1.txt'

# Sensitivity analyses
# White european ethnicity -> HELIX & GXII have participants removed, GenR & ALPSAC all whites
# EX:  3_rlm.result.cauc.model.1.CT.B.txt' 
# N <- c(1014,711,892,390) # cauc 
# Healthy diet -> model 4 

# Extra: HELIX subcohorts data 
# Prefixes: 'BIB','EDEN','INMA','RHEA','KANC','MOBA'
# N <- c(201,146,185,196,211,199,732,892,390) 
# artype <- c('450K', '450K', '450K', '450k', '450k', '450k') 


########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to metaanalysis folder -> DON'T 
# setwd("<path to metaanalysis folder>/metaanalysis")

# Files used in QC, needed in meta-analysis to plot ForestPlot
files <- c('results/preQC/HELIX/3_rlm.result.model.3.CT.B.txt', 
           'results/preQC/GXXI/3_rlm.result.model.3.CT.B.txt',
           'results/preQC/ALSPAC/3_rlm.result.model.3.CT.B.txt',
           'results/preQC/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Result folder
results_folder <- 'results/QC_EASIER/MODEL_·CT/' #if doesn't exist, it creates it

# Prefixes for each file
#prefixes <- c('HELIX', 'GXII', 'ALSPAC', 'GenR')
prefixes <- c('HELIX', 'GXXI','ALSPAC' , 'GenR')

# Exclude - MASK snp5 -> we have ethnicity exclusion in the ewas not need it here 
# ethnic <- c('EUR','GMAF1p', 'EUR', 'SAS', 'EUR', 'EAS', 'EUR', 'SAS', 'EUR', 'EUR', 'EUR', 'EAS')

# Exclude  
exclude <-c('control_probes','noncpg_probes','MASK_mapping','MASK_sub30_copy',
           'MASK_extBase','MASK_typeINextBaseSwitch','Unrel_450_EPIC_blood','Sex')

# Array type, used : EPIC or 450K 
artype <- c('450K','EPIC','450K','450K')
#artype <- c('450K','EPIC','450K')

N <- c(1138,732,892,390) 
#N <- c(1138,732,390)

#N <- c(1014,711,892,390)
n <- c(NA) 

# Minimum sample representation percentage required for CpGs
# Filter minimum percentage of missings for each CpG in cohort
# We need to define two parameters,
#  - colname_NforProbe:
#        Column name with Number of individuals per probe, this variable only needs
#           to be defined if you want to filter CpGs with low representation.
#         If defined value in colname_NforProbe not exists, no filter will be applied
#  - pcMissingSamples :
#        Máximum percent of missing samples allowed,

colname_NforProbe <- 'N_for_probe'
pcMissingSamples <- 0.9 # it is not used here but in the next step that is the meta-analysis 

# Venn diagrams
venn_diagrams <- list(
  #c('HELIX, GXII, ALSPAC','GenR'))
  c('HELIX, GXXI','GenR'))
  


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



## ###################### ##
##  QC - Quality Control  ##
## ###################### ##

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

