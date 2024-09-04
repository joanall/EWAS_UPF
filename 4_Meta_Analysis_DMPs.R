# Title: Meta Analysis of EWAS results 
# Author:  JLLP (based on Dolors Pelegrí EASIER package)
# Date: 03/01/2023

# Description: Do a Meta Analysis of EWAS results using GWAMA software
# To use GWAMA it has to be installed in the cluster where analysis are running


#==============================================================================#
#                          CHOOSE MODEL TO META-ANALYZE                        #
#==============================================================================#

# options: model_1, model_2, model_3, healthy_diet, ethnicity,
# out_helix, out_alspac, out_gxxi, out_genr

which_model <- 'model_3'

#==============================================================================#
#                           SET WORKING ENVIRONMENT                            #
#==============================================================================#

# load package
library(EASIER)

#===============================================================================


#==============================================================================#
#                                                                              #
#                       META ANALYSIS PARAMETERS                               #----
#                                                                              #
#==============================================================================#

#------------------------------------------------------------------------------#
#                                 MODEL 1                                      #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'model_1') {
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.1.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.1.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.1.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.1.CT.txt')
  
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(1138,732,892,390)
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Model_1' = c('HELIX', 'GXXI', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/MODEL_1'

}

#------------------------------------------------------------------------------#
#                                 MODEL 2                                      #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'model_2') {
  
  # Input files to perform QC 
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.2.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.2.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.2.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.2.CT.txt')
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(1138,732,892,390)
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Model_2' = c('HELIX', 'GXXI', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/MODEL_2'

}

#------------------------------------------------------------------------------#
#                                 MODEL 3                                      #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'model_3') {
  
  # Input files to perform QC 
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(1138,732,892,390)
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Model_3' = c('HELIX', 'GXXI', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/MODEL_3'

}



#------------------------------------------------------------------------------#
#                   SENSITIVITY 1 (healthy diet)                               #----                                      
#------------------------------------------------------------------------------#


if (which_model == 'healthy_diet') {
  
  # Input files to perform QC 
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.4.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.4.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.4.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.4.CT.txt')
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(1138,732,892,390)
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Model_Healthy_Diet' = c('HELIX', 'GXXI', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/HEALTHY_DIET'

}

#------------------------------------------------------------------------------#
#                       SENSITIVITY 2 (ethnicity)                              #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'ethnicity') {
  
  # Input files to perform QC 
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.cauc.model.3.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.cauc.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')
  
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(1014,711,892,390) 
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Model_Ethnicity' = c('HELIX', 'GXXI', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/ETHNICITY'

}



#------------------------------------------------------------------------------#
#                             Out HELIX                                        #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'out_helix') {
  
  # Input files to perform QC 
  files <- c('new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')
  
  
  # Prefixes for each file
  prefixes <- c('GXXI', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(732,892,390) # sample size x cohort
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Out_HELIX' = c('GXXI', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_HELIX'

}



#------------------------------------------------------------------------------#
#                              Out GXXI                                        #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'out_gxxi') {
  
  # Input files to perform QC 
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'ALSPAC', 'GenR')
  
  # Samples in original files used in QC
  N <- c(732,892,390) # sample size x cohort
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Out_GXXI' = c('HELIX', 'ALSPAC', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_GXXI'
}

#------------------------------------------------------------------------------#
#                             Out ALSPAC                                       #----                                      
#------------------------------------------------------------------------------#


if (which_model == 'out_alspac') {
  
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/GenR/3_rlm.result.cauc.model.3.CT.txt')
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'GenR')
  
  # Samples in original files used in QC
  N <- c(1138,732,390) # sample size x cohort
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Out_ALSPAC' = c('HELIX', 'GXXI', 'GenR'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_ALSPAC'
  
}



#------------------------------------------------------------------------------#
#                             Outt Gen R                                       #----                                      
#------------------------------------------------------------------------------#

if (which_model == 'out_genr') {
  
  files <- c('new_version_14_02_2024/results/QC_28082024/HELIX/3_rlm.result.model.3.CT.B.txt', 
             'new_version_14_02_2024/results/QC_28082024/GXXI/3_rlm.result.cauc.model.3.CT.B.txt',
             'new_version_14_02_2024/results/QC_28082024/ALSPAC/3_rlm.result.model.3.CT.B.txt')
  
  # Prefixes for each file
  prefixes <- c('HELIX', 'GXXI', 'ALSPAC')
  
  # Samples in original files used in QC
  N <- c(1138,732,892) # sample size x cohort
  
  # Define data for each meta-analysis
  metafiles <- list('Meta_Out_GenR' = c('HELIX', 'GXXI', 'ALSPAC'))
  
  # Array type, used in each meta-analysis : EPIC, 450K or MIX 
  artype <- c('MIX') #
  
  # Define maximum percent missing for each CpG 
  pcentMissing <- 0.75 #CpGs at least 3 out of 4 cohorts
  
  # Paths with QCResults and path to store GWAMA results
  results_folder <- 'new_version_14_02_2024/results/QC_28082024/Out_GenR'

}


#==============================================================================#
#                                                                              #
#                             DO META ANALYSIS                                 #----
#                                                                              #
#==============================================================================#



#------------------------------------------------------------------------------#
#                               RUN GWAMA                                      #----                                      
#------------------------------------------------------------------------------#

#GWAMA binary path  (GWAMA IsGlobal Server installation)
#gwama.dir <- paste0(Sys.getenv("HOME"), "/data/software/GWAMA/")
gwama.dir <- "/PROJECTES/HELIX_OMICS/software/GWAMA/"
results_gwama <- 'new_version_14_02_2024/results/' # it creates from this path another folder named ' GWAMA_results'

## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(paste(results_gwama, "GWAMA", sep="/") )))
  suppressWarnings(dir.create(file.path( paste(results_gwama, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_gwama, "/GWAMA_Results")
if(!dir.exists(file.path( outputfolder )))
  suppressWarnings(dir.create(file.path(outputfolder)))


# Create hapmap files for the different artypes that we cab use (450K and EPIC)
# map file is used in Manhattan plots
hapmapfile_450K <- paste(results_gwama,"GWAMA", "hapmap_450K.map" ,sep = "/")
generate_hapmap_file("450K", hapmapfile_450K)
hapmapfile_EPIC <- paste(results_gwama,"GWAMA", "hapmap_EPIC.map" ,sep = "/")
generate_hapmap_file("EPIC", hapmapfile_EPIC)
hapmapfile_MIX <- paste(results_gwama,"GWAMA", "hapmap_MIX.map" ,sep = "/")
generate_hapmap_file("MIX", hapmapfile_MIX)



for( metf in 1:length(metafiles))
{
  
  list.lowCpGs <- NULL
  
  # Create folder for a meta-analysis in GWAMA folder, here we store the GWAMA input files for each meta-analysis,
  # We create one for complete meta-analysis
  if(!dir.exists(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf] ,sep="/") )))
    suppressWarnings(dir.create(file.path( paste(results_gwama,"GWAMA", names(metafiles)[metf], sep="/"))))
  # We create another for meta-analysis without filtered CpGs with low percentage (sufix _Filtr)
  if(!dir.exists(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr") )))
    suppressWarnings(dir.create(file.path( paste0(results_gwama,"/GWAMA/", names(metafiles)[metf],"_Filtr"))))
  
  # GWAMA File name base
  inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf])
  
  modelfiles <- unlist(metafiles[metf])
  
  runs <- c('Normal', 'lowcpgs') # Execution with all CpGs and without filtered CpGs
  lowCpGs = FALSE;
  outputfiles <- list()
  
  outputgwama <- paste(outputfolder,names(metafiles)[metf],sep = '/')
  
  for(j in 1:length(runs))
  {
    if(runs[j]=='lowcpgs') {
      lowCpGs = TRUE
      # Get low presence CpGs in order to exclude this from the new meta-analysis
      list.lowCpGs <- get_low_presence_CpGs(outputfiles[[j-1]], pcentMissing)
      inputfolder <- paste0(results_gwama,"/GWAMA/",  names(metafiles)[metf], "_Filtr")
      outputgwama <- paste0(outputgwama,"_Filtr")
    }
    
    # Create GWAMA files for each file in meta-analysis and execute GWAMA
    for ( i in 1:length(modelfiles) )
      create_GWAMA_files(file.path(results_folder,modelfiles[i]),  modelfiles[i], inputfolder, N[which(prefixes==modelfiles[i])], list.lowCpGs )
    
    
    # Get hapmapfile attending to current metaanalysis artype
    hapmapfile <- hapmapfile_450K
    if(artype[metf]=='EPIC'){
      hapmapfile <- hapmapfile_EPIC
    } else if(artype[metf]=='MIX'){
      hapmapfile <- hapmapfile_MIX
    }
    
    #.Original.#outputfiles[[runs[j]]] <- execute_GWAMA_MetaAnalysis(prefixgwama, names(metafiles)[metf])
    outputfiles[[runs[j]]] <- run_GWAMA_MetaAnalysis(inputfolder, outputgwama, names(metafiles)[metf], gwama.dir, hapmapfile)
    
    # Post-metha-analysis QC --- >>> adds BN and FDR adjustment
    dataPost <- get_descriptives_postGWAMA(outputgwama, outputfiles[[runs[j]]], modelfiles, names(metafiles)[metf], artype[metf], N[which(prefixes %in% modelfiles)] )
    
    # Forest-Plot
    # plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf], files, outputgwama, nsignificatives = 30  )
  
  }
  
}


#===============================================================================















