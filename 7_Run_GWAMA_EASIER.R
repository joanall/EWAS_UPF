# Author:  Dolors Pelegrí
# Title: Quality Control Script to use with EASIER package
# Date: 03/01/2023


################################################################################
#                      SET WORKING ENVIRONMENT                                 #
################################################################################

# Note: RUN IT IN INTERACTIVE SESSION AS GWAMA IS ISNTALLED IN CLUSTER
# PATH SHOUL DBE THE ROOT . OF THE PROJECT AS PATH ARE WRITTEN FROM THERE

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

# load package
library(EASIER)

## ######################################### ##
##  Meta-Analysis to use with EASIER package ##
## ######################################### ##


########## ----------  VARIABLES DEFINED BY USER  ----------  ##########

# Set working directory to metaanalysis folder
# setwd("<path to metaanalysis folder>/metaanalysis")

# Files used in QC, needed in meta-analysis to plot ForestPlot 

files <- c('results/preQC/HELIX/3_rlm.result.model.3.CT.B.txt', 
           'results/preQC/GXXI/3_rlm.result.model.3.CT.B.txt',
           'results/preQC/ALSPAC/3_rlm.result.model.3.CT.B.txt',
           'results/preQC/GenR/3_rlm.result.cauc.model.3.CT.txt')

# Prefixes for each file
prefixes <- c('HELIX', 'GXXI', 'ALSPAC', 'GenR')

# Samples in original files used in QC
# <- c(1014,711,892,390) #cauc
N <- c(1138,732,892,390)


# Define data for each meta-analysis
metafiles <- list(
  'Meta_Model_3CT' = c('HELIX', 'GXXI', 'ALSPAC', 'GenR'))



# Array type, used in each meta-analysis : EPIC, 450K or MIX (to be used when in meta-analyses we have 450K and EPIC arrays)
artype <- c('MIX') #


# Define maximum percent missing for each CpG (in the example they use 0.8)
#     if pcenMissin = 0 only runs meta-analysis with all data
# CpGs with precense lower than pcentMissing after GWAS meta-analysis will be deleted from the study.

# When running meta anlaysis with 4 cohorts, we'll consider those CpGs at least 3 out of 4 cohorts, so 75%.
# When running meta nalaysis sub helix, following same criteria 8/9 -> 0.667
pcentMissing <- 0.75 

# Paths with QCResults and path to store GWAMA results
results_folder <- 'results/QC_EASIER/MODEL_3CT'
results_gwama <- 'results/' # it creates from this path another folder named ' GWAMA_results'


# Venn diagrams ==> IMPORTANT : maximum 5 meta-analysis in each venn diagram
venndiag_threshold <- 0.05
venn_diagrams <- list(c("Model_3CT"))


########## ----------  END VARIABLES DEFINED BY USER  ----------  ##########



# GWAMA binary path  (GWAMA IsGlobal Server installation)
#gwama.dir <- paste0(Sys.getenv("HOME"), "/data/software/GWAMA/")
gwama.dir <- "/PROJECTES/HELIX_OMICS/software/GWAMA/"


## Create directory for GWAMA configuration files and GWAMA_Results
if(!dir.exists(file.path(paste(results_gwama, "GWAMA", sep="/") )))
  suppressWarnings(dir.create(file.path( paste(results_gwama, "GWAMA", sep="/"))))

## Create directory for GWAMA_Results
outputfolder <- paste0(results_gwama, "/GWAMA_Results")
if(!dir.exists(file.path( outputfolder )))
  suppressWarnings(dir.create(file.path( outputfolder)))


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
    plot_ForestPlot( dataPost, metafiles[[metf]], runs[j], inputfolder, names(metafiles)[metf], files, outputgwama, nsignificatives = 30  )
    
  }
  
}


# Venn_Diagrams for for meta-analysis with fixed effects

for (i in 1:length(venn_diagrams)){
  if(length(venn_diagrams[[i]])>1){
    plot_venndiagram(venn_diagrams[[i]], qcpath = outputfolder, plotpath =  paste0(results_gwama, "/GWAMA_Results"), pattern = '_Fixed_Modif.out',bn='Bonferroni', fdr='FDR', venndiag_threshold)
  }
}

if(dir.exists(file.path( paste(results_gwama, "GWAMA", sep="/") )))
  unlink(file.path(results_gwama, "GWAMA"), recursive=TRUE)