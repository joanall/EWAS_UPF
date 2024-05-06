# Created: Code created by M. Suderman.
# Adapted:  D. Caramaschi and Andrea Cortes & Sara Sammallahti on March 1st 2019
# Updated:  by Silvia Alemany and Jolien Rijlaarsdam on January 15th 2020


# Title: Creating a correlation matrix for DMRFF
# Date: 13/01/2023

# Use: Calculates the correlations between each CpG site and the next 20 CpG 
# sites in the genome. The resulting information cannot be used to reveal any 
# information about individual samples in the dataset so it is safe to share 
# with external users.

# It will be necessary for performing a meta-analysis of DMRs while taking into 
# account dependencies between CpG sites in order to avoid inflating association statistics.


# Note: Run it with main model = model 3 ct

################################################################################
#                          SET WORKING ENVIRONMENT                             #
################################################################################

# Load libraries - 

library(matrixStats)
library(parallel)
#library(meffil) #to install meffil https://github.com/perishky/meffil/wiki/Installation
# library(data.table)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")     
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(dplyr)
library(tidyverse)
library(minfi)
# library(EASIER)


################################################################################
#                       DEFINE PARAMETERS FROM USER                            #
################################################################################

# Define name of cohort 

cohort <- "HELIX" 

rlm_results_directory <- "results/EWAS_RESULTS/HELIX/3_rlm.result.model.3.CT.B.rds" 

methylation_data_directory <- "db/" 

#  Define output/results directory to save output files

output_directory <- "results/CorrMatrix_DMRs"

# Check if directory exists, if not create it

if (!dir.exists(output_directory)) { dir.create(output_directory)
  
} else { print("Directory already exists")}

#------------------------------------------------------------------------------#
#           DEFINE FUNCTIONS TO  BUILT CORRELATION MATRIX FOR dmrff            #
#------------------------------------------------------------------------------#

dmrff.pre <- function(estimate, se, methylation, chr, pos, maxsize=20, verbose=T) {
  stopifnot(is.vector(estimate))
  stopifnot(is.vector(se))
  stopifnot(is.matrix(methylation))
  stopifnot(is.vector(chr))
  stopifnot(is.vector(pos))
  stopifnot(length(estimate) == length(se))
  stopifnot(length(estimate) == nrow(methylation))
  stopifnot(length(estimate) == length(chr))
  stopifnot(length(estimate) == length(pos))
  
  # sort input by chromosomal position
  idx <- order(chr,pos)
  sorted <- identical(idx, 1:length(idx))
  if (!sorted) {
    estimate <- estimate[idx]
    se <- se[idx]
    chr <- chr[idx]
    pos <- pos[idx]
    methylation <- methylation[idx,,drop=F]
  }
  
  ## calculate rho (CpG correlation matrix)
  m <- methylation - rowMeans(methylation, na.rm=T)
  mm <- rowSums(m^2, na.rm=T)
  ss <- sqrt(mm)
  rho <- do.call(cbind, parallel::mclapply(1:maxsize, function(size) {
    sapply(1:length(chr), function(i) {
      if (i + size <= length(idx)) {
        numer <- sum(m[i,] * m[i+size,])
        denom <- ss[i] * ss[i+size]
        numer/denom
      }
      else NA
    })
  }))
  n <- rowSums(!is.na(m))
  sd <- sqrt(mm/(n-1))
  
  sites <- paste(chr, pos, sep=":")
  list(sites=sites,
       chr=chr,
       pos=pos,
       estimate=estimate,
       se=se,
       rho=rho,
       sd=sd)
}


#------------------------------------------------------------------------------#
#                     LOAD EWAS RESULTS AND METHYLATION DATA                   #
#------------------------------------------------------------------------------#

# Load methylation data (name: meth_vals)
load(methylation_data_directory) # 1138 386518

# Load results after running RLM model 
rlm.result <- readRDS(rlm_results_directory)
rlm.result$probeid <- rownames(rlm.result)

#  Rows are CpG sites & columns are samples. Transpose data. 
methylation <- t(meth_vals)
rm(meth_vals)

# Load annotation 
annotation450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Create annotation dataset 
annotation <- annotation450k %>% 
  as.data.frame() %>% 
  droplevels()  %>%
  dplyr::mutate(CPGs=rownames(annotation450k))

annotation$chr <- gsub("chr","", as.character(annotation$chr))
annotation$chr <- as.integer(annotation$chr)

rm(annotation450k)


#------------------------------------------------------------------------------#
#                            MATCH METHYLATION AND EWAS DATA                   #
#------------------------------------------------------------------------------#

# Match rlm results - meth data (ewas_res = rlm.result )

print('start match data')
methylation <- methylation[na.omit(match(rlm.result$probeid,rownames(methylation))),] 
rlm.result_match <- rlm.result[na.omit(match(rownames(methylation), rlm.result$probeid)),] 
ifelse(all(rlm.result_match$probeid==rownames(methylation)), 
       "meth and results data successfully matched :) ","Data not matched :(")

# Match  meth data - annotation 

annotation <- annotation[na.omit(match(rownames(methylation), annotation$Name)),]
methylation <- methylation[na.omit(match(annotation$Name,rownames(methylation))),]
ifelse(all(annotation$name==rownames(methylation)), "meth and annotation data successfully matched :) ","Data not matched :(")


#------------------------------------------------------------------------------#
#                            MATCH METHYLATION AND EWAS DATA                   #
#------------------------------------------------------------------------------#

print('start coorr')

# Run dmrff.pre function -> this step may take some time
pre <- dmrff.pre(estimate = rlm.result_match$estimate, #coef for each CpG site
                 se = rlm.result_match$std.error, 
                 methylation = methylation, 
                 chr = annotation$chr, 
                 pos = annotation$pos)  

# Save matrix
write_rds(pre, file = paste0(output_directory,"/5_CorrMatrix_", cohort, "model.1.B.rds"))



