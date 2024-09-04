# Title: Meta Analysis of DMRs
# Author:  JLLP 
# Date: 13/01/2023
# Description: Perform a Meta Analysis by combining the correlation matrix of DMRs
# of the different cohorts after running EWAS for DMPs. 



#==============================================================================#
#                           DEFINE USER PARAMETERS                             #
#==============================================================================#

# Load the pre-dmrff object 
corrmat_helix <- "new_version_14_02_2024/results/DMRs/dmrff/CorrMatrix_HELIX_model.3.CT.rds"
corrmat_GXXI <- "new_version_14_02_2024/results/DMRs/dmrff/CorrMatrix_GXXI_model.3.CT.rds"
corrmat_ALSPAC <- "new_version_14_02_2024/results/DMRs/dmrff/CorrMatrix_ALSPAC_model.3.CT.rds"
corrmat_GENR <- "new_version_14_02_2024/results/DMRs/dmrff/5_CorrMatrix_GenRmodel.3.CT.rds"

output_folder <- 'new_version_14_02_2024/results/ANALYZE_DMRs_28082024/MODEL_3/'
model <- 'MODEL_3'

if (!dir.exists(output_folder)) {dir.create(output_folder, recursive = TRUE)
} else { print("Directory already exists")}


#==============================================================================#
#                           SET WORKING ENVIRONMENT                            #----
#==============================================================================#

# install.packages("remotes")
#remotes::install_github("perishky/dmrff")

library(dmrff)
library(readr)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(bumphunter)
library(data.table)
library(GenomicRanges)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

#==============================================================================#
#                           META ANALYSIS DMRFF                                #----
#==============================================================================#

dmr_HELIX <- readRDS(corrmat_helix)
dmr_GXXI <- readRDS(corrmat_GXX)
dmr_ALSPAC <- readRDS(corrmat_ALSPAC)
dmr_GENR <- readRDS(corrmat_GENR)


meta <- dmrff.meta(list(dmr_ALSPAC, dmr_HELIX ,dmr_GXXI, dmr_GENR), 
                   verbose= T,
                   p.cutoff= 1E-05, # Unadjusted p-value cutoff for membership in a candidate DMR
                   maxgap = 500)



# Note canidates: here we havent filtered either for p value or number cpgs per region 
# The output  of dmrff.meta() has 2 objects in it:
#       -   $dmrsdata frame of the meta-analysed regions
#       -   $ewas  EWAS meta-analysis on which it was based.
# [dmrff.candidates] Found  12267  candidate regions.


#------------------------------------------------------------------------------#
#                          FILTERING DETECTED DMRs                             #----    
#------------------------------------------------------------------------------#

# Keep only those regions with more than 1 CpG & p adjusted < 0.05
dmrs <- meta$dmrs[which(meta$dmrs$p.adjust <= 0.05 & meta$dmrs$n >= 2), ]
write_rds(dmrs, paste0(output_folder, model,'_DMRs_FDR.rds'))

# Keep only those regions with more than 1 CpG & p value (suggestive) < 1E-05
dmrs <- meta$dmrs[which(meta$dmrs$p.value <= 1E-05 & meta$dmrs$n >= 2), ]
write_rds(dmrs, paste0(output_folder, model, '_DMRs_1E05.rds'))


#------------------------------------------------------------------------------#
#                        MAP DMRs TO GENES                                     #----    
#------------------------------------------------------------------------------#

# Match gene to region 
dmrs$seqnames <- paste0("chr", as.character(dmrs$chr))

dmrs.GR <- GenomicRanges::makeGRangesFromDataFrame(dmrs,
                                    keep.extra.columns=TRUE, 
                                    start.field="start",
                                    end.field="end", 
                                    seqnames.field="seqnames")

genome <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene) #annotation used to map genes, no tocar

dmrs_genes <- bumphunter::matchGenes(dmrs.GR,
                          genome, 
                          type = "any",  
                          skipExons = TRUE,  #+ rapid if exons skip
                          verbose = TRUE)


dmrs$gene <- dmrs_genes$name
write_rds(dmrs, paste0(output_folder, model,'_genemap_DMRs_1E05.rds'))

#------------------------------------------------------------------------------#
#                          MAP DMRs TO CPGs                                    #----   
#------------------------------------------------------------------------------#

# Load annotation
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) #485512     33
ann <- data.table(as.data.frame(ann))
ann$chr <- gsub("chr","", as.character(ann$chr))

# Keep only those chromosomes for which we have found a region  (ex. we exclude X and Y chr)
cpgs_to_join <- ann[ann$chr %in% dmrs$chr , c("chr", "pos", "Name", "UCSC_RefGene_Name")] # 473864


# Note: Taken into account that in each chromosome positions restart 

# List subset by chromosome 
subset_list_ann <- split(cpgs_to_join, cpgs_to_join$chr)
dmrs <- data.table(dmrs) # !!! L Compte hem de passar a data-table igual que en ann la manera d'iterar segons data.frmae i data.table es diferent!!!!!
subset_list_dmr <- split(dmrs, dmrs$chr)
N_chr <- length(unique(cpgs_to_join$chr)) # Get number chromosomes


for (n in 1:N_chr){ # go through each sublist of dmr /cpgs corresponding to  chr
  
  cat("Regions of Chromosome", n)
  #temp1 <- toy_subset_list_dmr[[n]] #ex: temp1 is the first dataframe of the sublist of dmr, with the dmrs found in chr1
  temp1 <- subset_list_dmr[[n]] # ex: here we only have those regions found in chr 1
  temp2 <- subset_list_ann[[n]] #ex: here we onlu have the anntotaiton for those cpgs in chr1
  
  cpgs_within_chr <- list()
  
  # go through each row (=region) in selected sublist of chr x
  for (i in 1:nrow(temp1)) { 
    start <- temp1[i]$start
    end <- temp1[i]$end
    
    cpgs_within_region <- list()
    #cat("Region",i)
    
    # go through each row (=annotated cpg) of chr x 
    for (j in 1:nrow(temp2)){
      pos <- temp2[j]$pos
      is_within_region <- data.table::between(pos,start,end) #check if cpg is in region x by overlapping positions
      #length(is_within_region) > 0 && is_within_region
      if (is_within_region == "TRUE") {
        cpg <- temp2[j]$Name
        cpgs_within_region <- unlist(append(cpgs_within_region, cpg))
      }
      
    }
    cat("Region", i, "of chr", n, "\n")
    #print(cpgs_within_region)
    cpgs_within_chr <- c(cpgs_within_chr, list(cpgs_within_region))
    
  }
  #print(cpgs_within_chr)
  temp1$cpgs <- cpgs_within_chr
  subset_list_dmr[[n]] <- temp1
  #subset_list_dmr[[n]] <- temp1
}


dmr_results_completed <- rbindlist(subset_list_dmr, fill=TRUE)

write_rds(dmr_results_completed, paste0(output_folder, model,'_cpgmap_DMRs_1E05.rds'))
write_csv(dmr_results_completed, paste0(output_folder, model,'_cpgmap_DMRs_1E05.csv'))















