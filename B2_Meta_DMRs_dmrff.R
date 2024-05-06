# Author:  JLLP
# Title: MEta analaysis DMRs detected with dmrff 
# Date: 03/01/2023

#Note: Run in yamabuki interactive section. Here no aviabale version of dmrff package

################################################################################
#                      SET WORKING ENVIRONMENT                                 #
################################################################################

# install.packages("remotes")
#remotes::install_github("perishky/dmrff")

library(dmrff)
library(readr)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(bumphunter)
library(data.table)
library(GenomicRanges)

# Load the pre-dmrff object 
dmr_HELIX <- readRDS(file="new_version_14_02_2024/results/DMRs/dmrff/CorrMatrix_HELIX_model.3.CT.rds")
dmr_GXXI <- readRDS(file="new_version_14_02_2024/results/DMRs/dmrff/CorrMatrix_GXXI_model.3.CT.rds")
dmr_ALSPAC <- readRDS(file="new_version_14_02_2024/results/DMRs/dmrff/CorrMatrix_ALSPAC_model.3.CT.rds")
dmr_GENR <- readRDS(file="new_version_14_02_2024/results/DMRs/dmrff/5_CorrMatrix_GenRmodel.3.CT.rds")


# -----------------------------------------------------------------------------#
#                             META ANALYSIS DMRFF                              #    
# -----------------------------------------------------------------------------#

print("Meta dmr default")

meta <- dmrff.meta(list(dmr_ALSPAC, dmr_HELIX ,dmr_GXXI, dmr_GENR), 
                   verbose= T,
                   p.cutoff= 1E-05, # Unadjusted p-value cutoff for membership in a candidate DMR
                   maxgap = 500)



# Note canidates: here we havent filtered either for p value or number cpgs per region 
# The output  of dmrff.meta() has 2 objects in it:
#       -   $dmrsdata frame of the meta-analysed regions
#       -   $ewas  EWAS meta-analysis on which it was based.
# [dmrff.candidates] Found  12267  candidate regions.

# Raw output (a lot of regions have only 1 Cpg)
readr::write_rds(meta,"new_version_14_02_2024/results/DMRs/dmrff/dmrff_meta_4c.rds")

# -----------------------------------------------------------------------------#
#                          FILTERING DETECTED DMRs                             #    
# -----------------------------------------------------------------------------#

# > head(metaresults_default$ewas)
# estimate         se           z   p.value chr    pos
# cg18147296  0.021473879 0.02594520  0.82766287 0.4078615   1 812539
# cg13938959 -0.010854420 0.02725891 -0.39819711 0.6904849   1 834183

head(meta$dmrs[which(meta$dmrs$p.value <= 1E-05), ])

# Keep only those regions with more than 1 CpG & p adjusted < 0.05
dmrs <- meta$dmrs[which(meta$dmrs$p.adjust <= 0.05 & meta$dmrs$n >= 2), ]
print(nrow(dmrs)) #0 
write_rds(dmrs, "results/DMRs/dmrff/suggestive_dmrs_fdr.rds")

# Keep only those regions with more than 1 CpG & p value (suggestive) < 1E-05
dmrs <- meta$dmrs[which(meta$dmrs$p.value <= 1E-05 & meta$dmrs$n >= 2), ]
write_rds(dmrs, "results/DMRs/dmrff/suggestive_dmrs_1E05.rds")
# data.table::fwrite(dmrs, file = "results/CorrMatrix_DMRs/metadmr.model.4.CT.default.v2.csv",sep=";")


# -----------------------------------------------------------------------------#
#                        MAP DMRs TO GENES                                     #    
# -----------------------------------------------------------------------------#

# Match gene to region 
dmrs$seqnames <- paste0("chr", as.character(dmrs$chr))

dmrs.GR <- GenomicRanges::makeGRangesFromDataFrame(dmrs,
                                    keep.extra.columns=TRUE, 
                                    start.field="start",
                                    end.field="end", 
                                    seqnames.field="seqnames")

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genome <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene) #annotation used to map genes, no tocar

dmrs_genes <- bumphunter::matchGenes(dmrs.GR,
                          genome, 
                          type = "any",  
                          skipExons = TRUE,  #+ rapid if exons skip
                          verbose = TRUE)


dmrs$gene <- dmrs_genes$name
write_rds(dmrs, "results/DMRs/dmrff/genemap_suggestive_dmrs_1E05.rds")

# -----------------------------------------------------------------------------#
#                          MAP DMRs TO CPGs                                    #    
# -----------------------------------------------------------------------------#

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
write_rds(dmr_results_completed, "results/DMRs/dmrff/cpgmap_suggestive_dmrs_1E05.rds")
write_csv(dmr_results_completed, "results/DMRs/dmrff/cpgmap_suggestive_dmrs_1E05.csv")













