# Title: Enrichment of DMPs
# Author:  JLLP 
# Date: 13/01/2023
# Description: Perform Functional Enrichment Analysis of the selected CpGs based on the
# specified threshold in order to gain biological insight 



#==============================================================================#
#                           DEFINE USER PARAMETERS                             #
#==============================================================================#


# Set working directory to enrichment folder
#setwd("<path to metaanalysis folder>/Enrichment")

# Files with CpG data to enrich may be a CpGs list or annotated GWAMA output
filename <- c('new_version_14_02_2024/results/ANALYZE_DMPs_28082024/MODEL_3/Meta_Model_3_Fixed_Modif.rds')

cutoff_enrich <- 1E-04

# Array type, used : EPIC or 450K
artype <- '450K'

# Result paths definition for QC, Meta-Analysis and Enrichment
output_folder <- 'new_version_14_02_2024/results/ENRICHMENT_DMPs/'


if (!dir.exists(output_folder)) {dir.create(output_folder, recursive = TRUE)
} else { print("Directory already exists")}



#==============================================================================#
#                           SET WORKING ENVIRONMENT                            #----
#==============================================================================#

# load package
library(EASIER)
#library(missMethyl)
#install.packages("devtools")
#library(devtools)
#install_github("isglobal-brge/brgeEnrich")

library(brgeEnrich) # Parse web Consensus Path
library(msigdbr) # Load collection Molec.Sign DB
library(dplyr)

#==============================================================================#
#                              FUNCTIONAL ENRICHMENT                           #----
#==============================================================================#

#------------------------------------------------------------------------------#
#                               Prepare the data                               #----    
#------------------------------------------------------------------------------#

results <- readRDS(filename) 

enrich_cpgs <- results[results$p.value <= cutoff_enrich, c('CpGId','UCSC_RefGene_Name')]

save_to <- paste0(output_folder, '/', sub(".rds", 
                                          paste0('_genes_enrichment',cutoff_enrich,".csv"), 
                                          basename(filename)))

write.table(enrich_cpgs, save_to, sep=";")


#------------------------------------------------------------------------------#
#                                 GO & KEGG                                    #----    
#------------------------------------------------------------------------------#
  
# Enrichment with missMethyl - GO and KEGG 


# Get all the CpG sites used in the analysis to form the background (gene universe)
allcpgs <- results$CpGId    # Alternative: allcpgs <- rownames(annotation)

# Get CpGs to enrich 
enrich_cpgs <- results[results$p.value <= cutoff_enrich, c('CpGId')]

# Run GO enrichment
go_enrich_result <- tryCatch({
  
  enrich <-  missMethyl::gometh(sig.cpg = enrich_cpgs, all.cpg = allcpgs, collection = "GO",
                                plot.bias = TRUE, prior.prob = TRUE)
  save_to <- paste0(output_folder, '/', sub(".rds", 
                                             paste0('_go_enrichment', cutoff_enrich, ".csv"), 
                                             basename(filename)))
  write.table(enrich, save_to, sep=";")
  
}, error = function(e) {
  message("Error in GO enrichment for file ", basename(filename), ": ", e$message)
})


# Run KEGG enrichment
kegg_enrich_result <- tryCatch({
  enrich <- missMethyl::gometh(sig.cpg = enrich_cpgs, all.cpg = allcpgs, collection = "KEGG",
                               plot.bias = FALSE, prior.prob = TRUE)
  
  save_to <- paste0(output_folder, '/', sub(".rds", 
                                             paste0('_kegg_enrichment', cutoff_enrich, ".csv"), 
                                             basename(filename)))
  write.table(enrich, save_to, sep=";")
  
}, error = function(e) {
  message("Error in KEGG enrichment for file ", basename(filename), ": ", e$message)
})



#------------------------------------------------------------------------------#
#                           ConsensusPathDB                                    #----    
#------------------------------------------------------------------------------#

# Consensus path http://cpdb.molgen.mpg.de/ (gene-set analysis – over-representation analysis)

# Available FSet types :
# 1 P     manually curated pathways from pathway databases
# 2 N     interaction network neighborhood-based functional sets
# 3 G2    Gene Ontology-based sets, GO level 2
# 4 G3    Gene Ontology-based sets, GO level 3
# 5 G4    Gene Ontology-based sets, GO level 4
# 6 G5    Gene Ontology-based sets, GO level 5
# 7 C     protein complex-based sets
  
# Select Genes to be enriched 
enrich_genes <- results[results$p.value <= cutoff_enrich, "UCSC_RefGene_Name"]
enrich_genes <- unique(unlist(sapply(enrich_genes, function(x) strsplit(x, ";"))))
enrich_genes <- na.omit(as.character(enrich_genes))
enrich_genes <- unlist(lapply(enrich_genes, function(gene_name) paste0(gene_name, '_HUMAN')))

# ORA for protein complex-based set (C)
CPDB_enrich <- cpdbOverRepresentationAnalysis(entityType = 'genes', 
                                              fsetType = 'C',
                                              accNumbers = enrich_genes, 
                                              accType = 'uniprot')

# Write data to a file
save_to <- paste0(output_folder, '/', sub(".rds", 
                                           paste0("CPDB_enrichment_ProteinComplex",cutoff_enrich,".csv"), 
                                           basename(file)))

write.table(CPDB_enrich, save_to, sep=";")

# ORA for manually curated pathways from pathway databases (P)  
CPDB_enrich <- cpdbOverRepresentationAnalysis(entityType = 'genes', 
                                              fsetType = 'P',
                                              accNumbers = enrich_genes, 
                                              accType = 'uniprot')


# Write data to a file
save_to <- paste0(output_folder, '/', sub(".rds", 
                                           paste0("CPDB_enrichment_Pathways",cutoff_enrich,".csv"), 
                                           basename(file)))

write.table(CPDB_enrich, save_to, sep=";")



#------------------------------------------------------------------------------#
#                   Molecular Signatures Database enrichment                   #----
#------------------------------------------------------------------------------#

# Load the 'C2' collection for Homo sapiens
Hs.c2 <- msigdbr(species = "Homo sapiens", category = "C2")

# Transform the data frame to a list format as required by gsameth
Hs.c2_list <- split(Hs.c2$entrez_gene, Hs.c2$gs_name)

# Get all the CpG sites used in the analysis to form the background (gene universe)
allcpgs <- results$CpGId    # Alternative: allcpgs <- rownames(annotation)

# Get CpGs to enrich 
enrich_cpgs <- results[results$p.value <= cutoff_enrich, 'CpGId']

msd_enrich <- missMethyl::gsameth(sig.cpg = enrich_cpgs,
                                  all.cpg = allcpgs,
                                  collection = Hs.c2_list,
                                  array.type = artype
)

# Write data to a file
save_to <- paste0(output_folder, '/', sub(".rds", 
                                           paste0("MSigDB_enrichment_",cutoff_enrich,".csv"), 
                                           basename(filename)))

write.table(msd_enrich, save_to, sep=";")













    
