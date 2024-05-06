#! /bin/bash
#SBATCH --job-name=dmrs
#SBATCH --partition=long
#SBATCH --mail-type=begin       
#SBATCH --mail-type=end         
#SBATCH --mail-user=joana.llaurado@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=140gb
#SBATCH --output=/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_UPF_JLLP/results/meta_dmrs.txt

module purge > /dev/null 2>&1
module load lang/R

#R CMD BATCH 3_Build_Run_Toy_Model.R $model

Rscript script/5B_Meta_dmrff.R
