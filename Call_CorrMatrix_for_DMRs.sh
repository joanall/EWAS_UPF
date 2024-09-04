#! /bin/bash
#SBATCH --job-name=CORRMATRIX
#SBATCH --partition=long
#SBATCH --mail-type=begin       
#SBATCH --mail-type=end         
#SBATCH --mail-user=joana.llaurado@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=140gb
#SBATCH --output=/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_UPF_JLLP/results/CORR_MATRIX.txt

module purge > /dev/null 2>&1
module load lang/R



Rscript script/6_CorrMatrix_for_DMRs.R



