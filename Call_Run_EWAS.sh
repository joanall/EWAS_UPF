#! /bin/bash
#SBATCH --job-name=sensitivity
#SBATCH --partition=long
#SBATCH --mail-type=begin       
#SBATCH --mail-type=end         
#SBATCH --mail-user=joana.llaurado@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=140gb
#SBATCH --output=/PROJECTES/HELIX_OMICS/analyses/ATH_EWAS_UPF_JLLP/results/CLUSTER_RUN_ALLMODEL.txt

module purge > /dev/null 2>&1
module load lang/R


# Define array with model names

models_array=("model.1" "model.2" "model.3" "model.4" "model.1.CT" "model.2.CT" "model.3.CT" "model.4.CT")

i=$(($SLURM_ARRAY_TASK_ID -1))

model=${models_array[i]} 

# metadata="db/final_INMA_metadata.RDS" 
# methylation="db/final_INMA_bVals.RDS"
# output=".INMA.rds"

#R CMD BATCH 3_Build_Run_Toy_Model.R $model

Rscript script/2_Run_EWAS.R $model 

# Rscript/3_Run_Model.R $model $metadata $methylation $output
