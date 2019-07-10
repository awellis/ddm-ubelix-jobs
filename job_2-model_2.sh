#!/bin/bash
#SBATCH --mail-user=andrew.ellis@psy.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="model_2_final"
#SBATCH --ntasks=4
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=2G

module load GSL
module load vital-it
module load R/3.5.1
Rscript job_2-model_2.R
