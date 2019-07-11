#!/bin/bash
#SBATCH --mail-user=andrew.ellis@psy.unibe.ch
#SBATCH --mail-type=end
#SBATCH --job-name="model_4_no_bias"
#SBATCH --ntasks=4
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=2G

module load GSL
module load vital-it
module load R/3.5.1
Rscript job_4-model_4.R
