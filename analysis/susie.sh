#!/bin/bash

#SBATCH --job-name=fit_susie_asthma
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=100:00:00
#SBATCH --mem=450gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/susie.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/susie.err


module load gcc/12.1.0
module load R/4.3.1

# Specify the region as a command-line argument
region="$1"

#Rscript susie.all.R "$region"
#Rscript susie.aoa.R "$region"
Rscript susie.coa.R "$region"
