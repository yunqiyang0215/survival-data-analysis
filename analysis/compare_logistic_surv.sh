#!/bin/bash

#SBATCH --job-name=compare_result
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --mem=180gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/compare_res.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/compare_res.err

module load gcc/12.1.0
module load R/4.2.1


Rscript compare_logistic_surv.R


