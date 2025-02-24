#!/bin/bash

#SBATCH --job-name=fit_susie_asthma
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=90:00:00
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_result/asthma/ukb.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_result/asthma/ukb.err


module load gcc/12.1.0
module load R/4.3.1
Rscript fit_susie.R $CHR


# sbatch --export=ALL,CHR=chr2 your_sbatch_script.sh

