#!/bin/bash

#SBATCH --job-name=get_geno_with_more_snps
#SBATCH --partition=tier2q
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=450G

module load gcc/12.1.0
module load R/4.3.1

Rscript coxphsusie_ukbb.R
