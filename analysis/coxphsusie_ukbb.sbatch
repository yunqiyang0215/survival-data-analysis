#!/bin/bash

#SBATCH --job-name=coxphsusie_ukbb
#SBATCH --partition=tier2q
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

# I increased the memory from 24 GB to 64 GB when using 10 CPUs.

module load gcc/12.1.0
module load R/4.3.1

Rscript coxphsusie_ukbb.R
