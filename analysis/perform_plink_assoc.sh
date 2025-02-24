#!/bin/bash

#SBATCH --job-name=perform_link_assoc
#SBATCH --partition=tier2q
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G

# srun --time=2:00:00 --mem=96G -c 4 --partition=tier2q --pty /bin/bash

PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2

# Linear regression association analysis.
$PLINK2_EXEC --pfile ../data/geno_regions_more_snps/chr10_6600001_12200000 \
  --glm hide-covar no-x-sex omit-ref --1 --vif 100 \
  --threads 4 --memory 90000 \
  --pheno ../data/gwas_logistic/pheno_coa_std.txt \
  --covar ../data/gwas_logistic/covar_coa.txt \
  --out chr10_6600001_12200000

# Logistic regression association analysis.
$PLINK2_EXEC --pfile ../data/geno_regions_more_snps/chr10_6600001_12200000 \
  --glm hide-covar no-x-sex omit-ref --1 --vif 100 \
  --threads 4 --memory 90000 \
  --pheno ../data/gwas_logistic/pheno_coa.txt \
  --covar ../data/gwas_logistic/covar_coa.txt \
  --out chr10_6600001_12200000
