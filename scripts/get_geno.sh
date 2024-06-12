#!/bin/bash

#SBATCH --job-name=get_geno
#SBATCH --partition=tier2q
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/gwas/get_geno.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/gwas/get_geno.err


# Command to submit the script 
#for CHR in {1..22}; do
#  sbatch your_script.sh $CHR
#done



# STEP1: filter SNPs with imputation score >0.6
# zcat /gpfs/data/xhe-lab/uk-biobank-genotypes/ukb_mfi_chr${CHR}_v3.txt.gz | awk '{if ($8 > 0.6) {print $2}}' | sort | uniq > /gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/raw/chr${CHR}.txt 

# Receive chromosome number as a command line argument
CHR=$1

PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
idfile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id.txt
genofile=/scratch/yunqiyang/ukb_asthma_self_report/data/genotypes_maf001_info6/chr${CHR}

# Get genotype 
$PLINK2_EXEC --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample \
       --bgen /gpfs/data/xhe-lab/uk-biobank-genotypes/ukb_imp_chr${CHR}_v3.bgen \
       --chr $CHR \
       --keep $idfile \
       --extract /gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/raw/chr${CHR}.txt \
       --maf 0.001 minor \
       --snps-only --max-alleles 2 --rm-dup exclude-all \
       --make-pgen --threads 20 --memory 64000 \
       --out $genofile

