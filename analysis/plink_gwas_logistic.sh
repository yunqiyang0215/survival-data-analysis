#!/bin/bash

#SBATCH --job-name=gwas
#SBATCH --partition=tier2q
#SBATCH --time=24:00:00
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/gwas/gwas.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/gwas/gwas.err

# Command to submit the script 
#for CHR in {1..22}; do
#  sbatch your_script.sh $CHR
#done


# Define the regions of interest
declare -A regions
regions[1]="chr1_150600001_155100000"
regions[2]="chr2_102100001_105300000 chr2_143400001_147900000 chr2_236400001_242193529"
regions[6]="chr6_30500001_32100000"
regions[10]="chr10_6600001_12200000"
regions[11]="chr11_75500001_77400000"
regions[12]="chr12_46000001_48700000 chr12_54500001_56200000"
regions[15]="chr15_59000001_63400000"
regions[17]="chr17_33500001_39800000"
regions[19]="chr19_31900001_35100000"


# Receive chromosome number as a command line argument
CHR=$1


# GWAS all asthma
PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
pheno=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/pheno.txt
covar=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar.txt
idfile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id.txt


# Run GWAS for specified regions in the current chromosome
for region in ${regions[$CHR]}; do
  genofile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/${region}

  $PLINK2_EXEC --pfile $genofile --glm hide-covar no-x-sex omit-ref --1\
       --covar $covar --pheno $pheno --vif 100 \
       --threads 20 --memory 64000 \
       --out /gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic_out/all/gwas_${region}
done


# GWAS coa
PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
pheno=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/pheno_coa.txt
covar=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_coa.txt
idfile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id_coa.txt


# Run GWAS for specified regions in the current chromosome
for region in ${regions[$CHR]}; do
  genofile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/${region}

  $PLINK2_EXEC --pfile $genofile --glm hide-covar no-x-sex omit-ref --1\
       --covar $covar --pheno $pheno --vif 100 \
       --threads 20 --memory 64000 \
       --out /gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic_out/coa/gwas_${region}
done


# GWAS aoa
PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
pheno=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/pheno_aoa.txt
covar=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_aoa.txt
idfile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id_aoa.txt


# Run GWAS for specified regions in the current chromosome
for region in ${regions[$CHR]}; do
  genofile=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/${region}

  $PLINK2_EXEC --pfile $genofile --glm hide-covar no-x-sex omit-ref --1\
       --covar $covar --pheno $pheno --vif 100 \
       --threads 20 --memory 64000 \
       --out /gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic_out/aoa/gwas_${region}
done




