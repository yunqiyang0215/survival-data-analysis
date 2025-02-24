#!/bin/bash

#SBATCH --job-name=get_geno_with_more_snps
#SBATCH --partition=tier2q
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G

# srun --time=24:00:00 --mem=96G -c 4 --partition=tier2q --pty /bin/bash

CHR=1
START=150600001
END=155100000
# CHR=2
# START=102100001
# END=105300000
# CHR=2
# START=143400001
# END=147900000
# CHR=2
# START=236400001
# END=242193529
# CHR=10
# START=6600001
# END=12200000
# CHR=11
# START=75500001
# END=77400000
# CHR=12
# START=46000001
# END=48700000
# CHR=12
# START=54500001
# END=56200000
# CHR=15
# START=59000001
# END=63400000
# CHR=17
# START=33500001
# END=39800000
# CHR=19
# START=31900001
# END=35100000

PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
BGEN_FILE=/gpfs/data/xhe-lab/uk-biobank-genotypes/ukb_imp_chr${CHR}_v3.bgen
GENO_SAMPLE_FILE=/gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample
ID_FILE=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id.txt
OUTFILE=chr${CHR}_${START}_${END}

$PLINK2_EXEC --sample $GENO_SAMPLE_FILE --bgen $BGEN_FILE \
  --chr $CHR --from-bp $START --to-bp $END --keep $ID_FILE \
  --maf 0.001 --snps-only --max-alleles 2 --rm-dup exclude-all \
  --make-pgen --threads 4 --memory 90000 --out $OUTFILE
