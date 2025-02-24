#!/bin/bash

#SBATCH --job-name=filter_ukb
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb-%j.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb-%j.err

PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
ID=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id.txt
GENO_SAMPLE_FILE=/gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample
PLINK_OUTPUT=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions

# To run this file:
# for i in {1..22}; do
#    sbatch your_script.sh $i
# done

# Define regions: COA specific
declare -a Regions=(
"chr3:188200001-192600000"
"chr4:35800001-41200000"
"chr18:61300001-63900000"
)


# Receive chromosome number as a command line argument
CHROMOSOME=$1

echo $CHROMOSOME
 
# Process each region
for region in "${Regions[@]}";do
    chr=$(echo $region | cut -d':' -f1 | sed 's/chr//')
    if [[ "$chr" -eq "$CHROMOSOME" ]]; then
        start=$(echo $region | cut -d':' -f2 | cut -d'-' -f1)
        end=$(echo $region | cut -d':' -f2 | cut -d'-' -f2)
        BGEN_FILE=/gpfs/data/xhe-lab/uk-biobank-genotypes/ukb_imp_chr${CHROMOSOME}_v3.bgen
        OUTPUT=${PLINK_OUTPUT}/chr${CHROMOSOME}_${start}_${end}
        
        echo "Running PLINK2 for $region"        
        $PLINK2_EXEC --sample $GENO_SAMPLE_FILE \
        --bgen $BGEN_FILE \
        --chr $CHROMOSOME \
        --from-bp $start \
        --to-bp $end \
        --keep $ID \
        --maf 0.01 \
        --snps-only \
        --max-alleles 2 \
        --rm-dup exclude-all \
        --make-pgen \
        --threads 8 \
        --memory 100000 \
        --out $OUTPUT
    fi
done


declare -a Regions=(
"chr1:150600001-155100000"
"chr12:54500001-56200000"
"chr17:33500001-39800000"
"chr2:143400001-147900000"
"chr2:102100001-105300000"
"chr2:236400001-242193529"
"chr6:30500001-32100000"
"chr10:6600001-12200000"
"chr11:75500001-77400000"
"chr12:46000001-48700000"
"chr15:59000001-63400000"
"chr19:31900001-35100000"
)

