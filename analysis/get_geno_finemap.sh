#!/bin/bash

#SBATCH --job-name=get_finemap_region
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --mem=200gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/finemap_region.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/finemap_region.err

module load gcc/12.1.0
module load R/4.2.1


# Define an array of regions
regions=(
    "chr10_6600001_12200000"
    "chr11_75500001_77400000"
    "chr12_46000001_48700000"
    "chr12_54500001_56200000"
    "chr15_59000001_63400000"
    "chr17_33500001_39800000"
    "chr19_31900001_35100000"
    "chr1_150600001_155100000"
    "chr2_102100001_105300000"
    "chr2_143400001_147900000"
    "chr2_236400001_242193529"
    "chr6_30500001_32100000"
)



# Loop through each region and run the R script
for region in "${regions[@]}"; do
    Rscript get_geno_finemap.R  "$region"
done



regions=(
    "chr10_6600001_12200000"
    "chr11_75500001_77400000"
    "chr12_46000001_48700000"
    "chr12_54500001_56200000"
    "chr15_59000001_63400000"
    "chr17_33500001_39800000"
    "chr19_31900001_35100000"
    "chr1_150600001_155100000"
    "chr2_102100001_105300000"
    "chr2_143400001_147900000"
    "chr2_236400001_242193529"
    "chr6_30500001_32100000"
)


