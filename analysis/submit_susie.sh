#!/bin/bash

#SBATCH --job-name=submit_jobs_for_regions
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=500Mb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/submit_jobs.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/submit_jobs.err


regions=(
    "chr6_30500001_32100000"
)


# Loop through each region and submit the job
for region in "${regions[@]}"; do
    sbatch susie.sh "$region"
done

regions=(
    "chr10_6600001_12200000"
    "chr11_75500001_77400000"
    "chr12_46000001_48700000"
    "chr15_59000001_63400000"
    "chr17_33500001_39800000"
    "chr1_150600001_155100000"
    "chr2_102100001_105300000"
    "chr2_143400001_147900000"
    "chr6_30500001_32100000"
)

