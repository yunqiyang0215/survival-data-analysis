#!/bin/bash

#SBATCH --job-name=convert_geno
#SBATCH --partition=tier2q
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb


data_dir="/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions"
PLINK2_EXEC=/gpfs/data/xhe-lab/uk-biobank/tools/plink2
# Change to the directory containing the files
cd "$data_dir"

# Loop through each .pgen file and convert it to .raw format
for pgen_file in *.pgen; do

    # Extract the base name (without extension) for output file naming
    base_name=$(basename "$pgen_file" .pgen)
    # Construct the output file name
    output_file="${base_name}"
    
    # Run PLINK2 command
    $PLINK2_EXEC --pfile "$base_name" --export A --out "$output_file"    
done

