#!/bin/bash
#SBATCH --job-name=SPAcox
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mem=1gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/gwas.cond.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/gwas.cond.err


regions=("chr11_75500001_77400000" "chr17_33500001_39800000" "chr2_102100001_105300000" "chr6_30500001_32100000")

for region in "${regions[@]}"; do
  if [ "$region" == "chr11_75500001_77400000" ]; then
    snps=("rs55646091_A" "rs11236797_A")
  elif [ "$region" == "chr17_33500001_39800000" ]; then
    snps=("rs4795400_T")
  elif [ "$region" == "chr2_102100001_105300000" ]; then
    snps=("rs72823641_A")
  elif [ "$region" == "chr6_30500001_32100000" ]; then
    snps=("rs2428494_A")
  fi

  for snp in "${snps[@]}"; do
    sbatch <<EOF
#!/bin/bash

#SBATCH --job-name=SPAcox_${region}_${snp}
#SBATCH --partition=tier2q
#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --mem=200gb
#SBATCH --output=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/gwas_${region}_${snp}.out
#SBATCH --error=/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/analysis/gwas_${region}_${snp}.err

module load gcc/12.1.0
module load R/4.2.1

Rscript gwas_survival_conditional.R  "$region" "$snp"

EOF

  done
done

