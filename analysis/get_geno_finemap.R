# Description: cut off the region into smaller regions and save the raw genotype for finemapping 
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
region <- args[1]  # e.g., 'chr11_1113000_1750000'
geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/%s.raw", region)
geno <- fread(geno.path)
snp_list <- read.table(sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma_finemap_snplist/sub_%s.txt", region))
# filter snps
geno = as.data.frame(geno)
geno.sub  = geno[, which(colnames(geno) %in% snp_list$x)]
geno = cbind(geno$FID, geno$IID,  geno.sub)
colnames(geno)[1:2] = c("FID", "IID")

saveRDS(geno, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_finemap_regions/", region, ".rds"))



args <- commandArgs(trailingOnly = TRUE)
region <- args[1]  # e.g., 'chr11_1113000_1750000'
geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/%s.raw", region)
geno <- fread(geno.path)
snp_list <- read.table(sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma_finemap_snplist202408/sub_%s.txt", region))
# filter snps
geno = as.data.frame(geno)
geno.sub  = geno[, which(colnames(geno) %in% snp_list$x)]
geno = cbind(geno$FID, geno$IID,  geno.sub)
colnames(geno)[1:2] = c("FID", "IID")

saveRDS(geno, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_finemap_region202408/", region, ".rds"))


