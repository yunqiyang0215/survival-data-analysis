args <- commandArgs(trailingOnly = TRUE)
region <- args[1]
snp <- args[2]

library("SPACox", lib="~/R/library")
library(data.table)
library(survival)

pheno <- readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_all_asthma.rds")
covar <- read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar.txt",
                    header = TRUE)

geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_finemap_regions/%s.rds", region)
geno <- readRDS(geno.path)
geno <- geno[order(geno$IID), ]


Phen.mtx <- cbind(pheno, covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                   "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                   "pc_genetic8", "pc_genetic9", "pc_genetic10")], geno[, snp])

snp_columns <- grep("^rs", colnames(geno), value = TRUE)
Geno.mtx <- as.matrix(geno[, snp_columns])
Geno.mtx <- apply(Geno.mtx, 2, as.numeric)
rownames(Geno.mtx) <- geno$IID

obj.null <- SPACox_Null_Model(Surv(time, event) ~ . - IID, data = Phen.mtx,
                              pIDs = as.character(Phen.mtx$IID), gIDs = rownames(Geno.mtx))

SPACox.res <- SPACox(obj.null, Geno.mtx, min.maf = 0.01)

saveRDS(SPACox.res, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_surv_conditional/all_gwas_", region, "_", snp, ".rds"))

