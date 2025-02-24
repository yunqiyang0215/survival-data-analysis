# Run three sets of survival gwas: all asthma, coa and aoa 
library("SPACox", lib="~/R/library")
library(data.table)
library(survival)

# all asthma
pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_all_asthma.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar.txt",
                   header = TRUE)


args <- commandArgs(trailingOnly = TRUE)
region <- args[1]  # e.g., 'chr19_31900001_35100000'
geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/%s.raw", region)
geno <- fread(geno.path)
geno = geno[order(geno$IID), ]

# Make sure data rows align with each other.
all(pheno$IID == covar$IID & covar$IID == geno$IID)

Phen.mtx = cbind(pheno, covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                  "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                  "pc_genetic8", "pc_genetic9", "pc_genetic10")])

snp_columns <- grep("^rs", colnames(geno), value = TRUE)
Geno.mtx <- as.matrix(geno[, ..snp_columns])
rownames(Geno.mtx) = geno$IID


obj.null = SPACox_Null_Model(Surv(time,event)~ .-IID, data=Phen.mtx,
                             pIDs=as.character(Phen.mtx$IID), gIDs=rownames(Geno.mtx))
SPACox.res = SPACox(obj.null, Geno.mtx, min.maf = 0.01)

saveRDS(SPACox.res, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_surv/all_"
                           ,"gwas_", region, ".rds"))

rm(pheno, covar, Geno.mtx, Phen.mtx)

# coa
pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_coa.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_coa.txt",
                   header = TRUE)

# sort data
pheno = pheno[order(pheno$IID), ]
covar = covar[order(covar$IID),	]
geno = geno[order(geno$IID), ]

indx = which(geno$IID %in% pheno$IID)
all(pheno$IID == covar$IID & covar$IID == geno$IID[indx])

# Create geno matrix
Geno.mtx <- as.matrix(geno[indx, ..snp_columns])
rownames(Geno.mtx) = geno$IID[indx]
Phen.mtx = cbind(pheno, covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                  "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                  "pc_genetic8", "pc_genetic9", "pc_genetic10")])


obj.null = SPACox_Null_Model(Surv(time,event)~ .-IID, data=Phen.mtx,
                             pIDs=as.character(Phen.mtx$IID), gIDs=rownames(Geno.mtx))
SPACox.res = SPACox(obj.null, Geno.mtx, min.maf = 0.01)
saveRDS(SPACox.res, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_surv/coa_"
                           ,"gwas_", region, ".rds"))
rm(pheno, covar, Geno.mtx, Phen.mtx)


# aoa
pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_aoa.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_aoa.txt",
                   header = TRUE)

# sort data         
pheno = pheno[order(pheno$IID), ]
covar = covar[order(covar$IID), ]
geno = geno[order(geno$IID), ]

indx = which(geno$IID %in% pheno$IID)
# check data
all(pheno$IID == covar$IID & covar$IID == geno$IID[indx])


# Create geno matrix
Geno.mtx <- as.matrix(geno[indx, ..snp_columns])
rownames(Geno.mtx) = geno$IID[indx]

Phen.mtx = cbind(pheno, covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                  "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                  "pc_genetic8", "pc_genetic9", "pc_genetic10")])


obj.null = SPACox_Null_Model(Surv(time,event)~ .-IID, data=Phen.mtx,
                             pIDs=as.character(Phen.mtx$IID), gIDs=rownames(Geno.mtx))
SPACox.res = SPACox(obj.null, Geno.mtx, min.maf = 0.01)

saveRDS(SPACox.res, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_surv/aoa_"
                           ,"gwas_", region, ".rds"))
rm(pheno, covar, Geno.mtx, Phen.mtx)

