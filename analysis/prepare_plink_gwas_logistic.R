# prepare plink input file for case-control gwas. 
# (all asthma, coa asthma, aoa asthma)

# all asthma
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")

out.id.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id.txt"
out.pheno.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/pheno.txt"
out.covar.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar.txt"

IID = dat$f.eid
FID = dat$f.eid

# prepare id file
id.table = cbind(FID, IID)
write.table(id.table, file= out.id.file, quote=FALSE, row.names = FALSE, col.names = FALSE)

# prepare pheno file
pheno.table = cbind(FID, IID, dat$status)
colnames(pheno.table) = c("FID", "IID", "asthma")
write.table(pheno.table, file= out.pheno.file, quote=FALSE, row.names = FALSE)

# prepare covariates
cols = which(colnames(dat) %in% c("age",paste0("pc_genetic", 1:10)))
Z = scale(apply(dat[, cols], 2, as.numeric))
Z = as.data.frame(Z)
Z$age2 = Z$age^2
Z$sex = factor(dat$sex)
cov.table = cbind(FID, IID, Z)
write.table(cov.table, file= out.covar.file, quote=FALSE, row.names = FALSE)


# coa asthma
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")
out.id.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id_coa.txt"
out.pheno.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/pheno_coa.txt"
out.covar.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_coa.txt"


# coa: onset age <= 12.
# cases: people who got asthma before 12
# controls: people who never got asthma & people who got asthma after age 26
indx_case = which(dat$AgeAsthma <= 12 & dat$status == 1)
indx_control = c(which(is.na(dat$AgeAsthma)), which(dat$AgeAsthma >= 26))

# set all control value to 0, and subset to dat_coa
dat$status[indx_control] = 0
dat = dat[c(indx_case, indx_control), ]

IID = dat$f.eid
FID = dat$f.eid

# prepare id file
id.table = cbind(FID, IID)
write.table(id.table, file= out.id.file, quote=FALSE, row.names = FALSE, col.names = FALSE)

# prepare pheno file
pheno.table = cbind(FID, IID, dat$status)
colnames(pheno.table) = c("FID", "IID", "asthma_coa")
write.table(pheno.table, file= out.pheno.file, quote=FALSE, row.names = FALSE)

# prepare covariates
cols = which(colnames(dat) %in% c("age",paste0("pc_genetic", 1:10)))
Z = scale(apply(dat[, cols], 2, as.numeric))
Z = as.data.frame(Z)
Z$age2 = Z$age^2
Z$sex = factor(dat$sex)
cov.table = cbind(FID, IID, Z)
write.table(cov.table, file= out.covar.file, quote=FALSE, row.names = FALSE)


# aoa asthma
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")
out.id.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id_aoa.txt"
out.pheno.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/pheno_aoa.txt"
out.covar.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_aoa.txt"


# aoa: asthma onset age between (25, 65).
# cases: people who got asthma at age 25-65.
# controls: people who never got asthma
indx_case = which(dat$AgeAsthma >=25 & dat$AgeAsthma <= 65 & dat$status == 1)
indx_control = which(is.na(dat$AgeAsthma))

# set all control value to 0, and subset to dat_aoa
dat$status[indx_control] = 0
dat = dat[c(indx_case, indx_control), ]

IID = dat$f.eid
FID = dat$f.eid

# prepare id file
id.table = cbind(FID, IID)
write.table(id.table, file= out.id.file, quote=FALSE, row.names = FALSE, col.names = FALSE)

# prepare pheno file
pheno.table = cbind(FID, IID, dat$status)
colnames(pheno.table) = c("FID", "IID", "asthma_aoa")
write.table(pheno.table, file= out.pheno.file, quote=FALSE, row.names = FALSE)

# prepare covariates
cols = which(colnames(dat) %in% c("age",paste0("pc_genetic", 1:10)))
Z = scale(apply(dat[, cols], 2, as.numeric))
Z = as.data.frame(Z)
Z$age2 = Z$age^2
Z$sex = factor(dat$sex)
cov.table = cbind(FID, IID, Z)
write.table(cov.table, file= out.covar.file, quote=FALSE, row.names = FALSE)


