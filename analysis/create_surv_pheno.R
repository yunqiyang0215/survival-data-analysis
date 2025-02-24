# Create survival data: all asthma, childhood onset, adulthood onset asthma

# All asthma
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")
pheno = data.frame(cbind(dat$f.eid, dat$status))
pheno$time = dat$AgeAsthma
# Replace NA in age asthma with last age of visit
indx.control = which(is.na(pheno$time))
pheno$time[indx.control] = dat$AgeVisit[indx.control]
colnames(pheno) = c("IID", "event", "time")
saveRDS(pheno, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_all_asthma.rds")

# coa
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")
id = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id_coa.txt")

# Create survival data for childhood onset, defined as onset age <= 12.
# Censor age is at 12. Other people who get asthma later are all treated as censored. 
dat_coa = dat[dat$f.eid %in% id$V1, ]
indx.case <- which(dat_coa$status == 1 & dat_coa$AgeAsthma <= 12)
pheno = data.frame(cbind(dat_coa$f.eid, rep(0, nrow(dat_coa)), dat_coa$AgeAsthma))
colnames(pheno) = c("IID", "event", "time")
pheno$event[indx.case] = 1
# people who never got asthma
pheno$time[is.na(pheno$time)] = 12 
pheno$time[which(pheno$time > 12)]  = 12
saveRDS(pheno, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_coa.rds")


# aoa
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")
id = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/id_aoa.txt")

# Create survival data for adult onset, defined as onset age between 26-65.
# censor age is at 65
dat_aoa = dat[dat$f.eid %in% id$V1, ]
indx.case <- which(dat_aoa$status == 1 & dat_aoa$AgeAsthma > 25 &
                   dat_aoa$AgeAsthma <= 65)
pheno = data.frame(cbind(dat_aoa$f.eid, rep(0, nrow(dat_aoa)), dat_aoa$AgeAsthma))
colnames(pheno) = c("IID", "event", "time")

pheno$event[indx.case] = 1
# people who never got asthma, set censor time to current age or 65 if older. 
pheno$time[is.na(pheno$time)] = pmin(65, dat_aoa$age[is.na(pheno$time)])
pheno$time = as.numeric(pheno$time)
saveRDS(pheno, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_aoa.rds")







