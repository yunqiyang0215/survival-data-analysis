# Compare top signals of logistic regression and coxph
library(data.table)
library(survival)


# all asthma
pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_all_asthma.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar.txt",
                   header = TRUE)

top_signals = read.csv("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_top_signal/top_signal_AA.csv",
    header = TRUE, row.names = NULL)

loglik.ratio = matrix(NA, ncol = 2, nrow = nrow(top_signals))
colnames(loglik.ratio) = c("survival", "logistic")

for (i in 1:nrow(top_signals)){
    region = top_signals$region[i]
    geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/%s.raw", region)
    geno <- fread(geno.path)
    geno = geno[order(geno$IID), ]
    class(geno) = "data.frame"

    # Make sure data rows align with each other.
    all(pheno$IID == covar$IID & covar$IID == geno$IID)

    snp_columns <- grep("^rs", colnames(geno), value = TRUE)
    geno2 <- geno[, snp_columns]
    colnames(geno2) <- sub("_.*", "", snp_columns)
    dat = cbind(pheno[, c("event", "time")], geno2[, top_signals$snp[i]], covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                  "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                  "pc_genetic8", "pc_genetic9", "pc_genetic10")])

    colnames(dat)[3] = "top_signal"
    mod.full = coxph(Surv(time,event) ~ ., data = dat)
    mod.reduced = coxph(Surv(time,event) ~ . -top_signal, data = dat)
    loglik.ratio[i, 1] <- mod.full$loglik[2] - mod.reduced$loglik[2]

    mod.full <- glm(event ~ .-time, data = dat, family = binomial)
    mod.reduced <- glm(event ~ .-time - top_signal, data = dat, family = binomial)
    loglik.ratio[i, 2] <- logLik(mod.full) - logLik(mod.reduced)
    rm(geno, geno2, dat)  
}

res = cbind(top_signals, loglik.ratio)
write.csv(res, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_top_signal/top_signal_LR_AA.csv")

rm()

# coa
pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_coa.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_coa.txt",
                   header = TRUE)

pheno = pheno[order(pheno$IID), ]
covar = covar[order(covar$IID), ]

top_signals = read.csv("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_top_signal/top_signal_COA.csv",
    header = TRUE, row.names = NULL)

loglik.ratio = matrix(NA, ncol = 2, nrow = nrow(top_signals))
colnames(loglik.ratio) = c("survival", "logistic")

for (i in 1:nrow(top_signals)){
    region = top_signals$region[i]
    geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/%s.raw", region)
    geno <- fread(geno.path)
    geno = geno[order(geno$IID), ]
    class(geno) = "data.frame"

    indx = which(geno$IID %in% pheno$IID)
    geno <- geno[indx, ]

    # Make sure data rows align with each other.
    all(pheno$IID == covar$IID & covar$IID == geno$IID)

    snp_columns <- grep("^rs", colnames(geno), value = TRUE)
    geno2 <- geno[, snp_columns]
    colnames(geno2) <- sub("_.*", "", snp_columns)
    dat = cbind(pheno[, c("event", "time")], geno2[, top_signals$snp[i]], covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                  "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                  "pc_genetic8", "pc_genetic9", "pc_genetic10")])

    colnames(dat)[3] = "top_signal"
    mod.full = coxph(Surv(time,event) ~ ., data = dat)
    mod.reduced = coxph(Surv(time,event) ~ . -top_signal, data = dat)
    loglik.ratio[i, 1] <- mod.full$loglik[2] - mod.reduced$loglik[2]

    mod.full <- glm(event ~ .-time, data = dat, family = binomial)
    mod.reduced <- glm(event ~ .-time - top_signal, data = dat, family = binomial)
    loglik.ratio[i, 2] <- logLik(mod.full)-logLik(mod.reduced)
    rm(geno, geno2, dat)  
}

res = cbind(top_signals, loglik.ratio)
write.csv(res, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_top_signal/top_signal_LR_COA.csv")

rm()




# aoa
pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/surv_pheno/surv_aoa.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_logistic/covar_aoa.txt",
                   header = TRUE)

pheno = pheno[order(pheno$IID), ]
covar = covar[order(covar$IID), ]

top_signals = read.csv("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_top_signal/top_signal_AOA.csv",
    header = TRUE, row.names = NULL)

loglik.ratio = matrix(NA, ncol = 2, nrow = nrow(top_signals))
colnames(loglik.ratio) = c("survival", "logistic")

for (i in 1:nrow(top_signals)){
    region = top_signals$region[i]
    geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/geno_regions/%s.raw", region)
    geno <- fread(geno.path)
    geno = geno[order(geno$IID), ]
    class(geno) = "data.frame"

    indx = which(geno$IID %in% pheno$IID)
    geno <- geno[indx, ]

    # Make sure data rows align with each other.
    all(pheno$IID == covar$IID & covar$IID == geno$IID)

    snp_columns <- grep("^rs", colnames(geno), value = TRUE)
    geno2 <- geno[, snp_columns]
    colnames(geno2) <- sub("_.*", "", snp_columns)
    dat = cbind(pheno[, c("event", "time")], geno2[, top_signals$snp[i]], covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3",
                                  "pc_genetic4", "pc_genetic5", "pc_genetic6", "pc_genetic7",
                                  "pc_genetic8", "pc_genetic9", "pc_genetic10")])

    colnames(dat)[3] = "top_signal"
    mod.full = coxph(Surv(time,event) ~ ., data = dat)
    mod.reduced = coxph(Surv(time,event) ~ . -top_signal, data = dat)
    loglik.ratio[i, 1] <- mod.full$loglik[2] - mod.reduced$loglik[2]

    mod.full <- glm(event ~ .-time, data = dat, family = binomial)
    mod.reduced <- glm(event ~ .-time - top_signal, data = dat, family = binomial)
    loglik.ratio[i, 2] <- logLik(mod.full) - logLik(mod.reduced)
    rm(geno, geno2, dat)  
}

res = cbind(top_signals, loglik.ratio)
write.csv(res, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/gwas_top_signal/top_signal_LR_AOA.csv")
rm()

