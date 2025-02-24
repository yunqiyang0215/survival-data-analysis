library(survival)
library(susieR)
library(tictoc)
library(data.table)
devtools::load_all("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/logisticsusie")


#### Helper functions ######
# Function to calculate log of approximate BF based on Wakefield approximation
# @param z: zscore of the regression coefficient
# @param s: standard deviation of the estimated coefficient
compute_lbf <- function(z, s, prior_variance){
  abf <- sqrt(s^2/(s^2+prior_variance))
  lbf <- log(sqrt(s^2/(s^2+prior_variance))) + z^2/2*(prior_variance/(s^2+prior_variance))
  return(lbf)
}


compute_approx_post_var <- function(z, s, prior_variance){
  post_var <- 1/(1/s^2 + 1/prior_variance)
  return(post_var)
}

# @param post_var: posterior variance
# @param s: standard deviation of the estimated coefficient
# @param bhat: estimated beta effect
compute_approx_post_mean <- function(post_var, s, bhat){
  mu <- post_var/(s^2)*bhat
  return(mu)
}

surv_uni_fun <- function(x, y, o, prior_variance, estimate_intercept = 0, ...){
  fit <- coxph(y~ x + offset(o))
  bhat <- summary(fit)$coefficients[1, 1] # bhat = -alphahat
  sd <- summary(fit)$coefficients[1, 3]
  zscore <- bhat/sd
  lbf <- compute_lbf(zscore, sd, prior_variance)
  lbf.corr <- lbf - bhat^2/sd^2/2+ summary(fit)$logtest[1]/2
  var <- compute_approx_post_var(zscore, sd, prior_variance)
  mu <- compute_approx_post_mean(var, sd, bhat)
  return(list(mu = mu, var=var, lbf=lbf.corr, intercept=0))
}


library(data.table)
library(survival)
dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/gwas/covar.txt",
                   header = TRUE)
# read in genotype
args <- commandArgs(trailingOnly = TRUE)
region <- args[1]  # e.g., 'chr11_1113000_1750000'
geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/geno_finemap/geno/%s.raw", region)
geno <- fread(geno.path)

sum(is.na(geno))

# create survival phenotype
# those who are controls, use last age of visit as censoring time.
indx_control = which(is.na(dat$AgeAsthma))
dat$AgeAsthma[indx_control] = dat$AgeVisit[indx_control]
pheno = data.frame(cbind(dat$f.eid, dat$AgeAsthma, dat$status))
colnames(pheno) = c("IID", "time", "status")

# Order IDs
geno <- geno[order(geno$IID), ]
covar <- covar [order(covar$IID), ]
pheno <- pheno[order(pheno$IID), ]
  
# Subsample controls for susie fit
case_indx = which(pheno$status == 1)
ncase = length(case_indx)
control_indx = sample(which(pheno$status == 0), ncase, replace = FALSE)
sample_indx = c(case_indx, control_indx)
sample_indx = c(1, nrow(geno))

y <- Surv(pheno[sample_indx, ]$time, pheno[sample_indx, ]$status)
snp_columns <- grep("^rs", colnames(geno), value = TRUE)
X <- geno[sample_indx, ..snp_columns]
Z <- covar[sample_indx, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3", "pc_genetic4", "pc_genetic5", 
                                 "pc_genetic6", "pc_genetic7","pc_genetic8", "pc_genetic9", "pc_genetic10")]



### Fit susie
fit_coxph <- ser_from_univariate(surv_uni_fun)
maxiter = 10
L = 10
p = ncol(X)
X = as.matrix(X)
fit.susie <- ibss_from_ser(X, y, Z = Z, L = L, prior_variance = 1., prior_weights = rep(1/p, p), tol = 1e-3, maxit = maxiter, 
                           estimate_intercept = TRUE, ser_function = fit_coxph)

res = list(fit.susie, X)
saveRDS(res, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/res_susie/all/fit.susie.", region, ".rds"))
                           
                           
                           
                           
                           

