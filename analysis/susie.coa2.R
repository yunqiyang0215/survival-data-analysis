library(survival)
library(susieR)
library(tictoc)
devtools::load_all("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/logisticsusie.vp")

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

surv_uni_fun <- function(x, y, o, prior_variance, estimate_intercept = 0, num_cores = num_cores, ...){
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



pheno = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/survival-data-analysis/data/surv_pheno/surv_coa.rds")
covar = read.table("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/survival-data-analysis/data/gwas_logistic/covar_coa.txt",
                   header = TRUE)


# read in genotype
args <- commandArgs(trailingOnly = TRUE)
region <- args[1]  # e.g., 'chr11_1113000_1750000'
geno.path <- sprintf("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/survival-data-analysis/data/geno_finemap_region202408/%s.rds", region)
geno <- readRDS(geno.path)


# Order IDs
geno <- geno[order(geno$IID), ]
covar <- covar [order(covar$IID), ]
pheno <- pheno[order(pheno$IID), ]

indx = which(geno$IID %in% pheno$IID)
geno = geno[indx, ]

all.equal(geno$IID, pheno$IID)
all.equal(geno$IID, covar$IID)

y <- Surv(pheno$time, pheno$event)
X <- geno[, 3:ncol(geno)]
Z <- covar[, c("sex", "pc_genetic1", "pc_genetic2", "pc_genetic3", "pc_genetic4", "pc_genetic5",
                          "pc_genetic6", "pc_genetic7","pc_genetic8", "pc_genetic9", "pc_genetic10")]

rm(geno)

### Fit susie
p = ncol(X)
X = as.matrix(X)
fit_coxph <- ser_from_univariate(surv_uni_fun)
niter = 10
L = 10
num_cores = 10
t1 = proc.time()
fit <- ibss_from_ser(X, y,Z=Z, L = L, prior_variance = 1., prior_weights = rep(1/p, p), tol = 1e-3, maxit = niter,
                      estimate_intercept = TRUE, ser_function = fit_coxph, num_cores = num_cores)
t2 = proc.time()
res = list(fit, X, time = t2 - t1)
saveRDS(res, paste0("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/survival-data-analysis/result202504/coa/fit.susie.", region, ".rds"))




