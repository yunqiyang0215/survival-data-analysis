# module load gcc/12.1.0 R/4.3.1
library(survival)
library(susieR)
library(logisticsusie)

# Perform a (Bayesian) univariate association analysis using
# the CoxPH model.
surv_uni_fun <- function (x, y, e, prior_variance,
                          estimate_intercept = 0, ...) {
  v0   <- prior_variance                         
  fit  <- coxph(y ~ x + offset(e))
  out  <- summary(fit)$coefficients
  bhat <- out[1,"coef"]
  s    <- out[1,"se(coef)"]
  z    <- bhat/s
  lbf  <- log(s^2/(v0 + s^2))/2 + z^2/2*v0/(v0 + s^2)
  lbf  <- lbf - z^2/2 + summary(fit)$logtest["test"]/2
  v1   <- 1/(1/v0 + 1/s^2)
  mu1  <- v1*bhat/s^2
  return(list(mu = mu1,var = v1,lbf = lbf,
              prior_variance = mu1^2 + v1,
              intercept = 0))
}

# Load the phenotype and genotype data.
pheno <- readRDS("../data/surv_pheno/surv_coa.rds")
geno <- readRDS("../data/geno_finemap_region202408/chr10_6600001_12200000.rds")
rownames(pheno) <- pheno$IID
rownames(geno) <- geno$IID
geno <- geno[,-(1:2)]
geno <- as.matrix(geno)

# Align the two data sets and prepare the genotype and phenotype data
# for the CoxPH-SuSiE analysis.
ids   <- rownames(pheno)
geno  <- geno[ids,]
pheno <- Surv(pheno$time,pheno$event)

### Fit susie
p <- ncol(geno)
fit_coxph <- ser_from_univariate(surv_uni_fun)
num_cores <- 1
print(num_cores)
t0 <- proc.time()
fit <- ibss_from_ser(geno,pheno,L = 10,tol = 0.001,maxit = 10,
                     ser_function = ser_from_univariate(surv_uni_fun),
                     num_cores = num_cores)
t1 <- proc.time()
print(t1 - t0)
save(list = "fit",file = "coxphsusie_ukbb.RData")

