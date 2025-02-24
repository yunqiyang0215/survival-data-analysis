# Creat asthma age and binary case-control status.
# Merge additional covariates with pheno file. 

dat = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma_age_info.rds")
covar = readRDS("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/data/covar.rds") 
out.file <- "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma.rds"

AgeAsthma <- apply(dat[, c("f.3786.0.0", "f.3786.1.0", "f.3786.2.0", "f.3786.3.0", "f.22147.0.0")], 1, function(x) {
  if (all(is.na(x))) {
    return(NA)  
  } else {
    x1 = x[!is.na(x)]
    if (length(x1[x1 > 0]) == 0) {  
      return(-1)  # this is to remove people coded as -1, -3. (Don't know/Don't wish to answer)
    } else {
      return(min(x1[x1 > 0]))  # Use min of reported age as asthma age onset
    }
  }
})

status = rep(0, length(AgeAsthma))
status[!is.na(AgeAsthma)] = 1
status[which(AgeAsthma == -1)] = -1

dat$status = status
dat$AgeAsthma = AgeAsthma

cat(sprintf("Check number of NAs in age at last visit"))
sum(is.na(dat$AgeVisit))
  
cat(sprintf("Remove %d individuals who don't know/don't wish to answer asthma age onset. \n",sum(status == -1)))
dat2 = dat[which(dat$status != -1), ]

# Merge pheno with covar
dat3 = merge(dat2, covar, by.x = "f.eid", by.y = "id")
saveRDS(dat3, out.file)

# To get censoring age, simply replace NA with age at last visit. 
# dat$AgeAsthma[dat$status == 0] = dat$AgeVisit[dat$status == 0]

