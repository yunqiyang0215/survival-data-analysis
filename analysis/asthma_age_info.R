library(data.table)
dat1 = fread("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes/7-jul-2023/ukb676950.tab")
dat2 <- fread("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes/10-apr-2024/ukb678533.tab.gz")
withdrawn<-read.csv("/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/data/withdraw27386_258_20230608.txt",header=F)[,1]

class(dat1) = "data.frame"
class(dat2) = "data.frame"


# Remove withdrawn individuals
dat1 = dat1[!(dat1$f.eid %in% withdrawn), ]

# Step 1: select columns

# Create birth day
colbirth<-which(colnames(dat1) %in% c("f.34.0.0","f.52.0.0"))
birth<-dat1[,colbirth]

# treat birthday as first day of the month
birth<-cbind(birth,1)
birth <- data.frame(lapply(birth, as.numeric))
colnames(birth) = c("year", "month", "day")


# Age when attended assessment centre
colVisit<- grep("21003", colnames(dat1))
dat_Visit<-dat1[,colVisit]

# Find age of most recent visit
AgeVisit <- apply(dat_Visit, 1, function(x) {
  if (all(is.na(x))) {NA  # Return NA if all values are NA
  } else {max(x, na.rm = TRUE)}})

#Age at death
colDeath <- grep("40007", colnames(dat1))
dat_Death_age<-dat1[,colDeath]
AgeDeath <- apply(dat_Death_age, 1, function(x) {
  if (all(is.na(x))) {NA  # Return NA if all values are NA
  } else {min(x, na.rm = TRUE)}})

dat3 <- cbind(dat1[, "f.eid"], birth, AgeVisit, AgeDeath)
colnames(dat3) = c("f.eid", "birth.year", "birth.month", "birth.day", "AgeVisit", "AgeDeath")

# Select self-report asthma onset age
cols <- c("f.eid", "f.3786.0.0", "f.3786.1.0", "f.3786.2.0", "f.3786.3.0", "f.22147.0.0")
dat2 = dat2[, cols]

dat = merge(dat2, dat3, by = "f.eid")
saveRDS(dat, "/gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/ukb_asthma_self_report/data/asthma_age_info.rds")











