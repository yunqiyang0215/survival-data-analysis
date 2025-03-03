# Used to fill out the "top SNP" column for the table summarizing the
# fine-mapping regions considered.
regions <- c("chr1_150600001_155100000",
             "chr10_6600001_12200000",
             "chr11_75500001_77400000",
             "chr15_59000001_63400000",
             "chr17_33500001_39800000",
             "chr2_102100001_105300000",
             "chr12_46000001_48700000",
             "chr2_143400001_147900000")
analysis.prefer <- c("coa","aa","aa","aa",
                     "coa","coa","aa","aoa")

for (i in 1:length(regions)) {
  reg = regions[i]
  if (analysis.prefer[i] == "aa") {
    path.surv = paste0("../data/gwas_surv/all_gwas_",reg,".rds")
    path.logit = paste0("../data/gwas_logistic_out/all/gwas_",reg,
      ".asthma.glm.logistic")
  }
  if (analysis.prefer[i] == "coa") {
    path.surv = paste0("../data/gwas_surv/coa_gwas_",reg,".rds")
    path.logit = paste0("../data/gwas_logistic_out/coa/gwas_",reg,
      ".asthma_coa.glm.logistic")
  }
  if (analysis.prefer[i] == "aoa") {
    path.surv = paste0("../data/gwas_surv/aoa_gwas_",reg,".rds")
    path.logit = paste0("../data/gwas_logistic_out/aoa/gwas_",reg,
      ".asthma_aoa.glm.logistic")
  }
  gwas.surv <- data.frame(readRDS(path.surv))
  gwas.logit <- read.csv(path.logit,header = TRUE,sep = "\t")
  gwas.surv$ID <-
    sapply(1:nrow(gwas.surv),
           function(i) strsplit(rownames(gwas.surv)[i],"_")[[1]][1])
  gwas.surv$variant <- rownames(gwas.surv)
  dat <- merge(gwas.surv,gwas.logit,by = "ID")
  top <- which(dat$p.value.spa == min(dat$p.value.spa,na.rm = TRUE))
  print(reg)
  print(dat[top,])
}
