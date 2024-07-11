
region = "chr17_33500001_39800000"
gwas.survival = readRDS(paste0("/Users/nicholeyang/downloads/survivalsusie/result/gwas_surv/all_gwas_", region, ".rds"))
gwas.logistic = read.csv(paste0("/Users/nicholeyang/Downloads/survivalsusie/data/gwas_logistic_out/all/gwas_", region,
                                ".asthma.glm.logistic"), sep = "\t")




# merge two datasets
gwas.survival = data.frame(gwas.survival)
snp_ids = sapply(rownames(gwas.survival), function(x) unlist(strsplit(x, "_"))[1])
gwas.survival$ID = snp_ids
dat = merge(gwas.logistic, gwas.survival, by = "ID")
dat2 = dat[, c("ID", "X.CHROM", "POS", "REF", "ALT", "A1", "p.value.spa")]
colnames(dat2) = c("ID", "CHROM", "POS", "REF", "ALT", "A1", "p.value.spa")
write.table(dat2, row.names = FALSE, sep = "\t", quote = FALSE,
            paste0("/Users/nicholeyang/Downloads/survival-data-analysis/data/gwas_", region, ".txt"))


