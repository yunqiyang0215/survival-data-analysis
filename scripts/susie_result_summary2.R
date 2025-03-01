# This script generates CSV file ukbiobank_asthma_cs_summary.csv
# summarizing the CoxPH-SuSiE fine-mapping results on the UK Biobank
# asthma data.
library(susieR)
# library(SPACox)
analyses <-
  data.frame(region = c("chr1_150600001_155100000","chr2_102100001_105300000",
                        "chr2_143400001_147900000","chr10_6600001_12200000",
                        "chr11_75500001_77400000","chr12_46000001_48700000",
                        "chr15_59000001_63400000","chr17_33500001_39800000"),
             name = c("1q21.3","2q12.1","2q22.3","10p14","11q13.5",
               "12q13.1","15q22.2","17q12"),
             trait = c("COA","AA","AOA","AA","AA","AOA","AA","COA"))

# Repeat for each analysis.
out <- NULL
n <- nrow(analyses)
for (i in 1:n) {
  cat(i,"")
  region <- analyses[i,"region"]
  trait <- analyses[i,"trait"]
  if (trait == "COA") {
    gwas <- data.frame(readRDS(paste0("gwas_surv/coa_gwas_",region,".rds")))
    res  <- readRDS(paste0("result202408/coa/fit.susie.",region,".rds"))
  }
  else if (trait == "AOA") {
    gwas <- data.frame(readRDS(paste0("gwas_surv/aoa_gwas_",region,".rds")))
    res  <- readRDS(paste0("result202408/aoa/fit.susie.",region,".rds"))
  } else {
    gwas <- data.frame(readRDS(paste0("../data/gwas_surv/all_gwas_",
                                      region,".rds")))
    res  <- readRDS(paste0("../result202408/all/fit.susie.",region,".rds"))
  }
  fit <- res[[1]]
  X <- res[[2]]
  class(fit) <- c("susie","list")
  rm(res)

  # Extract the key results for the table.
  cs   <- susie_get_cs(fit,X)
  pips <- susie_get_pip(fit)
  ids  <- colnames(X)
  names(pips) <- ids
  if (!is.null(cs$cs)) {
    num_cs <- length(cs$cs)
    for (j in 1:num_cs) {
      cs_size <- length(cs$cs[[j]])
      snps    <- cs$cs[[j]]
      snps    <- ids[snps]
      pips1   <- pips[snps]
      top_snp <- which.max(pips1)
      top_snp <- names(pips1[top_snp])
      out <- rbind(out,
                   data.frame(region = region,
                              trait = trait,
                              CS = j,
                              CS_size = cs_size,
                              purity  = cs$purity[j,"min.abs.corr"],
                              SNP     = top_snp,
                              pval    = gwas[top_snp,"p.value.spa"],
                              PIP     = pips[top_snp],
                              stringsAsFactors = FALSE))
    }
  }
}
cat("\n")

# Save the results to a CSV file.
write.csv(format(out,digits = 6,trim = TRUE),
          "ukbiobank_asthma_cs_summary.csv",
          quote = FALSE,row.names = FALSE)
