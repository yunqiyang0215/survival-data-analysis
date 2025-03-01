# Script to extract the asthma fine-mapping results for all the
# regions and traits.
library(susieR)
# library(SPACox)
regions <-
  data.frame(filename = c("chr1_150600001_155100000",
                          "chr2_102100001_105300000",
                          "chr2_143400001_147900000",
                          "chr10_6600001_12200000",
                          "chr11_75500001_77400000", 
                          "chr12_46000001_48700000",
                          "chr15_59000001_63400000",
                          "chr17_33500001_39800000"),
             name = c("1q21.3","2q12.1","2q22.3","10p14","11q13.5",
               "12q13.1","15q22.2","17q12"),
             chr = c(1,2,2,10,11,12,15,17))
traits <- c("COA","AOA","AA")

# Repeat for each region and for each trait.
out <- NULL
n <- nrow(regions)
for (i in 1:n) {
  region <- regions[i,"filename"]
  pvar <- read.table(sprintf("../data/geno_regions/%s.pvar",region),sep = "\t",
                     header = TRUE,stringsAsFactors = FALSE,comment.char = "")
  rownames(pvar) <- pvar$ID
  for (trait in traits) {
    cat("region = ",regions[i,"name"],", trait = ",trait,"\n",sep = "")
    if (trait == "COA") {
      gwas <- data.frame(readRDS(paste0("../data/gwas_surv/coa_gwas_",
                                        region,".rds")))
      res  <- readRDS(paste0("../result202408/coa/fit.susie.",region,".rds"))
    }
    else if (trait == "AOA") {
      gwas <- data.frame(readRDS(paste0("../data/gwas_surv/aoa_gwas_"
                                        ,region,".rds")))
      res  <- readRDS(paste0("../result202408/aoa/fit.susie.",region,".rds"))
    } else {
      gwas <- data.frame(readRDS(paste0("../data/gwas_surv/all_gwas_",
                                        region,".rds")))
      res  <- readRDS(paste0("../result202408/all/fit.susie.",region,".rds"))
    }
    fit <- res[[1]]
    X   <- res[[2]]
    class(fit) <- c("susie","list")
    rm(res)

    # Get the GWAS results.
    ids <- sapply(strsplit(colnames(X),"_"),"[[",1)
    gwas <-
      cbind(gwas,
            data.frame(id = sapply(strsplit(rownames(gwas),"_"),"[[",1),
                       counted = sapply(strsplit(rownames(gwas),"_"),"[[",2),
                       stringsAsFactors = FALSE))
    rownames(gwas) <- gwas$id
    gwas <- gwas[ids,]
    pvar <- pvar[ids,]
    
    # Extract the susie results.
    cs   <- susie_get_cs(fit,X)
    pips <- susie_get_pip(fit)
    if (!is.null(cs$cs)) {
      num_cs <- length(cs$cs)
      for (j in 1:num_cs) {
        snps <- cs$cs[[j]]
        dat <- data.frame(region  = regions[i,"name"],
                          trait   = trait,
                          CS      = paste0("L",j),
                          id      = gwas[snps,"id"],
                          pos     = pvar[snps,"POS"],
                          counted = pvar[snps,"REF"], # Same as gwas$counted.
                          alt     = pvar[snps,"ALT"],
                          AF      = gwas[snps,"MAF"],
                          pip     = pips[snps],
                          p.value.spa = gwas[snps,"p.value.spa"],
                          Stat    = gwas[snps,"Stat"],
                          Var     = gwas[snps,"Var"],
                          z       = gwas[snps,"z"])
        out <- rbind(out,dat)
      }
    }
  }
}

# Save the results to a CSV file.
write.csv(format(out,digits = 6,trim = TRUE),
          "ukbiobank_asthma_coxphsusie.csv",
          quote = FALSE,row.names = FALSE)
