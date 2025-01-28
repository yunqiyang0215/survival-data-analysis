# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
library(susieR)
# 10p14
gwas <- readRDS(file.path("../output/gwas_surv",
                          "all_gwas_chr10_6600001_12200000.rds"))
res <- readRDS(file.path("../output/result202408/all",
                         "fit.susie.chr10_6600001_12200000.rds"))
fit <- res[[1]]
X <- res[[2]]
class(fit) <- c("susie","list")
cs <- susie_get_cs(fit,X)
pvals <- (-log10(gwas[,"p.value.spa"]))
plot(pvals,col = colors,xlab = "SNP",ylab = "-log10(p-value)",cex = 0.8,
     pch = 20,ylim = c(0,max(pvals,na.rm = TRUE) + 0.1))
for (j in 1:length(cs$cs)) {
  snps <- colnames(X)[cs$cs[[j]]]
  indx <- which(rownames(gwas) %in% snps)
  points(indx,-log10(gwas[indx,"p.value.spa"]),pch = 20)
}
