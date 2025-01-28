# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
library(susieR)
region <- regions[i]
res <- readRDS(paste0("/Users/nicholeyang/Downloads/survivalsusie/result/asthma_self_report/result202408/all/fit.susie.",region,".rds"))
gwas <- readRDS(paste0("/Users/nicholeyang/downloads/survivalsusie/result/gwas_surv/all_gwas_",region,".rds"))
fit <- res[[1]]
X <- res[[2]]
class(fit) <- c("susie","list")
cs <- susie_get_cs(fit,X)
colors <- rep("#C3C3C3",nrow(gwas))
pvals <- (-log10(gwas[,"p.value.spa"]))
plot(pvals,col = colors,xlab = "SNP",ylab = "-log10(p-value)",cex = 0.8,pch = 20,main = paste0(region_name[i],": AA"),ylim = c(0,max(pvals,na.rm = TRUE)+0.1))
for (j in 1:length(cs$cs)){
  snps <- colnames(X)[cs$cs[[j]]]
  indx <- which(rownames(gwas) %in% snps)
  colors[indx] <- cols[j]
  points(indx,-log10(gwas[indx,"p.value.spa"]),col = cols[j],pch = 20)
}
