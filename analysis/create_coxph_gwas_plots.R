# TO DO: Explain here what this script is for, and how to use it.
library(ggplot2)
library(ggrepel)
library(cowplot)
region <- "chr1_150600001_155100000"
aoa  <- readRDS(sprintf("../output/gwas_surv/aoa_gwas_%s.rds",region))
coa  <- readRDS(sprintf("../output/gwas_surv/coa_gwas_%s.rds",region))
res  <- data.frame(id = rownames(aoa),
                   aoa = -log10(aoa[,"p.value.spa"]),
                   coa = -log10(coa[,"p.value.spa"]))
pvar <- read.table(sprintf("../data/geno_regions/%s.pvar.gz",region),
                   sep = "\t",header = TRUE,stringsAsFactors = FALSE,
                   comment.char = "")
ids  <- res$id
ids  <- sapply(strsplit(ids,"_"),"[[",1)
rows <- match(ids,pvar$ID)
pvar <- pvar[rows,]
res  <- cbind(res,pvar["POS"])
rownames(res) <- NULL
pdat <- res
pdat <- transform(pdat,POS = POS/1e6)
rows <- c(which.max(pdat$coa),which.max(pdat$aoa))
pdat[-rows,"id"] <- ""
p <- ggplot(pdat) +
  geom_point(mapping = aes(x = POS,y = coa),color = "orangered",size = 1) +
  geom_point(mapping = aes(x = POS,y = -aoa),color = "dodgerblue",size = 1) + 
  geom_text_repel(mapping = aes(x = POS,y = coa,label = id),
                  color = "black",max.overlaps = Inf,size = 2.5,
                  min.segment.length = 0) +
  geom_text_repel(mapping = aes(x = POS,y = -aoa,label = id),
                  color = "black",max.overlaps = Inf,size = 2.5,
                  min.segment.length = 0) +
  geom_hline(yintercept = 0,color = "black",linetype = "dotted") +
  scale_x_continuous(breaks = seq(0,200,0.5)) +
  scale_y_continuous(breaks = seq(-100,100,5)) +
  labs(x = "position (Mb)",y = "log10(p-value)") +
  theme_cowplot(font_size = 10)
ggsave("coxph_gwas_plot.pdf",p,height = 3,width = 5)
