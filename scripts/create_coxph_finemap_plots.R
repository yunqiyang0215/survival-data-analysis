# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(ggrepel)
library(cowplot)
library(susieR)
#
# TESTING: 10p14, AA
#
# Load the SPACox association testing results and the fine-mapping
# results, and align the two sets of results.
cs_colors <- c("dodgerblue","limegreen","#33a02c","orange")
pvar <- read.table("../data/geno_regions/chr10_6600001_12200000.pvar.gz",
                   sep = "\t",header = TRUE,stringsAsFactors = FALSE,
                   comment.char = "")
gwas <- readRDS(file.path("../output/gwas_surv",
                          "all_gwas_chr10_6600001_12200000.rds"))
res  <- readRDS(file.path("../output/result202408/all",
                          "fit.susie.chr10_6600001_12200000.rds"))
fit  <- res[[1]]
X    <- res[[2]]
rownames(pvar) <- pvar$ID
ids  <- sapply(strsplit(rownames(gwas),"_"),"[[",1)
rownames(gwas) <- ids
ids  <- sapply(strsplit(colnames(X),"_"),"[[",1)
pvar <- pvar[ids,]
gwas <- gwas[ids,]

# Set up the plotting data structure.
pdat <- data.frame(id    = ids,
                   label = "",
                   pos   = pvar$POS,
                   pval  = gwas[,"p.value.spa"],
                   CS    = "none",
                   stringsAsFactors = FALSE)
class(fit) <- c("susie","list")
pips <- susie_get_pip(fit)
cs <- susie_get_cs(fit,X)
num_cs <- length(cs$cs)
for (i in 1:num_cs) {
  snps <- cs$cs[[i]]
  j <- which.max(pips[snps])
  j <- snps[j]
  pdat[snps,"CS"] <- paste0("L",i)
  pdat[j,"label"] <- paste0(pdat[j,"id"]," (L",i,")")
}
pdat <- transform(pdat,
                  pos  = pos/1e6,
                  pval = -log10(pval),
                  CS   = factor(CS,c("none","L1","L2","L3","L4")))

# Create the plot.
p <- ggplot(mapping = aes(x = pos,y = pval,color = CS,label = label)) +
  geom_point(data = subset(pdat,CS == "none"),shape=20,size=2,color="black") +
  geom_point(data = subset(pdat,CS != "none"),shape=20,size=5) +
  geom_point(data = subset(pdat,CS != "none"),shape=20,size=2,color="black") +
  geom_text_repel(data = pdat,color = "black",max.overlaps = Inf,size = 3.5,
                  min.segment.length = 0) +
  scale_color_manual(values = cs_colors) +
  scale_y_continuous(breaks = seq(0,100,5)) +
  labs(x = "base-pair position on chromosome (Mb)",y = "-log10 p-value",
       title = "10p14, AA") +
  theme_cowplot(font_size = 12)

# Save the plots in a PDF.
ggsave("plot.pdf",p,height = 3.5,width = 6)
