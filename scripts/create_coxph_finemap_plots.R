# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(ggrepel)
library(cowplot)
library(susieR)

# We will create plots summarizing the following fine-mapping
# analyses:
analyses <-
  data.frame(name = c("1q21.3","1q21.3","2q12.1","2q12.1","2q12.1","2q22.3",
                      "2q22.3","10p14","10p14","10p14","11q13.5","11q13.5",
                       "11q13.5","12q13.1","15q22.2","15q22.2","17q12",
                       "17q12"),
             region = c("chr1_150600001_155100000","chr1_150600001_155100000",
                        "chr2_102100001_105300000","chr2_102100001_105300000",
                        "chr2_102100001_105300000","chr2_143400001_147900000",
                        "chr2_143400001_147900000","chr10_6600001_12200000",
                        "chr10_6600001_12200000","chr10_6600001_12200000",
                        "chr11_75500001_77400000","chr11_75500001_77400000",
                        "chr11_75500001_77400000","chr12_46000001_48700000",
                        "chr15_59000001_63400000","chr15_59000001_63400000",
                        "chr17_33500001_39800000","chr17_33500001_39800000"),
             trait = c("COA","all","COA","AOA","all","AOA","all","COA","AOA",
                       "all","COA","AOA","all","AOA","COA","all","COA","all"),
             stringsAsFactors = FALSE)

# Repeat for each of the fine-mapping analyses.
n <- nrow(analyses)
plots <- vector("list",n)
for (analysis in 1:n) {
  region <- analyses[analysis,"region"]
  trait  <- analyses[analysis,"trait"]
  cat("region = ",analyses[analysis,"name"],", trait = ",trait,"\n",sep="")
  
  # Load the SPACox association testing results and the fine-mapping
  # results,and align the two sets of results.
  cs_colors <- c("dodgerblue","limegreen","#33a02c","orange")
  pvar <- read.table(sprintf("../data/geno_regions/%s.pvar.gz",region),
                     sep = "\t",header = TRUE,stringsAsFactors = FALSE,
                     comment.char = "")
  gwas <- readRDS(sprintf("../output/gwas_surv/%s_gwas_%s.rds",
                          tolower(trait),region))
  res  <- readRDS(sprintf("../output/result202408/%s/fit.susie.%d.rds",
                          tolower(trait),region))
  next
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
  plots[[analysis]] <-
    ggplot(mapping = aes(x = pos,y = pval,color = CS,label = label)) +
    geom_point(data=subset(pdat,CS == "none"),shape=20,size=2,color="black") +
    geom_point(data=subset(pdat,CS != "none"),shape=20,size=5) +
    geom_point(data=subset(pdat,CS != "none"),shape=20,size=2,color="black") +
    geom_text_repel(data = pdat,color = "black",max.overlaps = Inf,size = 3.5,
                  min.segment.length = 0) +
    scale_color_manual(values = cs_colors) +
    scale_y_continuous(breaks = seq(0,100,5)) +
    labs(x = "base-pair position on chromosome (Mb)",
         y = "-log10 p-value",
         title = paste0(region,", ",trait)) +
    theme_cowplot(font_size = 12)
}

# Save the plots in two PDFs.
# ggsave("coxph_finemap_plots1.pdf",plots,height = 3.5,width = 6)
# ggsave("coxph_finemap_plots2.pdf",plots,height = 3.5,width = 6)
