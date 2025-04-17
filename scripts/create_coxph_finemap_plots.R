# This was the script used to generate the main fine-mapping plots for
# the asthma fine-mapping section of the paper.
library(ggplot2)
library(ggrepel)
library(cowplot)
library(susieR)

# We will create plots summarizing the following fine-mapping
# analyses:
analyses <-
  data.frame(name = c("1q21.3","1q21.3","2q12.1","2q12.1","2q12.1","2q22.3",
                      "10p14","10p14","10p14","11q13.5","11q13.5",
                      "12q13.1","15q22.2","15q22.2","17q12","17q12",
                      "2q12.1","11q13.5"),
             region = c("chr1_150600001_155100000","chr1_150600001_155100000",
                        "chr2_102100001_105300000","chr2_102100001_105300000",
                        "chr2_102100001_105300000","chr2_143400001_147900000",
                        "chr10_6600001_12200000","chr10_6600001_12200000",
                        "chr10_6600001_12200000","chr11_75500001_77400000",
                        "chr11_75500001_77400000","chr12_46000001_48700000",
                        "chr15_59000001_63400000","chr15_59000001_63400000",
                        "chr17_33500001_39800000","chr17_33500001_39800000",
                        "chr2_102100001_105300000","chr11_75500001_77400000"),
             trait = c("COA","all","COA","AOA","all","AOA","COA","AOA",
                       "all","COA","all","AOA","COA","all","COA","all",
                       "all","all"),
             gwas_file = c("coa_gwas_chr1_150600001_155100000.rds",
                           "all_gwas_chr1_150600001_155100000.rds",
                           "coa_gwas_chr2_102100001_105300000.rds",
                           "aoa_gwas_chr2_102100001_105300000.rds",
                           "all_gwas_chr2_102100001_105300000.rds",
                           "aoa_gwas_chr2_143400001_147900000.rds",
                           "coa_gwas_chr10_6600001_12200000.rds",
                           "aoa_gwas_chr10_6600001_12200000.rds",
                           "all_gwas_chr10_6600001_12200000.rds",
                           "coa_gwas_chr11_75500001_77400000.rds",
                           "all_gwas_chr11_75500001_77400000.rds",
                           "aoa_gwas_chr12_46000001_48700000.rds",
                           "coa_gwas_chr15_59000001_63400000.rds",
                           "all_gwas_chr15_59000001_63400000.rds",
                           "coa_gwas_chr17_33500001_39800000.rds",
                           "all_gwas_chr17_33500001_39800000.rds",
                           "all_gwas_chr2_102100001_105300000_rs72823641_A.rds",
                           "all_gwas_chr11_75500001_77400000_rs11236797_A.rds"),
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
  gwas <- readRDS(sprintf("../output/gwas_surv/%s",analyses[analysis,"gwas_file"]))
  res  <- readRDS(sprintf("../result202504/%s/fit.susie.%s.rds",
                          tolower(trait),region))
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
    scale_y_continuous(breaks = seq(0,100,
        ifelse(max(pdat$pval,na.rm = TRUE) > 15,5,2))) +
    labs(x = "base-pair position on chromosome (Mb)",
         y = "-log10 p-value",
         title = paste0(analyses[analysis,"name"],", ",trait)) +
    theme_cowplot(font_size = 12) 
}

# Save the plots in a PDF.
ggsave("coxph_finemap_plots.pdf",
       plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],
                 plots[[6]],plots[[7]],plots[[8]],plots[[9]],plots[[10]],
                 plots[[11]],plots[[12]],plots[[13]],plots[[14]],
                 plots[[15]],plots[[16]],plots[[17]],plots[[18]],
                 nrow = 5,ncol = 4),
       height = 15,width = 16)
