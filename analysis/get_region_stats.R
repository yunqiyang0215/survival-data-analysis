# cd /gpfs/data/stephens-lab/uk-biobank-survival/yunqi/survpheno/
# cd ukb_asthma_self_report/data
#
# Set region to one of:
#
#   1 chr1_150600001_155100000
#   2 chr2_102100001_105300000
#   3 chr2_143400001_147900000
#   4 chr6_30500001_32100000
#   5 chr10_6600001_12200000
#   6 chr11_75500001_77400000
#   7 chr12_46000001_48700000
#   8 chr15_59000001_63400000
#   9 chr17_33500001_39800000
#
region <- "chr1_150600001_155100000"
geno <- readRDS(sprintf("geno_finemap_region202408/%s.rds",region))
geno <- geno[,-(1:2)]
geno <- as.matrix(geno)
pvar <- read.table(sprintf("geno_regions/%s.pvar",region),sep = "\t",
                   header = TRUE,stringsAsFactors = FALSE,
                   comment.char = "")
ids  <- colnames(geno)
ids  <- sapply(strsplit(ids,"_"),"[[",1)
rows <- match(ids,pvar$ID)
pvar <- pvar[rows,]
print(diff(range(pvar$POS))/1000)
print(range(pvar$POS))
print(length(ids))
