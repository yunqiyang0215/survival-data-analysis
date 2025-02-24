# Run this script to prepare the phenotype data for the
# perform_plink_assoc.sh script.
#
# module load gcc/12.1.0 R/4.3.1
#
infiles <- c("pheno","pheno_aoa","pheno_coa")
for (infile in infiles) {
  outfile <- paste0(infile,"_std.txt")
  pheno <- read.csv(sprintf("../data/gwas_logistic/%s.txt",infile),
                    sep = " ",header = TRUE,stringsAsFactors = FALSE)
  y <- pheno[,3]
  y <- drop(scale(y,center = TRUE,scale = TRUE))
  y <- round(y,digits = 6)
  pheno[,3] <- y
  write.table(pheno,outfile,sep = " ",quote = FALSE,row.names = FALSE,
              col.names = TRUE)
}
