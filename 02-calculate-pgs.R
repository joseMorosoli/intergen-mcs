#########################################
# Calculate polygenic scores
# Author: Jose J. Morosoli
# Date: 23-01-2025
#########################################

# --- Setup (for reproducibility) ---
knitr::opts_chunk$set(echo = TRUE)

#########################################
# Index:
# 1. Preliminary steps
# 2. Using LDpred2 (pipeline)
# Based on: https://github.com/AndreAllegrini/LDpred2
#########################################


#########################################
# Preliminary steps
#########################################

# --- Extract HapMap3 variants from dataset containing only trios ---
# bash:
# plink --bfile ${IN}mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios \
#       --extract /lustre/scratch/scratch/ucju659/hm3plus-pos.list \
#       --make-bed \
#       --out ${OUT}mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios-hapmap-only


# --- Impute missing values with bigsnpr ---
library(bigstatsr)
library(foreach)
library(bigsnpr)

geno_path = '/lustre/scratch/scratch/ucju659/MCS-Jun2024/'
out_path = '/lustre/scratch/scratch/ucju659/MCS-Jun2024/out/'

cat('Reading .bed...')
(rds <- snp_readBed2(
  bedfile = paste0(geno_path,'mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios-hapmap-only.bed'),
  backingfile = paste0(out_path,"mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios-hapmap-only"), 
  ncores = 8))
obj.bigsnp <- snp_attach(rds)
str(obj.bigsnp, max.level = 2)

# Impute missing values
G <- obj.bigsnp$genotypes
counts <- big_counts(G)
counts[, 1:8]

cat('Imputing SNPs...')
G2 <- snp_fastImputeSimple(G, method = "mean2", ncores = 8)

big_counts(G2, ind.col = 1:8)
big_counts(G, ind.col = 1:8)

obj.bigsnp$genotypes <- G2
snp_save(obj.bigsnp)


#########################################
# Run LDpred2 pipeline
#########################################

# bash:
# sumstat=$(awk -F ',' 'NR==1 {print $1}' ${WORKdir}/sumstats_list.csv)
# type=$(awk -F ',' 'NR==1 {print $2}' ${WORKdir}/sumstats_list.csv)                                         
# ldref=/lustre/scratch/scratch/ucju659/MCS-Jun2024/misc/hapmap3plus/
# genofile=/lustre/scratch/scratch/ucju659/MCS-Jun2024/out/mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios-hapmap-only.rds
# outfile=/lustre/scratch/scratch/ucju659/MCS-Jun2024/out/
# sumstatdir=/lustre/scratch/scratch/ucju659/MCS-Jun2024/sumstats/
# Rscript --vanilla /lustre/scratch/scratch/ucju659/MCS-Jun2024/R/ldpred2_auto_inf_qc.R \
#         -s $sumstat -t $type -m $ldref -g $genofile -o $outfile -c 8 -d $sumstatdir


#########################################
# Inside the LDpred2 pipeline
#########################################

library(optparse)
library(ggplot2)
library(bigsnpr)
library(bigreadr)
library(data.table)
library(RhpcBLASctl)
library(future)

# Define options
option_list = list(
  make_option(c("-s", "--sumstats"), type="character", default=NULL, 
              help="Name of GWAS summary statistics. Required columns:
                    case/control: CHR BP A2 A1 NCAS NCON BETA SE
                    continuous: CHR BP A2 A1 N BETA SE", metavar="character"),
  make_option(c("-g", "--geno"), type="character", 
              default="/cluster/projects/p471/people/andrea/LDpred2/geno_data/genoHapMap3plus_N200k.rds", 
              help="Path to genetic dataset in .rds format.", metavar="character"),
  make_option(c("-t", "--type"), type="logical", 
              default = TRUE,
              help="Whether GWAS trait is case/control.", metavar="logical"),
  make_option(c("-o", "--out"), type="character", 
              default="/cluster/projects/p471/people/andrea/LDpred2/out/", 
              help="Output directory.", metavar="character"),
  make_option(c("-d", "--Sdir"), type="character", 
              default="/cluster/p/p471/cluster/people/andrea/LDpred2/sumstats/", 
              help="Sumstats directory.", metavar="character"),
  make_option(c("-m", "--misc"), type="character", 
              default= "/cluster/projects/p471/people/andrea/LDpred2/misc/hapmap3plus/", 
              help="Directory including LD reference info and matrices.", metavar="character"),
  make_option(c("-c", "--cores"), type="integer", 
              default= 32, 
              help="Number of cores.", metavar="number")
) 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Print selected options
print_help(opt_parser)
print(opt)

out_path = opt$out
misc_path = opt$misc
NCORES = opt$cores
geno = opt$geno
sumstatDir = opt$Sdir

start <- Sys.time()
file_log <- paste0(opt$s,".log") 
file.create(file_log)

cat("Running LDpred2 pipeline for: ", opt$s, "\n", file=file_log, append=TRUE)

# --- Load and QC summary statistics ---
sumstats <- bigreadr::fread2(input = paste0(sumstatDir,opt$s))
cat("Loaded ", dim(sumstats)[1], " rows and ", dim(sumstats)[2], " columns.\n", file=file_log, append=TRUE)

# Basic QC: adjust MAF, filter INFO, rename columns, compute n_eff
# (Same as your original pipeline: see bigsnpr documentation)
# ...

# --- Match with LD reference and test data, remove duplicates ---
map_ldref <- readRDS(paste0(misc_path,"map_hm3_plus.rds"))
(info_snp <- tibble::as_tibble(snp_match(sumstats, map_ldref, return_flip_and_rev = T)))
# ...

# --- Compute LD Score regression (heritability) ---
(ldsc <- with(info_snp, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                 chi2 = (beta / beta_se)^2,
                                 sample_size = n_eff,
                                 ncores = NCORES)))
h2_est <- ldsc[["h2"]]
cat("h2 estimate: ", h2_est, "\n", file=file_log, append=TRUE)

# --- Prepare sparse LD matrix ---
for (chr in 1:22) {
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
  corr_chr <- readRDS(paste0(misc_path,"LDref/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, paste0(out_path,"corr_chr"), compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

# --- LDpred2-inf (infinitesimal model) ---
beta_inf <- snp_ldpred2_inf(corr, info_snp, h2 = h2_est)
# Save betas and predicted scores
write.table(cbind(info_snp[,1:4], beta_inf),
            paste0(out_path,opt$sumstats,"_beta_inf.txt"), col.names=T,row.names=F,quote=F)

# --- LDpred2-auto ---
multi_auto <- snp_ldpred2_auto(corr, info_snp, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               allow_jump_sign = FALSE, shrink_corr = 0.95,
                               ncores = NCORES)

saveRDS(multi_auto, file = paste0(out_path,opt$sumstats,"_multi_beta_auto.rds"))
# Filter chains and average betas
(range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))
keep <- (range > (0.95 * quantile(range, 0.95)))
final_beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
write.table(cbind(info_snp[,1:4], final_beta_auto),
            paste0(out_path,opt$sumstats,"_final_beta_auto.txt"), col.names=T,row.names=F,quote=F)

# --- Compute PRS for individuals ---
obj.bigsnp <- snp_attach(geno)
G <- obj.bigsnp$genotypes
pred_auto <- big_prodVec(G, final_beta_auto, ind.col = info_snp[["_NUM_ID_"]], ncores = NCORES)
final_pred_auto <- cbind(obj.bigsnp$fam, pred_auto)
write.table(final_pred_auto, paste0(out_path,opt$sumstats,"_pred_auto.txt"),
            col.names=T, row.names=F, quote=F)

end <- Sys.time()
cat("Analysis complete. Duration: ", round(as.numeric(difftime(end, start , units="mins")),1)," minutes.\n", file=file_log, append=TRUE)
