#########################################
# Follow-up QC for MCS
# Author: Biyao Wang, Jose J. Morosoli
# Date: 23-01-2025
#########################################

# --- Load libraries ---
library(tidyverse)
library(haven)
library(foreign)

#########################################
# Index:
# 1. Genetic quality control
# 2. Principal component analysis
# 3. Identify and extract complete trios
#########################################

#########################################
# Quality Control
#########################################

# --- Preliminary QC in plink ---
# (Run in terminal)
# DIR="/lustre/scratch/scratch/ucjues9/mcs/genotype/output/"
# module load plink/1.90b3.40
# plink --bfile ${DIR}mcs_imputed_allchr \
#       --mind 0.05  --geno 0.05 --maf 0.01 --hwe 0.000001 \
#       --extract /home/ucjues9/Scratch/mcs/genotype/rsq/70.rsq \
#       --make-bed --out ${DIR}mcs_imputed_allchr_minQC

# --- Follow-up QC (bash checks) ---
# cd /lustre/projects/CMAP/data/MCS/genotype/plink
# wc -l mcs_imputed_allchr.bim         # 38909200 SNPs
# wc -l mcs_imputed_allchr.fam         # 21349 individuals
# wc -l mcs_imputed_allchr_minQC.bim   # 6948515 SNPs
# wc -l mcs_imputed_allchr_minQC.fam   # 21349 individuals

#########################################
# Update FID and IID
#########################################

# bash:
# cd /lustre/projects/CMAP/data/MCS/genotype/imputed_QCed/output
# awk '{print $1}' mcs_imputed_allchr.fam > fam_onlyIID_raw.txt
# cut -d_ -f1,2 fam_onlyIID_raw.txt > fam_onlyIID_new.txt
# paste fam_onlyIID_raw.txt fam_onlyIID_new.txt > fam_onlyIID_rawnew.txt

# --- Map IDs using phenotype file ---
setwd("/lustre/projects/CMAP/data/MCS/phenotype/MDAC-2020-0016-05A-PINGAULT")
df <- read_spss("MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav") 
fam_pheno <- df %>% select(FID_pheno = 2, IID_pheno = 1)

setwd("/lustre/projects/CMAP/data/MCS/genotype/imputed_QCed/output")
fam_geno <- read_tsv("fam_onlyIID_rawnew.txt", col_names = F)
colnames(fam_geno) <- c("IID_geno_raw", "IID_geno_new")
fam_geno$FID_geno_new <- fam_pheno$FID_pheno[match(fam_geno$IID_geno_new, fam_pheno$IID_pheno)]

fam_update <- fam_geno[,c(1,1,3,2)]
colnames(fam_update) <- c("FID_old", "IID_old", "FID_new", "IID_new")
write.table(fam_update, file="fam_update.txt", col.names = T, row.names = F, quote=F, sep="\t")

# bash:
# plink --bfile mcs_imputed_allchr_minQC \
#       --update-ids fam_update.txt \
#       --make-bed \
#       --out mcs_imputed_allchr_minQC_updateID

#########################################
# Pruning
#########################################

# bash:
# plink --bfile ${DIR}/mcs_imputed_allchr \
#       --indep-pairwise 250 50 0.1 \
#       --out ${DIR}/mcs_imputed_allchr_pruned

#########################################
# Heterozygosity QC
#########################################

# bash:
# plink --bfile ${DIR}/mcs_imputed_allchr_minQC_updateID \
#       --extract ${DIR}/mcs_imputed_allchr_pruned.prune.in \
#       --het \
#       --out ${DIR}/mcs_imputed_allchr_minQC_updateID_checkhet

# --- List heterozygosity outliers ---
het <- read.table("mcs_imputed_allchr_minQC_updateID_checkhet.het", head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | 
                    (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)))
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE)
write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)

# bash:
# sed 's/"// g' fail-het-qc.txt | awk '{print $1, $2}' | tail -n +2 > het_fail_list.txt
# plink --bfile mcs_imputed_allchr_minQC_updateID \
#       --remove het_fail_list.txt \
#       --make-bed \
#       --out mcs_imputed_allchr_minQC_updateID_hetQC

#########################################
# Relatedness QC
#########################################

fam_pheno <- read_spss("MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav") 
df <- fam_pheno %>% select(FID = 2, IID = 1, MFC) 
df %>% count(FID, MFC) %>% count(MFC, n, name = "total")

#########################################
# Separate children and parents
#########################################

df <- read_spss("MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav") 
fam_pheno <- df %>% select(IID = 1, MFC) 
fam_geno <- read.table("mcs_imputed_allchr_minQC_updateID.fam", sep = " " , header = F)
colnames(fam_geno) <- c("FID", "IID", "X1", "X2", "X3", "X4")
fam_geno$MFC <- fam_pheno$MFC[match(fam_geno$IID, fam_pheno$IID)]

child  <- fam_geno %>% filter(MFC == "C") 
parent <- fam_geno %>% filter(MFC != "C") 

write.table(child[,-7], file="children.txt", col.names = F, row.names = F, quote=F, sep = "\t")
write.table(parent[,-7], file="parents.txt", col.names = F, row.names = F, quote=F, sep = "\t")

# bash:
# plink --bfile ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC \
#       --keep ${DIR}/children.txt \
#       --make-bed \
#       --out ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC_children
# plink --bfile ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC \
#       --keep ${DIR}/parents.txt \
#       --make-bed \
#       --out ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC_parents

#########################################
# Get unrelated children and parents
#########################################
# bash:
# king -b ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC_children.bed \
#      --unrelated --degree 3 --cpus 36 \
#      --prefix ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC_children_3rd_
# king -b ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC_parents.bed \
#      --unrelated --degree 3 --cpus 36 \
#      --prefix ${DIR}/mcs_imputed_allchr_minQC_updateID_hetQC_parents_3rd_

# bash:
# cat mcs_imputed_allchr_minQC_updateID_hetQC_children_3rd_unrelated.txt mcs_imputed_allchr_minQC_updateID_hetQC_parents_3rd_unrelated.txt > mcs_imputed_allchr_minQC_updateID_hetQC_unrelated.txt
# plink --bfile mcs_imputed_allchr_minQC_updateID_hetQC \
#       --keep mcs_imputed_allchr_minQC_updateID_hetQC_unrelated.txt \
#       --make-bed \
#       --out mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC

#########################################
# Principal Component Analysis
#########################################

# bash:
# plink --bfile mcs_imputed_allchr_minQC_updateID_hetQC \
#       --extract mcs_imputed_allchr_pruned.prune.in \
#       --pca \
#       --out mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC

# --- Import and save PCs ---
PCs = read.table("output/mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC.eigenvec", header=FALSE)
colnames(PCs) = c("FID", "IID", paste0("PC",1:20))
write.table(PCs, file="output/mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_PCs.csv",col.names = T, row.names = F, quote=F, sep = "\t")
save(PCs, file = "output/mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_PCs.RData")

#########################################
# Get Complete Trios
#########################################

fam_pheno <- read_spss("MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav") 
df_pheno <- fam_pheno %>% select(FID = 2, IID = 1, MFC) 
fam_geno <- read.table("mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC.fam", sep = " " , header = F)
colnames(fam_geno) <- c("FID", "IID", "PID","MID","Sex","Pheno")
df_geno <- fam_geno %>% select(FID = 1, IID = 2)
df_geno$MFC <- df_pheno$MFC[match(df_geno$IID, df_pheno$IID)]

children <- subset(df_geno,MFC=="C")
mothers <- subset(df_geno,MFC=="M")
fathers <- subset(df_geno,MFC=="F")
commonFID <- Reduce(intersect, list(children$FID, mothers$FID, fathers$FID))
trios <- fam_geno[fam_geno$FID %in% commonFID,]
write.table(trios, file="mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_completetrios.txt", row=F, col=F, qu=F)

# bash:
# plink --bfile mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC \
#       --keep mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_completetrios.txt \
#       --make-bed \
#       --out mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios
