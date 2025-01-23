
#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=04:00:00
#$ -l mem=5G
#$ -l tmpfs=15G
#$ -N exclHapMap
#$ -pe smp 10
#$ -wd /lustre/scratch/scratch/ucju659/MCS-Jun2024/


module load plink/1.90b3.40
set -o errexit

IN="/lustre/projects/CMAP/data/MCS/genotype/imputed_QCed/output/"  # path to trios MCS data in CMAP folder
OUT="/lustre/scratch/scratch/ucju659/MCS-Jun2024/"

plink --bfile ${IN}mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios --extract /lustre/scratch/scratch/ucju659/hm3plus-pos.list --make-bed --out ${OUT}mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios-hapmap-only


