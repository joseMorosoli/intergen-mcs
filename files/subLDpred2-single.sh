#!/bin/bash -l
#$ -S /bin/bash
#$ -l h_rt=01:00:00
#$ -l mem=8G
#$ -l tmpfs=8G
#$ -pe smp 8
#$ -N subLDpred2-single
#$ -wd /lustre/scratch/scratch/ucju659/MCS-Jun2024/out/
#$ -m be

export R_LIBS=/lustre/home/ucju659/MyRLibs

## Set up job enviroment:
module -f unload compilers mpi gcc-libs
module load r/recommended

set -o errexit

WORKdir=/lustre/scratch/scratch/ucju659/MCS-Jun2024/R

cd $WORKdir

sumstat=$(awk -F ',' 'NR==1 {print $1}' ${WORKdir}/sumstats_list.csv)                                       #take first column of excel input file
type=$(awk -F ',' 'NR==1 {print $2}' ${WORKdir}/sumstats_list.csv)                                         #take second column of excel input file
ldref=/lustre/scratch/scratch/ucju659/MCS-Jun2024/misc/hapmap3plus/
genofile=/lustre/scratch/scratch/ucju659/MCS-Jun2024/out/mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_trios-hapmap-only.rds
outfile=/lustre/scratch/scratch/ucju659/MCS-Jun2024/out/
sumstatdir=/lustre/scratch/scratch/ucju659/MCS-Jun2024/sumstats/

Rscript --vanilla /lustre/scratch/scratch/ucju659/MCS-Jun2024/R/ldpred2_auto_inf_qc.R -s $sumstat -t $type -m $ldref -g $genofile -o $outfile -c 8 -d $sumstatdir

