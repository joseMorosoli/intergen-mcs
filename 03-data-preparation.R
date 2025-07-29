###############################################
# Script: Data preparation for MCS PGS analyses
# Author: [Your Name]
# Date: [Date]
# Description: Merge phenotypic & polygenic score data 
# for MCS trios, residualize for PCs, and standardize 
# variables for analysis.
###############################################

# Load required packages
library(tidyverse)
library(foreign)
library(psych)

# Define file path
file_path <- "C:/Users/Jose Morosoli/Documents/UCL/MCS/"

#-----------------------------
# Load Phenotypic Data
#-----------------------------
CMstructure <- read.spss(paste0(file_path,"GENDAC_PINGAULT_mcs_cm_structure_2021_05_12.sav"), to.data.frame = TRUE)
parentstructure <- read.spss(paste0(file_path,"GENDAC_PINGAULT_mcs_parent_structure_2021_05_12.sav"), to.data.frame = TRUE)
parentCMstructure <- read.spss(paste0(file_path,"GENDAC_PINGAULT_mcs_parent_cm_structure_2021_05_12.sav"), to.data.frame = TRUE)
matches <- read.spss(paste0(file_path,"MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav"), to.data.frame = TRUE)

#-----------------------------
# Load Polygenic Scores
#-----------------------------
EApgs  <- read.table(paste0(file_path,"EA3_no23.csv.gz_pred_auto.txt"), header=TRUE)
CPpgs  <- read.table(paste0(file_path,"CP_all.csv.gz_pred_auto.txt"), header=TRUE)
NCPpgs <- read.table(paste0(file_path,"NonCog_no23.csv.gz_pred_auto.txt"), header=TRUE)

#-----------------------------
# Load Principal Components
#-----------------------------
PS_PC <- read.table(paste0(file_path,"mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_PCs.txt"), sep="\t", header=TRUE)
PSPCID <- PS_PC %>% select(-FID)

#-----------------------------
# Helper: Clean Pingault_IDs
#-----------------------------
clean_id <- function(x) gsub(gsub('PNG_05A35003UD','',x[1]),'', x)

#-----------------------------
# SDQ Scores (Child)
#-----------------------------
SDQ <- CMstructure %>% select(Pingault_ID, BEBDTOT, CEBDTOT, DDDEBDTOT, FEBDTOT)
SDQ$Pingault_ID <- clean_id(SDQ$Pingault_ID)

#-----------------------------
# Merge Function for PGS + Phenotypes
#-----------------------------
merge_pgs <- function(PGS, matches, SDQ, label){
  merged <- merge(matches, PGS, by.x="Pingault_ID", by.y="sample.ID") %>% select(-c(9:10))
  child  <- merged %>% filter(MFC == 'C         ') %>% merge(SDQ, by="Pingault_ID")
  mother <- merged %>% filter(MFC == 'M         ')
  father <- merged %>% filter(MFC == 'F         ')
  
  # Rename columns
  child  <- child %>% rename_with(~paste0(label,"_PGS_C"), "final_pred_auto") %>%
    rename(CM_ID = Pingault_ID)
  mother <- mother %>% rename_with(~paste0(label,"_PGS_M"), "final_pred_auto") %>%
    rename(M_ID = Pingault_ID)
  father <- father %>% rename_with(~paste0(label,"_PGS_F"), "final_pred_auto") %>%
    rename(F_ID = Pingault_ID)
  
  # Merge into trio-level dataframe
  list(child, mother, father) %>% reduce(left_join, by="Pingault_FID")
}

#-----------------------------
# Create Trio-Level PGS Data
#-----------------------------
EA_data  <- merge_pgs(EApgs, matches, SDQ, "EA")
CP_data  <- merge_pgs(CPpgs, matches, SDQ, "CP")
NCP_data <- merge_pgs(NCPpgs, matches, SDQ, "NCP")

# Merge all PGS data
all_data <- list(
  EA_data %>% select(-c(PNUM.x, MFC.x, PNUM.y, MFC.y, PNUM, MFC)),
  CP_data %>% select(Pingault_FID, CP_PGS_C, CP_PGS_M, CP_PGS_F),
  NCP_data %>% select(Pingault_FID, NCP_PGS_C, NCP_PGS_M, NCP_PGS_F)
) %>% reduce(left_join, by="Pingault_FID")

#-----------------------------
# Add Child Principal Components
#-----------------------------
trio_ID <- matches %>% select(Pingault_ID, Pingault_FID, MFC)
trio_ID_PSPC <- merge(trio_ID, PSPCID, by.x="Pingault_ID", by.y="IID")
child_PSPC <- trio_ID_PSPC %>% filter(MFC == 'C         ') %>%
  select(-Pingault_ID, -MFC) %>% rename_with(~paste0(.,"_C"), starts_with("PC")) %>%
  mutate(across(starts_with("PC"), scale))
all_data_PSPC <- merge(all_data, child_PSPC, by="Pingault_FID")

#-----------------------------
# Compute Internalising & Externalising Scores
#-----------------------------
internal_raw <- CMstructure %>% select(Pingault_ID, BEMOTION,BPEER,CEMOTION,CPEER,DDEMOTION,DDPEER,FEMOTION,FPEER)
external_raw <- CMstructure %>% select(Pingault_ID, BCONDUCT,BHYPER,CCONDUCT,CHYPER,DDCONDUCT,DDHYPER,FCONDUCT,FHYPER)
internal_raw$Pingault_ID <- clean_id(internal_raw$Pingault_ID)
external_raw$Pingault_ID <- clean_id(external_raw$Pingault_ID)
internal <- internal_raw %>% mutate(across(-Pingault_ID, ~as.numeric(as.character(.))),
                                    BINT = BEMOTION + BPEER, CINT = CEMOTION + CPEER,
                                    DINT = DDEMOTION + DDPEER, FINT = FEMOTION + FPEER)
external <- external_raw %>% mutate(across(-Pingault_ID, ~as.numeric(as.character(.))),
                                    BEXT = BCONDUCT + BHYPER, CEXT = CCONDUCT + CHYPER,
                                    DEXT = DDCONDUCT + DDHYPER, FEXT = FCONDUCT + FHYPER)
INTEXT <- merge(internal, external, by="Pingault_ID")

#-----------------------------
# Merge with PGS + SDQ
#-----------------------------
myData <- all_data_PSPC %>%
  left_join(INTEXT, by=c("CM_ID"="Pingault_ID")) %>%
  left_join(SDQ, by=c("CM_ID"="Pingault_ID"))

#-----------------------------
# Sex Recoding & Clinical Cut-offs
#-----------------------------
myData <- myData %>%
  mutate(sex = as.numeric(SEX_C)-2,
         BCCT = (BEBDTOT >= 15)*1,
         CCCT = (CEBDTOT >= 15)*1,
         DCCT = (DDDEBDTOT >= 15)*1,
         FCCT = (FEBDTOT >= 15)*1)

#-----------------------------
# Add Age at Interview
#-----------------------------
addVars <- read.spss(paste0(file_path,"MDAC-2020-0016-05A-PINGAULT-v4_mcs_cm_structure_pheno_data_2023-03-08_10-37-46.sav"), to.data.frame=TRUE) %>%
  select(PINGAULT_SID,BCASAG00,CHCAGE00,DAGEDY000,FCMCS6AG) %>%
  rename(CM_ID=PINGAULT_SID,BAGE=BCASAG00,CAGE=CHCAGE00,DAGE=DAGEDY000,FAGE=FCMCS6AG) %>%
  mutate(across(BAGE:FAGE, as.numeric))
myData <- left_join(myData, addVars, by="CM_ID") %>%
  distinct(CM_ID,.keep_all=TRUE) %>%
  mutate(across(c(BAGE,FAGE,CAGE,DAGE), scale, .names="{.col}z"))

#-----------------------------
# Add Parent PCs
#-----------------------------
mother_PSPC <- trio_ID_PSPC %>% filter(MFC == 'M         ') %>%
  select(-Pingault_ID,-MFC) %>%
  rename_with(~paste0(.,"_M"), starts_with("PC")) %>%
  mutate(across(starts_with("PC"), scale))
father_PSPC <- trio_ID_PSPC %>% filter(MFC == 'F         ') %>%
  select(-Pingault_ID,-MFC) %>%
  rename_with(~paste0(.,"_F"), starts_with("PC")) %>%
  mutate(across(starts_with("PC"), scale))
myData <- myData %>%
  left_join(mother_PSPC, by="Pingault_FID") %>%
  left_join(father_PSPC, by="Pingault_FID")

#-----------------------------
# Add Parental Education
#-----------------------------
parentstructure$Pingault_ID <- clean_id(parentstructure$Pingault_ID)
parent_edu <- parentstructure %>% select(Pingault_ID,ADACAQ00,ADDNVQ00)
mother_edu <- mother_PSPC %>% left_join(parent_edu, by=c("Pingault_FID"="Pingault_ID")) %>%
  rename(ADACAQ00_M=ADACAQ00,ADDNVQ00_M=ADDNVQ00)
father_edu <- father_PSPC %>% left_join(parent_edu, by=c("Pingault_FID"="Pingault_ID")) %>%
  rename(ADACAQ00_F=ADACAQ00,ADDNVQ00_F=ADDNVQ00)
myData <- myData %>%
  left_join(mother_edu %>% select(Pingault_FID,ADACAQ00_M,ADDNVQ00_M), by="Pingault_FID") %>%
  left_join(father_edu %>% select(Pingault_FID,ADACAQ00_F,ADDNVQ00_F), by="Pingault_FID")

#-----------------------------
# Residualize PGS & SDQ for PCs
#-----------------------------
resid_var <- function(y, x, df){ residuals(lm(as.formula(paste(y,"~",paste(x,collapse="+"))), data=df)) }
pcs_child <- paste0("PC",1:20,"_C")
pcs_mother <- paste0("PC",1:20,"_M")
pcs_father <- paste0("PC",1:20,"_F")

for(trait in c("EA","NCP","CP")){
  for(role in c("C","M","F")){
    myData[[paste0(trait,"_PGS_",role,"_res")]] <- resid_var(paste0(trait,"_PGS_",role),
                                                             ifelse(role=="C",c(pcs_child,"sex"),
                                                                    ifelse(role=="M",pcs_mother,pcs_father)),
                                                             myData)
  }
}

#-----------------------------
# Standardize Residuals
#-----------------------------
std_vars <- function(...){ scale(c(...)) }
n <- nrow(myData)
for(trait in c("EA","NCP","CP")){
  all <- std_vars(myData[[paste0(trait,"_PGS_C_res")]],
                  myData[[paste0(trait,"_PGS_M_res")]],
                  myData[[paste0(trait,"_PGS_F_res")]])
  myData[[paste0(trait,"_PGS_C_zres")]] <- all[1:n]
  myData[[paste0(trait,"_PGS_M_zres")]] <- all[(n+1):(2*n)]
  myData[[paste0(trait,"_PGS_F_zres")]] <- all[(2*n+1):(3*n)]
}

#-----------------------------
# Save Final Dataset
#-----------------------------
save(myData, file="MCS_PGS_cleaned.RData")
