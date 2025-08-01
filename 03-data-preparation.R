#########################################
#   Polygenic Scores & Phenotypic Data  #
#   Millennium Cohort Study Processing  #
#   Author: Jose J Morosoli
#   Date: 31-07-2025
#########################################

# --- Load required packages ---
library(foreign)
library(tidyverse)
library(dplyr)
library(psych)

# --- Define file paths ---
file_path <- "C:/Users/Jose Morosoli/Documents/UCL/MCS/"

# --- Load SPSS data (phenotypic information) ---
CMstructure       <- read.spss(paste0(file_path, "GENDAC_PINGAULT_mcs_cm_structure_2021_05_12.sav"), to.data.frame = TRUE)
parentstructure   <- read.spss(paste0(file_path, "GENDAC_PINGAULT_mcs_parent_structure_2021_05_12.sav"), to.data.frame = TRUE)
parentCMstructure <- read.spss(paste0(file_path, "GENDAC_PINGAULT_mcs_parent_cm_structure_2021_05_12.sav"), to.data.frame = TRUE)

# --- Load Polygenic Score (PGS) data (generated via LDpred2) ---
EApgs   <- read.table(paste0(file_path, "EA3_no23.csv.gz_pred_auto.txt"), header = TRUE)
CPpgs   <- read.table(paste0(file_path, "CP_all.csv.gz_pred_auto.txt"), header = TRUE)
NCPpgs  <- read.table(paste0(file_path, "NonCog_no23.csv.gz_pred_auto.txt"), header = TRUE)

# --- Load matching IDs for trios ---
matches <- read.spss(paste0(file_path, "MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav"), to.data.frame = TRUE)

# --- Prepare Strengths and Difficulties Questionnaire (SDQ) data ---
SDQ <- subset(CMstructure, select = c(Pingault_ID, BEBDTOT, CEBDTOT, DDDEBDTOT, FEBDTOT))
SDQ$Pingault_ID <- gsub(gsub('PNG_05A35003UD', '', SDQ$Pingault_ID[1]), '', SDQ$Pingault_ID)  # Remove erroneous spaces

# ============================================================
#          CREATE TRIO DATASETS (Child, Mother, Father)
# ============================================================

# Helper function: Create child, mother, and father dataframes for each PGS
create_trio_data <- function(pgs_data, matches, SDQ, pgs_label) {
  matches_pgs <- merge(matches, pgs_data, by.x = "Pingault_ID", by.y = "sample.ID")
  matches_pgs[9:10] <- NULL  # Remove extra columns
  
  # --- Child data ---
  child_pgs <- subset(matches_pgs, MFC == 'C         ')
  child_pgs_sdq <- merge(SDQ, child_pgs, by = "Pingault_ID")
  child_data <- child_pgs_sdq %>%
    rename(
      CM_ID = Pingault_ID,
      SDQ2 = BEBDTOT, SDQ3 = CEBDTOT, SDQ4 = DDDEBDTOT, SDQ6 = FEBDTOT,
      !!paste0(pgs_label, "_C") := final_pred_auto,
      SEX_C = SEX, REGION_C = REGION, ETHNICITY_C = ETHNICITY
    ) %>%
    select(-family.ID)
  
  # --- Mother data ---
  mother_pgs <- subset(matches_pgs, MFC == 'M         ')
  mother_data <- mother_pgs %>%
    rename(
      M_ID = Pingault_ID,
      !!paste0(pgs_label, "_M") := final_pred_auto,
      SEX_M = SEX, REGION_M = REGION, ETHNICITY_M = ETHNICITY
    ) %>%
    select(-family.ID)
  
  # --- Father data ---
  father_pgs <- subset(matches_pgs, MFC == 'F         ')
  father_data <- father_pgs %>%
    rename(
      F_ID = Pingault_ID,
      !!paste0(pgs_label, "_F") := final_pred_auto,
      SEX_F = SEX, REGION_F = REGION, ETHNICITY_F = ETHNICITY
    ) %>%
    select(-family.ID)
  
  # --- Merge all ---
  trio_data <- list(child_data, mother_data, father_data) %>%
    reduce(left_join, by = "Pingault_FID")
  
  return(trio_data)
}

# Create datasets for each PGS
EA_data  <- create_trio_data(EApgs, matches, SDQ, "EA_PGS")
CP_data  <- create_trio_data(CPpgs, matches, SDQ, "CP_PGS")
NCP_data <- create_trio_data(NCPpgs, matches, SDQ, "NCP_PGS")

# ============================================================
#                 MERGE PGS DATASETS
# ============================================================
EA_data_subset  <- subset(EA_data,  select = -c(PNUM.x, MFC.x, PNUM.y, MFC.y, PNUM, MFC))
CP_data_subset  <- subset(CP_data,  select = c(Pingault_FID, CP_PGS_C, CP_PGS_M, CP_PGS_F))
NCP_data_subset <- subset(NCP_data, select = c(Pingault_FID, NCP_PGS_C, NCP_PGS_M, NCP_PGS_F))
all_data        <- list(EA_data_subset, CP_data_subset, NCP_data_subset) %>% reduce(left_join, by = "Pingault_FID")

# ============================================================
#     ADD POPULATION STRATIFICATION PRINCIPAL COMPONENTS
# ============================================================
PS_PC     <- read.table(paste0(file_path, "mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_PCs.txt"), sep = "\t", header = TRUE)
PSPCID    <- subset(PS_PC, select = -FID)
trio_ID   <- subset(matches, select = c(Pingault_ID, Pingault_FID, MFC))
trio_PSPC <- merge(trio_ID, PSPCID, by.x = "Pingault_ID", by.y = "IID")

# Child PCs (scaled)
child_PSPC <- subset(trio_PSPC, MFC == 'C         ', select = -c(Pingault_ID, MFC)) %>%
  rename_with(~ paste0(., "_C"), starts_with("PC")) %>%
  mutate_if(is.numeric, scale)

all_data_PSPC <- list(all_data, child_PSPC) %>% reduce(left_join, by = "Pingault_FID")

# ============================================================
#       EXTRACT & SUMMARIZE INTERNALISING / EXTERNALISING
# ============================================================
internal_raw <- subset(CMstructure, select = c(Pingault_ID, BEMOTION, BPEER, CEMOTION, CPEER, DDEMOTION, DDPEER, FEMOTION, FPEER))
external_raw <- subset(CMstructure, select = c(Pingault_ID, BCONDUCT, BHYPER, CCONDUCT, CHYPER, DDCONDUCT, DDHYPER, FCONDUCT, FHYPER))

# Clean IDs
internal_raw$Pingault_ID <- gsub(gsub('PNG_05A35003UD', '', internal_raw$Pingault_ID[1]), '', internal_raw$Pingault_ID)
external_raw$Pingault_ID <- gsub(gsub('PNG_05A35003UD', '', external_raw$Pingault_ID[1]), '', external_raw$Pingault_ID)

# Convert to numeric
internal <- as.data.frame(lapply(internal_raw[-1], function(x) as.numeric(as.character(x))))
external <- as.data.frame(lapply(external_raw[-1], function(x) as.numeric(as.character(x))))
internal$Pingault_ID <- internal_raw$Pingault_ID
external$Pingault_ID <- external_raw$Pingault_ID

# Compute summed scores
internal <- internal %>%
  mutate(BINT = BEMOTION + BPEER, CINT = CEMOTION + CPEER, DINT = DDEMOTION + DDPEER, FINT = FEMOTION + FPEER)
external <- external %>%
  mutate(BEXT = BCONDUCT + BHYPER, CEXT = CCONDUCT + CHYPER, DEXT = DDCONDUCT + DDHYPER, FEXT = FCONDUCT + FHYPER)

INTEXT <- merge(internal, external, by = "Pingault_ID")

# Merge with main dataset
myData <- merge(all_data_PSPC, INTEXT, by.x = "CM_ID", by.y = "Pingault_ID")
myData <- merge(myData, SDQ, by.x = "CM_ID", by.y = "Pingault_ID")

# --- Additional data processing (ages, parental education, standardization, residualization) ---
# (Section left as in original for brevity — see full script for details)

# --- Save final dataset ---
save(myData, file = 'NATCOMMS_R1.RData')
