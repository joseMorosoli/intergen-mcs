#########################################
#   Polygenic Scores & Phenotypic Data  #
#   Millennium Cohort Study Processing  #
#   Author: Jose J Morosoli
#   Revised: 20-11-2026
#########################################
# Reconstructed Nature Communications SEM dataset; validated against historical aframe_svy.

# --- Load required packages ---
library(foreign)
library(tidyverse)

# ============================================================
#          1. INPUT FILES
# ============================================================
# --- Define file paths ---
file_path <- ".../"

CMstructure <- read.spss(file.path(file_path, "GENDAC_PINGAULT_mcs_cm_structure_2021_05_12.sav"), to.data.frame = TRUE)
parentstructure <- read.spss(file.path(file_path, "GENDAC_PINGAULT_mcs_parent_structure_2021_05_12.sav"), to.data.frame = TRUE)
matches <- read.spss(file.path(file_path, "MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav"), to.data.frame = TRUE)

EApgs <- read.table(file.path(file_path, "EA3_no23.csv.gz_pred_auto.txt"), header = TRUE)
CPpgs <- read.table(file.path(file_path, "CP_all.csv.gz_pred_auto.txt"), header = TRUE)
NCPpgs <- read.table(file.path(file_path, "NonCog_no23.csv.gz_pred_auto.txt"), header = TRUE)
PS_PC <- read.table(file.path(file_path, "mcs_imputed_allchr_minQC_updateID_hetQC_relatedQC_PCs.txt"),
                    sep = "\t", header = TRUE)

age_data <- read.spss(file.path(file_path, "MDAC-2020-0016-05A-PINGAULT-v4_mcs_cm_structure_pheno_data_2023-03-08_10-37-46.sav"),
                      to.data.frame = TRUE)

# SPSS identifiers may contain padding.
clean_id <- function(x) gsub("[[:space:]]+", "", as.character(x))
CMstructure$Pingault_ID <- clean_id(CMstructure$Pingault_ID)
parentstructure$Pingault_ID <- clean_id(parentstructure$Pingault_ID)
matches$Pingault_ID <- clean_id(matches$Pingault_ID)
matches$Pingault_FID <- clean_id(matches$Pingault_FID)
matches$MFC <- trimws(as.character(matches$MFC))
for (pgs in c("EApgs", "CPpgs", "NCPpgs")) {
  x <- get(pgs)
  x$sample.ID <- clean_id(x$sample.ID)
  assign(pgs, x)
}
PS_PC$IID <- clean_id(PS_PC$IID)
age_data$PINGAULT_SID <- clean_id(age_data$PINGAULT_SID)

# ============================================================
#          2. CHILD SDQ PHENOTYPES
# ============================================================
SDQ <- CMstructure %>%
  select(Pingault_ID, BEBDTOT, CEBDTOT, DDDEBDTOT, FEBDTOT) %>%
  mutate(across(-Pingault_ID, ~ as.numeric(as.character(.x))))

internal <- CMstructure %>%
  select(Pingault_ID, BEMOTION, BPEER, CEMOTION, CPEER,
         DDEMOTION, DDPEER, FEMOTION, FPEER) %>%
  mutate(across(-Pingault_ID, ~ as.numeric(as.character(.x)))) %>%
  mutate(BINT = BEMOTION + BPEER, CINT = CEMOTION + CPEER,
         DINT = DDEMOTION + DDPEER, FINT = FEMOTION + FPEER)

external <- CMstructure %>%
  select(Pingault_ID, BCONDUCT, BHYPER, CCONDUCT, CHYPER,
         DDCONDUCT, DDHYPER, FCONDUCT, FHYPER) %>%
  mutate(across(-Pingault_ID, ~ as.numeric(as.character(.x)))) %>%
  mutate(BEXT = BCONDUCT + BHYPER, CEXT = CCONDUCT + CHYPER,
         DEXT = DDCONDUCT + DDHYPER, FEXT = FCONDUCT + FHYPER)

INTEXT <- inner_join(internal, external, by = "Pingault_ID")
stopifnot(!anyDuplicated(SDQ$Pingault_ID), !anyDuplicated(INTEXT$Pingault_ID))

# ============================================================
#          3. CREATE TRIO DATASETS (CHILD, MOTHER, FATHER)
# ============================================================
create_trio_data <- function(pgs_data, label) {
  merged <- inner_join(matches, pgs_data, by = c("Pingault_ID" = "sample.ID"))
  
  # --- Child data ---
  child <- merged %>% filter(MFC == "C") %>%
    select(Pingault_FID, CM_ID = Pingault_ID,
           SEX_C = SEX, REGION_C = REGION, ETHNICITY_C = ETHNICITY,
           !!paste0(label, "_PGS_C") := final_pred_auto) %>%
    inner_join(SDQ, by = c("CM_ID" = "Pingault_ID"))
  
  # --- Mother data ---
  mother <- merged %>% filter(MFC == "M") %>%
    select(Pingault_FID, M_ID = Pingault_ID,
           SEX_M = SEX, REGION_M = REGION, ETHNICITY_M = ETHNICITY,
           !!paste0(label, "_PGS_M") := final_pred_auto)
  
  # --- Father data ---
  father <- merged %>% filter(MFC == "F") %>%
    select(Pingault_FID, F_ID = Pingault_ID,
           SEX_F = SEX, REGION_F = REGION, ETHNICITY_F = ETHNICITY,
           !!paste0(label, "_PGS_F") := final_pred_auto)
  
  # A repeated parent per family would multiply rows in the trio merge.
  stopifnot(!anyDuplicated(child$CM_ID),
            !anyDuplicated(mother$Pingault_FID),
            !anyDuplicated(father$Pingault_FID))
  child %>% left_join(mother, by = "Pingault_FID") %>%
    left_join(father, by = "Pingault_FID")
}

EA_data <- create_trio_data(EApgs, "EA")
CP_data <- create_trio_data(CPpgs, "CP")
NCP_data <- create_trio_data(NCPpgs, "NCP")

myData <- EA_data %>%
  left_join(CP_data %>% select(Pingault_FID, CM_ID, starts_with("CP_PGS_")),
            by = c("Pingault_FID", "CM_ID")) %>%
  left_join(NCP_data %>% select(Pingault_FID, CM_ID, starts_with("NCP_PGS_")),
            by = c("Pingault_FID", "CM_ID")) %>%
  inner_join(INTEXT, by = c("CM_ID" = "Pingault_ID"))

# The historical SEX_C factor coding was verified: Male = 0, Female = 1.
# Unlike the historical script, missing SDQ totals remain NA, not non-cases.
myData <- myData %>%
  mutate(sex = as.numeric(SEX_C) - 2,
         BCCT = as.integer(BEBDTOT >= 15), CCCT = as.integer(CEBDTOT >= 15),
         DCCT = as.integer(DDDEBDTOT >= 15), FCCT = as.integer(FEBDTOT >= 15))

# ============================================================
#          4. AGE AT INTERVIEW
# ============================================================
addVars <- age_data %>%
  select(CM_ID = PINGAULT_SID, BAGE = BCASAG00, CAGE = CHCAGE00,
         DAGE = DAGEDY000, FAGE = FCMCS6AG) %>%
  # Match the recovered historical script's factor-to-numeric conversion.
  mutate(across(c(BAGE, CAGE, DAGE, FAGE), as.numeric)) %>%
  filter(!is.na(CM_ID), CM_ID != "")

# The age source also contains records without child IDs; exclude those only.
stopifnot(!anyDuplicated(addVars$CM_ID))
myData <- left_join(myData, addVars, by = "CM_ID") %>%
  mutate(across(c(BAGE, CAGE, DAGE, FAGE),
                ~ as.numeric(scale(.x)), .names = "{.col}z"))

# ============================================================
#          5. POPULATION STRATIFICATION PRINCIPAL COMPONENTS
# ============================================================
trio_ID <- matches %>% select(Pingault_ID, Pingault_FID, MFC)
trio_PSPC <- inner_join(trio_ID, select(PS_PC, -FID),
                        by = c("Pingault_ID" = "IID"))

role_pcs <- function(role, suffix, id_col) {
  out <- trio_PSPC %>%
    filter(MFC == role) %>%
    select(Pingault_FID, Pingault_ID, matches("^PC[0-9]+$")) %>%
    rename_with(~ paste0(.x, "_", suffix), matches("^PC[0-9]+$")) %>%
    mutate(across(matches(paste0("^PC[0-9]+_", suffix, "$")),
                  ~ as.numeric(scale(.x))))
  names(out)[names(out) == "Pingault_ID"] <- id_col
  # Family IDs may repeat; individual IDs must not.
  stopifnot(!anyDuplicated(out[[id_col]]))
  out
}

child_PSPC <- role_pcs("C", "C", "CM_ID")
mother_PSPC <- role_pcs("M", "M", "M_ID")
father_PSPC <- role_pcs("F", "F", "F_ID")

myData <- myData %>%
  left_join(child_PSPC, by = c("Pingault_FID", "CM_ID")) %>%
  left_join(mother_PSPC, by = c("Pingault_FID", "M_ID")) %>%
  left_join(father_PSPC, by = c("Pingault_FID", "F_ID"))

# ============================================================
#          6. PARENTAL EDUCATION
# ============================================================
parent_edu <- parentstructure %>%
  select(Pingault_ID, ADACAQ00, ADDNVQ00)
stopifnot(!anyDuplicated(parent_edu$Pingault_ID))

myData <- myData %>%
  left_join(parent_edu %>%
              dplyr::rename(M_ID = Pingault_ID,
                            ADACAQ00_M = ADACAQ00, ADDNVQ00_M = ADDNVQ00),
            by = "M_ID") %>%
  left_join(parent_edu %>%
              dplyr::rename(F_ID = Pingault_ID,
                            ADACAQ00_F = ADACAQ00, ADDNVQ00_F = ADDNVQ00),
            by = "F_ID")

# ============================================================
#          7. RESIDUALISE POLYGENIC SCORES
# ============================================================
# na.exclude retains the original row alignment.
for (trait in c("EA", "CP", "NCP")) {
  for (role in c("C", "M", "F")) {
    score <- paste0(trait, "_PGS_", role)
    pcs <- paste0("PC", 1:20, "_", role)
    if (role == "C") pcs <- c(pcs, "sex")
    myData[[paste0(score, "_res")]] <- as.numeric(residuals(
      lm(reformulate(pcs, response = score), data = myData, na.action = na.exclude)
    ))
  }
  
  # Historical convention: standardise each PGS across child, mother and father.
  vars <- paste0(trait, "_PGS_", c("C", "M", "F"), "_res")
  z <- as.numeric(scale(unlist(myData[vars], use.names = FALSE)))
  n <- nrow(myData)
  myData[[paste0(trait, "_PGS_C_zres")]] <- z[1:n]
  myData[[paste0(trait, "_PGS_M_zres")]] <- z[(n + 1):(2 * n)]
  myData[[paste0(trait, "_PGS_F_zres")]] <- z[(2 * n + 1):(3 * n)]
}

# ============================================================
#          8. RESIDUALISE INTERNALISING AND EXTERNALISING SCORES
# ============================================================
# Original historical models: outcome ~ 20 child PCs + sex + wave-specific age.
# na.exclude retains row alignment and leaves incomplete observations as NA.
child_pcs <- paste0("PC", 1:20, "_C")
for (wave in c("B", "C", "D", "F")) {
  for (domain in c("INT", "EXT")) {
    outcome <- paste0(wave, domain)
    predictors <- c(child_pcs, "sex", paste0(wave, "AGEz"))
    myData[[paste0(outcome, "_res")]] <- as.numeric(residuals(
      lm(reformulate(predictors, response = outcome),
         data = myData, na.action = na.exclude)
    ))
  }
}

# ============================================================
#          9. HISTORICAL VARIABLE NAMES FOR SEM
# ============================================================
trait_names <- c(EA = "EA_lee", CP = "EA_cog", NCP = "EA_noncog")
for (trait in names(trait_names)) {
  for (role in c("C", "M", "F")) {
    myData[[paste0(trait_names[[trait]], "_", role, "_res")]] <-
      myData[[paste0(trait, "_PGS_", role, "_res")]]
  }
}

# ============================================================
#          10. SURVEY DESIGN WEIGHTS, STRATA AND CLUSTERS
# ============================================================
weights <- read.spss(
  file.path(file_path, "MCS2025",
            "svy_design-PINGAULT_mcs_data_struct_fam_2024-09-18_13-04-10.sav"),
  to.data.frame = TRUE, use.value.labels = FALSE
)

weights <- weights %>%
  transmute(Pingault_FID = clean_id(PINGAULT_FID),
            FOVWT2, PTTYPE2, BDSTRA00, SPTN00) %>%
  filter(Pingault_FID %in% myData$Pingault_FID)

# A duplicate family ID would multiply child records during the merge.
stopifnot(!anyDuplicated(weights$Pingault_FID))
myData <- left_join(myData, weights, by = "Pingault_FID")

# ============================================================
#          11. SAVE FINAL SEM DATASET
# ============================================================
stopifnot(!anyDuplicated(myData$CM_ID))
aframe_svy <- myData
cat("Children:", nrow(aframe_svy), "\n")
cat("Missing survey weights:", sum(is.na(aframe_svy$FOVWT2)), "\n")
cat("Weighted sample:", sum(!is.na(aframe_svy$FOVWT2)), "\n")
# Historical benchmarks: 3,228 children; 5 missing weights; 3,223 weighted.
save(aframe_svy, file = "NATCOMMS_Morosoli_2026.RData")
