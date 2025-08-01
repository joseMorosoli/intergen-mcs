#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   Fully constrained model with phenotypes (Reviewer Request)
# Purpose:  Fit fully constrained SEM models with no 
#           parent-specific or time-specific coefficients,
#           including parental phenotypic predictors 
#           (e.g., education).
#
# Outputs:  - Survey-weighted model fits (with & without PGS)
#           - Standardized estimates with confidence intervals
#           - Excel export of parameter estimates
#
# Usage:    Run this script to estimate fully constrained 
#           models incorporating parental phenotypes. 
#           Outputs can be used for reporting and comparison 
#           with previous models.
#
# Author:   Jose J. Morosoli
# Date:     31-07-2025
#############################################################

### --- Load Required Libraries ---
library(lavaan)
library(lme4)
library(psych)
library(tidyverse)
library(data.table)
library(reshape)
library(knitr)
library(kableExtra)
library(foreign)
library(survey)
library(lavaan.survey) # archived; must be installed manually
library(writexl)

### --- Helper Function: Extract Chi-Square ---
extract_chisq <- function(myanova, modelname) {
  chisq.diff <- myanova[2, "Chisq diff"]
  df.diff <- myanova[2, "Df diff"]
  pvalue <- myanova[2, "Pr(>Chisq)"]
  model_summary <- as.data.frame(cbind(modelname, chisq.diff, df.diff, pvalue))
  model_summary$pvalue_num <- as.numeric(as.character(model_summary$pvalue))
  model_summary$pvalue_sig[model_summary$pvalue_num <= 0.05 & model_summary$pvalue_num > 0.01] <- "*"
  model_summary$pvalue_sig[model_summary$pvalue_num <= 0.01] <- "**"
  model_summary$pvalue_sig[model_summary$pvalue_num > 0.05] <- "ns"
  return(model_summary)
}

options(stringsAsFactors = FALSE)

### --- Load Data ---
load(file = "NATCOMMS_R1.RData")

### --- Survey Design ---
aframe_svy_na <- aframe_svy[!is.na(aframe_svy$FOVWT2), ]
svy_design <- svydesign(
  ids = ~SPTN00,
  strata = ~PTTYPE2,
  weights = ~FOVWT2,
  data = aframe_svy_na,
  nest = TRUE
)

#############################################################
# MODEL SETUP
#############################################################

### --- Define Shared Objects ---
outcomes <- c(
  "BEXT_res", "CEXT_res", "DEXT_res", "FEXT_res",
  "BINT_res", "CINT_res", "DINT_res", "FINT_res"
)
traitG <- "lee"        # Genetic predictor
traitP <- "ADACAQ00"   # Parental phenotype (education)
modelname <- paste(traitG, "_", traitP, "_EXTINT_ph", sep = "")

ind <- "C"
ind2 <- "M"
ind3 <- "F"

# PGS variables
pgs1 <- paste0("EA_", traitG, "_", ind, "_res")
pgs2 <- paste0("EA_", traitG, "_", ind2, "_res")
pgs3 <- paste0("EA_", traitG, "_", ind3, "_res")

# Parental phenotype variables
phe2 <- paste0(traitP, "_", ind2)
phe3 <- paste0(traitP, "_", ind3)

#############################################################
# MODEL 1: With PGS + Parental Phenotypes
#############################################################

### --- Create Model Formulas ---
myformula_pheGen <- c(
  paste(outcomes[1], "~ be11*", pgs1, "+ be12*", pgs2, "+ be13*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[2], "~ be21*", pgs1, "+ be22*", pgs2, "+ be23*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[3], "~ be31*", pgs1, "+ be32*", pgs2, "+ be33*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[4], "~ be41*", pgs1, "+ be42*", pgs2, "+ be43*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[5], "~ bi11*", pgs1, "+ bi12*", pgs2, "+ bi13*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[6], "~ bi21*", pgs1, "+ bi22*", pgs2, "+ bi23*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[7], "~ bi31*", pgs1, "+ bi32*", pgs2, "+ bi33*", pgs3, "+", phe2, "+", phe3),
  paste(outcomes[8], "~ bi41*", pgs1, "+ bi42*", pgs2, "+ bi43*", pgs3, "+", phe2, "+", phe3)
)

### --- Combine with Correlated Residuals ---
model_pheGen <- paste(
  paste(myformula_pheGen, collapse = "\n"),
  "BINT_res ~~ CINT_res + DINT_res + FINT_res",
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res",
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "EA_lee_M_res ~~ EA_lee_F_res",
  "ADACAQ00_M ~~ ADACAQ00_F",
  sep = "\n"
)
model_pheGen_eaTidy <- noquote(strsplit(model_pheGen, "\n")[[1]])

### --- Fit the Model ---
fit_model_pheGen_ea <- sem(model_pheGen_eaTidy, data = aframe_svy_na, estimator = "ML", missing = "ML")
summary(fit_model_pheGen_ea)

### --- Apply Survey Design ---
fit_model_pheGen_ea_svy <- lavaan.survey(
  lavaan.fit = fit_model_pheGen_ea,
  survey.design = svy_design,
  estimator = "MLM"
)
summary(fit_model_pheGen_ea_svy)

#############################################################
# MODEL 2: Parental Phenotypes Only (No PGS)
#############################################################

### --- Create Model Formulas ---
myformula_phe <- sapply(outcomes, function(o) paste(o, "~", phe2, "+", phe3))

### --- Combine with Correlated Residuals ---
model_phe <- paste(
  paste(myformula_phe, collapse = "\n"),
  "BINT_res ~~ CINT_res + DINT_res + FINT_res",
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res",
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "ADACAQ00_M ~~ ADACAQ00_F",
  sep = "\n"
)
model_phe_eaTidy <- noquote(strsplit(model_phe, "\n")[[1]])

### --- Fit the Model ---
phe_ea_fit <- sem(model_phe_eaTidy, data = aframe_svy_na, estimator = "ML", missing = "ML")

### --- Apply Survey Design ---
phe_ea_fit_svy <- lavaan.survey(
  lavaan.fit = phe_ea_fit,
  survey.design = svy_design,
  estimator = "MLM"
)
summary(phe_ea_fit_svy, standardize = TRUE, rsquare = TRUE)

#############################################################
# EXPORT RESULTS
#############################################################

### --- Extract Standardized Estimates ---
estimates_cistd <- standardizedSolution(fit_model_pheGen_ea_svy)

### --- Select & Rename Key Columns ---
table_data <- estimates_cistd[, c("lhs", "op", "rhs", "est.std", "ci.lower", "ci.upper", "se")]
names(table_data) <- c("outcome", "operator", "predictor", "beta.std", "ci.lower", "ci.upper", "se")

### --- Export to Excel ---
write_xlsx(table_data, path = "pheGen_ea_svy_cist.xlsx")
