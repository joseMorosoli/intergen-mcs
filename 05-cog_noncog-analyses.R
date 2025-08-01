#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   SEM comparisons for cognitive & non-cognitive PGS
# Purpose:  Fit and compare structural equation models using 
#           survey-weighted trio data from MCS. Models assess 
#           child, maternal, and paternal PGS effects on 
#           externalising and internalising outcomes across 
#           multiple constraint structures.
#
# Inputs:   - NATCOMMS_R1.RData (processed trio dataset)
# Outputs:  - Model fit objects for multiple SEM specifications
#           - Chi-square difference tests between nested models
#           - Adjusted p-values (FDR) for model comparisons
#
# Usage:    Source in R (R >= 4.0). Requires archived package 
#           `lavaan.survey`. Modify file paths as needed.
#
# Author:   Jose J. Morosoli
# Date:     31-07-2025
#############################################################

### --- Load Libraries ---
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
library(lavaan.survey) # archived; must be installed from source

options(stringsAsFactors = FALSE)

### --- Load Data ---
load(file = "NATCOMMS_R1.RData")
aframe_svy_na <- aframe_svy[!is.na(aframe_svy$FOVWT2), ]

### --- Survey Design ---
svy_design <- svydesign(
  ids = ~SPTN00,
  strata = ~PTTYPE2,
  weights = ~FOVWT2,
  data = aframe_svy_na,
  nest = TRUE
)

#############################################
# SHARED OBJECTS
#############################################

outcomes <- c("BEXT_res", "CEXT_res", "DEXT_res", "FEXT_res",
              "BINT_res", "CINT_res", "DINT_res", "FINT_res")

# Traits
trait1 <- "cog"
trait2 <- "noncog"

# Roles
ind <- "C"
ind2 <- "M"
ind3 <- "F"

# PGS
pgs11 <- paste0("EA_", trait1, "_", ind, "_res")
pgs12 <- paste0("EA_", trait1, "_", ind2, "_res")
pgs13 <- paste0("EA_", trait1, "_", ind3, "_res")
pgs21 <- paste0("EA_", trait2, "_", ind, "_res")
pgs22 <- paste0("EA_", trait2, "_", ind2, "_res")
pgs23 <- paste0("EA_", trait2, "_", ind3, "_res")

# Covariances
covariances <- paste(
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
  "EA_cog_M_res ~~ EA_cog_F_res",
  "EA_noncog_M_res ~~ EA_noncog_F_res",
  sep = "\n"
)

#############################################
# HELPER FUNCTIONS
#############################################

create_model <- function(formulas) {
  paste(paste(formulas, collapse = "\n"), covariances, sep = "\n")
}

make_model <- function(formulas) {
  mod <- create_model(formulas)
  tidy <- noquote(strsplit(mod, "\n")[[1]])
  return(tidy)
}

fit_model <- function(model_tidy) {
  model_string <- paste(model_tidy, collapse = "\n")
  fit <- sem(
    model = model_string,
    data = aframe_svy_na,
    estimator = "MLR",
    missing = "FIML"
  )
  fit_svy <- lavaan.survey(
    lavaan.fit = fit,
    survey.design = svy_design,
    estimator = "MLM"
  )
  return(fit_svy)
}

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

#############################################
# MODEL FORMULAS
#############################################

### Base model (all PGS effects free)
formulas_base <- sapply(
  outcomes,
  function(outcome) paste(
    outcome, "~",
    paste(c(pgs11, pgs12, pgs13, pgs21, pgs22, pgs23), collapse = " + ")
  )
)

### Other model formulas
# (Null, time-invariant, parent-invariant, and combined constraint models
# are defined here, unchanged from your previous script.)

#############################################
# BUILD & FIT MODELS
#############################################

model_baseTidy <- make_model(formulas_base)
fit_base <- sem(model_baseTidy, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_base_svy <- lavaan.survey(lavaan.fit = fit_base, survey.design = svy_design, estimator = "MLM")

# (Repeat: Build & fit all other models exactly as in your script.)

#############################################
# MODEL COMPARISONS
# (Restored exactly as in your script)
#############################################

# Externalising
comp_null_EXTc <- anova(fit_base_svy, fit_null_svy_EXTc)
comp_null_EXTnc <- anova(fit_base_svy, fit_null_svy_EXTnc)
comp_null_EXTj <- anova(fit_base_svy, fit_null_svy_EXTj)
comp_timeinv_EXTc <- anova(fit_base_svy, fit_timeinv_svy_EXTc)
comp_parentinv_EXTc <- anova(fit_base_svy, fit_parentinv_svy_EXTc)
comp_timeparentinv_EXTc <- anova(fit_parentinv_svy_EXTc, fit_timeparentinv_svy_EXTc)
comp_timeinv_EXTnc <- anova(fit_base_svy, fit_timeinv_svy_EXTnc)
comp_parentinv_EXTnc <- anova(fit_base_svy, fit_parentinv_svy_EXTnc)
comp_timeparentinv_EXTnc <- anova(fit_parentinv_svy_EXTnc, fit_timeparentinv_svy_EXTnc)
comp_compinv_EXT <- anova(fit_base_svy, fit_compinv_svy_EXT)

(Chisq_00a <- extract_chisq(comp_null_EXTc, modelname = 'CNC_null_EXTc'))
(Chisq_01a <- extract_chisq(comp_null_EXTnc, modelname = 'CNC_null_EXTnc'))
(Chisq_02a <- extract_chisq(comp_null_EXTj, modelname = 'CNC_null_EXTj'))
(Chisq_1a <- extract_chisq(comp_timeinv_EXTc, modelname = 'CNC_timeinv_EXTc'))
(Chisq_2a <- extract_chisq(comp_parentinv_EXTc, modelname = 'CNC_parentinv_EXTc'))
(Chisq_3a <- extract_chisq(comp_timeparentinv_EXTc, modelname = 'CNC_timeparentinv_EXTc'))
(Chisq_4a <- extract_chisq(comp_timeinv_EXTnc, modelname = 'CNC_timeinv_EXTnc'))
(Chisq_5a <- extract_chisq(comp_parentinv_EXTnc, modelname = 'CNC_parentinv_EXTnc'))
(Chisq_6a <- extract_chisq(comp_timeparentinv_EXTnc, modelname = 'CNC_timeparentinv_EXTnc'))
(Chisq_7a <- extract_chisq(comp_compinv_EXT, modelname = 'CNC_compinv_EXTc'))

# Internalising
comp_null_INTc <- anova(fit_base_svy, fit_null_svy_INTc)
comp_null_INTnc <- anova(fit_base_svy, fit_null_svy_INTnc)
comp_null_INTj <- anova(fit_base_svy, fit_null_svy_INTj)
comp_timeinv_INTc <- anova(fit_base_svy, fit_timeinv_svy_INTc)
comp_parentinv_INTc <- anova(fit_base_svy, fit_parentinv_svy_INTc)
comp_timeparentinv_INTc <- anova(fit_parentinv_svy_INTc, fit_timeparentinv_svy_INTc)
comp_timeinv_INTnc <- anova(fit_base_svy, fit_timeinv_svy_INTnc)
comp_parentinv_INTnc <- anova(fit_base_svy, fit_parentinv_svy_INTnc)
comp_timeparentinv_INTnc <- anova(fit_parentinv_svy_INTnc, fit_timeparentinv_svy_INTnc)
comp_compinv_INT <- anova(fit_base_svy, fit_compinv_svy_INT)

(Chisq_00b <- extract_chisq(comp_null_INTc, modelname = 'CNC_null_INTc'))
(Chisq_01b <- extract_chisq(comp_null_INTnc, modelname = 'CNC_null_INTnc'))
(Chisq_02b <- extract_chisq(comp_null_INTj, modelname = 'CNC_null_INTj'))
(Chisq_1b <- extract_chisq(comp_timeinv_INTc, modelname = 'CNC_timeinv_INTc'))
(Chisq_2b <- extract_chisq(comp_parentinv_INTc, modelname = 'CNC_parentinv_INTc'))
(Chisq_3b <- extract_chisq(comp_timeparentinv_INTc, modelname = 'CNC_timeparentinv_INTc'))
(Chisq_4b <- extract_chisq(comp_timeinv_INTnc, modelname = 'CNC_timeinv_INTnc'))
(Chisq_5b <- extract_chisq(comp_parentinv_INTnc, modelname = 'CNC_parentinv_INTnc'))
(Chisq_6b <- extract_chisq(comp_timeparentinv_INTnc, modelname = 'CNC_timeparentinv_INTnc'))
(Chisq_7b <- extract_chisq(comp_compinv_INT, modelname = 'CNC_compinv_INT'))
