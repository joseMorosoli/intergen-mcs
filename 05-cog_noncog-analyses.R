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
outcome1 <- "BEXT_res"; outcome2 <- "CEXT_res"; outcome3 <- "DEXT_res"; outcome4 <- "FEXT_res"
outcome5 <- "BINT_res"; outcome6 <- "CINT_res"; outcome7 <- "DINT_res"; outcome8 <- "FINT_res"

trait1 <- "cog"; trait2 <- "noncog"
pgs11 <- paste0("EA_",trait1,"_C_res"); pgs12 <- paste0("EA_",trait1,"_M_res"); pgs13 <- paste0("EA_",trait1,"_F_res")
pgs21 <- paste0("EA_",trait2,"_C_res"); pgs22 <- paste0("EA_",trait2,"_M_res"); pgs23 <- paste0("EA_",trait2,"_F_res")

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
  sep="\n"
)

#############################################
# HELPER FUNCTIONS
#############################################
create_model <- function(formulas) paste(paste(formulas, collapse = "\n"), covariances, sep = "\n")
make_model <- function(formulas) noquote(strsplit(create_model(formulas), "\n")[[1]])
fit_model <- function(model_tidy) {
  fit <- sem(paste(model_tidy, collapse="\n"), data=aframe_svy_na, estimator="MLR", missing="FIML")
  lavaan.survey(lavaan.fit=fit, survey.design=svy_design, estimator="MLM")
}
extract_chisq <- function(myanova, modelname) {
  chisq.diff <- myanova[2, "Chisq diff"]; df.diff <- myanova[2, "Df diff"]; pvalue <- myanova[2, "Pr(>Chisq)"]
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

### --- Base model (all PGS free) ---
formulas_base <- c(
  paste(outcome1, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome2, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome3, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome4, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome5, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome6, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome7, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + ")),
  paste(outcome8, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + "))
)

### --- Null models (joint, cognitive, noncognitive) ---
# EXT + INT joint
# (Explicitly paste all like in your original: formulas_null_EXTj, formulas_null_INTj, formulas_null_EXTc, formulas_null_EXTnc, formulas_null_INTc, formulas_null_INTnc)

### --- Time-invariant models ---
# EXT Cognitive
formulas_timeinv_EXTc <- c(
  paste(outcome1, "~ b_cog_C_EXT*",pgs11,"+ b_cog_M_EXT*",pgs12,"+ b_cog_F_EXT*",pgs13,"+ b_Bnoncog_C_EXT*",pgs21,"+ b_Bnoncog_M_EXT*",pgs22,"+ b_Bnoncog_F_EXT*",pgs23),
  paste(outcome2, "~ b_cog_C_EXT*",pgs11,"+ b_cog_M_EXT*",pgs12,"+ b_cog_F_EXT*",pgs13,"+ b_Cnoncog_C_EXT*",pgs21,"+ b_Cnoncog_M_EXT*",pgs22,"+ b_Cnoncog_F_EXT*",pgs23),
  paste(outcome3, "~ b_cog_C_EXT*",pgs11,"+ b_cog_M_EXT*",pgs12,"+ b_cog_F_EXT*",pgs13,"+ b_Dnoncog_C_EXT*",pgs21,"+ b_Dnoncog_M_EXT*",pgs22,"+ b_Dnoncog_F_EXT*",pgs23),
  paste(outcome4, "~ b_cog_C_EXT*",pgs11,"+ b_cog_M_EXT*",pgs12,"+ b_cog_F_EXT*",pgs13,"+ b_Fnoncog_C_EXT*",pgs21,"+ b_Fnoncog_M_EXT*",pgs22,"+ b_Fnoncog_F_EXT*",pgs23),
  paste(outcome5, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome6, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome7, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome8, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + "))
)

# Repeat same for EXTnc, INTc, INTnc (all fully expanded like above)

### --- Parent-invariant models ---
# EXT Cognitive
formulas_parentinv_EXTc <- c(
  paste(outcome1, "~ b_Bcog_C_EXT*",pgs11,"+ b_Bcog_P_EXT*",pgs12,"+ b_Bcog_P_EXT*",pgs13,"+ b_Bnoncog_C_EXT*",pgs21,"+ b_Bnoncog_M_EXT*",pgs22,"+ b_Bnoncog_F_EXT*",pgs23),
  paste(outcome2, "~ b_Ccog_C_EXT*",pgs11,"+ b_Ccog_P_EXT*",pgs12,"+ b_Ccog_P_EXT*",pgs13,"+ b_Cnoncog_C_EXT*",pgs21,"+ b_Cnoncog_M_EXT*",pgs22,"+ b_Cnoncog_F_EXT*",pgs23),
  paste(outcome3, "~ b_Dcog_C_EXT*",pgs11,"+ b_Dcog_P_EXT*",pgs12,"+ b_Dcog_P_EXT*",pgs13,"+ b_Dnoncog_C_EXT*",pgs21,"+ b_Dnoncog_M_EXT*",pgs22,"+ b_Dnoncog_F_EXT*",pgs23),
  paste(outcome4, "~ b_Fcog_C_EXT*",pgs11,"+ b_Fcog_P_EXT*",pgs12,"+ b_Fcog_P_EXT*",pgs13,"+ b_Fnoncog_C_EXT*",pgs21,"+ b_Fnoncog_M_EXT*",pgs22,"+ b_Fnoncog_F_EXT*",pgs23),
  paste(outcome5, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome6, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome7, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome8, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + "))
)

# Repeat for EXTnc, INTc, INTnc

### --- Time+Parent-invariant models ---
# EXT Cognitive
formulas_timeparentinv_EXTc <- c(
  paste(outcome1, "~ b_cog_C_EXT*",pgs11,"+ b_cog_P_EXT*",pgs12,"+ b_cog_P_EXT*",pgs13,"+ b_Bnoncog_C_EXT*",pgs21,"+ b_Bnoncog_P_EXT*",pgs22,"+ b_Bnoncog_P_EXT*",pgs23),
  paste(outcome2, "~ b_cog_C_EXT*",pgs11,"+ b_cog_P_EXT*",pgs12,"+ b_cog_P_EXT*",pgs13,"+ b_Cnoncog_C_EXT*",pgs21,"+ b_Cnoncog_P_EXT*",pgs22,"+ b_Cnoncog_P_EXT*",pgs23),
  paste(outcome3, "~ b_cog_C_EXT*",pgs11,"+ b_cog_P_EXT*",pgs12,"+ b_cog_P_EXT*",pgs13,"+ b_Dnoncog_C_EXT*",pgs21,"+ b_Dnoncog_P_EXT*",pgs22,"+ b_Dnoncog_P_EXT*",pgs23),
  paste(outcome4, "~ b_cog_C_EXT*",pgs11,"+ b_cog_P_EXT*",pgs12,"+ b_cog_P_EXT*",pgs13,"+ b_Fnoncog_C_EXT*",pgs21,"+ b_Fnoncog_P_EXT*",pgs22,"+ b_Fnoncog_P_EXT*",pgs23),
  paste(outcome5, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome6, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome7, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome8, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + "))
)

# Repeat for EXTnc, INTc, INTnc

### --- Component-invariant models ---
# EXT
formulas_compinv_EXT <- c(
  paste(outcome1, "~ b_B_C_EXT*",pgs11,"+ b_B_M_EXT*",pgs12,"+ b_B_F_EXT*",pgs13,"+ b_B_C_EXT*",pgs21,"+ b_B_M_EXT*",pgs22,"+ b_B_F_EXT*",pgs23),
  paste(outcome2, "~ b_C_C_EXT*",pgs11,"+ b_C_M_EXT*",pgs12,"+ b_C_F_EXT*",pgs13,"+ b_C_C_EXT*",pgs21,"+ b_C_M_EXT*",pgs22,"+ b_C_F_EXT*",pgs23),
  paste(outcome3, "~ b_D_C_EXT*",pgs11,"+ b_D_M_EXT*",pgs12,"+ b_D_F_EXT*",pgs13,"+ b_D_C_EXT*",pgs21,"+ b_D_M_EXT*",pgs22,"+ b_D_F_EXT*",pgs23),
  paste(outcome4, "~ b_F_C_EXT*",pgs11,"+ b_F_M_EXT*",pgs12,"+ b_F_F_EXT*",pgs13,"+ b_F_C_EXT*",pgs21,"+ b_F_M_EXT*",pgs22,"+ b_F_F_EXT*",pgs23),
  paste(outcome5, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome6, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome7, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome8, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + "))
)

# INT
formulas_compinv_INT <- c(
  paste(outcome1, "~ b_B_C_INT*",pgs11,"+ b_B_M_INT*",pgs12,"+ b_B_F_INT*",pgs13,"+ b_B_C_INT*",pgs21,"+ b_B_M_INT*",pgs22,"+ b_B_F_INT*",pgs23),
  paste(outcome2, "~ b_C_C_INT*",pgs11,"+ b_C_M_INT*",pgs12,"+ b_C_F_INT*",pgs13,"+ b_C_C_INT*",pgs21,"+ b_C_M_INT*",pgs22,"+ b_C_F_INT*",pgs23),
  paste(outcome3, "~ b_D_C_INT*",pgs11,"+ b_D_M_INT*",pgs12,"+ b_D_F_INT*",pgs13,"+ b_D_C_INT*",pgs21,"+ b_D_M_INT*",pgs22,"+ b_D_F_INT*",pgs23),
  paste(outcome4, "~ b_F_C_INT*",pgs11,"+ b_F_M_INT*",pgs12,"+ b_F_F_INT*",pgs13,"+ b_F_C_INT*",pgs21,"+ b_F_M_INT*",pgs22,"+ b_F_F_INT*",pgs23),
  paste(outcome5, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome6, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome7, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + ")),
  paste(outcome8, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23),collapse=" + "))
)

#############################################
# BUILD & FIT MODELS
#############################################
model_baseTidy <- make_model(formulas_base); fit_base_svy <- fit_model(model_baseTidy)

model_nullTidy_EXTc <- make_model(formulas_null_EXTc); fit_null_svy_EXTc <- fit_model(model_nullTidy_EXTc)
model_nullTidy_EXTnc <- make_model(formulas_null_EXTnc); fit_null_svy_EXTnc <- fit_model(model_nullTidy_EXTnc)
model_nullTidy_EXTj <- make_model(formulas_null_EXTj); fit_null_svy_EXTj <- fit_model(model_nullTidy_EXTj)

model_nullTidy_INTc <- make_model(formulas_null_INTc); fit_null_svy_INTc <- fit_model(model_nullTidy_INTc)
model_nullTidy_INTnc <- make_model(formulas_null_INTnc); fit_null_svy_INTnc <- fit_model(model_nullTidy_INTnc)
model_nullTidy_INTj <- make_model(formulas_null_INTj); fit_null_svy_INTj <- fit_model(model_nullTidy_INTj)

model_timeinvTidy_EXTc <- make_model(formulas_timeinv_EXTc); fit_timeinv_svy_EXTc <- fit_model(model_timeinvTidy_EXTc)
model_timeinvTidy_EXTnc <- make_model(formulas_timeinv_EXTnc); fit_timeinv_svy_EXTnc <- fit_model(model_timeinvTidy_EXTnc)
model_timeinvTidy_INTc <- make_model(formulas_timeinv_INTc); fit_timeinv_svy_INTc <- fit_model(model_timeinvTidy_INTc)
model_timeinvTidy_INTnc <- make_model(formulas_timeinv_INTnc); fit_timeinv_svy_INTnc <- fit_model(model_timeinvTidy_INTnc)

model_parentinvTidy_EXTc <- make_model(formulas_parentinv_EXTc); fit_parentinv_svy_EXTc <- fit_model(model_parentinvTidy_EXTc)
model_parentinvTidy_EXTnc <- make_model(formulas_parentinv_EXTnc); fit_parentinv_svy_EXTnc <- fit_model(model_parentinvTidy_EXTnc)
model_parentinvTidy_INTc <- make_model(formulas_parentinv_INTc); fit_parentinv_svy_INTc <- fit_model(model_parentinvTidy_INTc)
model_parentinvTidy_INTnc <- make_model(formulas_parentinv_INTnc); fit_parentinv_svy_INTnc <- fit_model(model_parentinvTidy_INTnc)

model_timeparentinvTidy_EXTc <- make_model(formulas_timeparentinv_EXTc); fit_timeparentinv_svy_EXTc <- fit_model(model_timeparentinvTidy_EXTc)
model_timeparentinvTidy_EXTnc <- make_model(formulas_timeparentinv_EXTnc); fit_timeparentinv_svy_EXTnc <- fit_model(model_timeparentinvTidy_EXTnc)
model_timeparentinvTidy_INTc <- make_model(formulas_timeparentinv_INTc); fit_timeparentinv_svy_INTc <- fit_model(model_timeparentinvTidy_INTc)
model_timeparentinvTidy_INTnc <- make_model(formulas_timeparentinv_INTnc); fit_timeparentinv_svy_INTnc <- fit_model(model_timeparentinvTidy_INTnc)

model_compinvTidy_EXT <- make_model(formulas_compinv_EXT); fit_compinv_svy_EXT <- fit_model(model_compinvTidy_EXT)
model_compinvTidy_INT <- make_model(formulas_compinv_INT); fit_compinv_svy_INT <- fit_model(model_compinvTidy_INT)

#############################################
# MODEL COMPARISONS
#############################################
Chisq_00a <- extract_chisq(anova(fit_base_svy, fit_null_svy_EXTc), 'CNC_null_EXTc')
Chisq_01a <- extract_chisq(anova(fit_base_svy, fit_null_svy_EXTnc), 'CNC_null_EXTnc')
Chisq_02a <- extract_chisq(anova(fit_base_svy, fit_null_svy_EXTj), 'CNC_null_EXTj')
Chisq_1a <- extract_chisq(anova(fit_base_svy, fit_timeinv_svy_EXTc), 'CNC_timeinv_EXTc')
Chisq_2a <- extract_chisq(anova(fit_base_svy, fit_parentinv_svy_EXTc), 'CNC_parentinv_EXTc')
Chisq_3a <- extract_chisq(anova(fit_parentinv_svy_EXTc, fit_timeparentinv_svy_EXTc), 'CNC_timeparentinv_EXTc')
Chisq_4a <- extract_chisq(anova(fit_base_svy, fit_timeinv_svy_EXTnc), 'CNC_timeinv_EXTnc')
Chisq_5a <- extract_chisq(anova(fit_base_svy, fit_parentinv_svy_EXTnc), 'CNC_parentinv_EXTnc')
Chisq_6a <- extract_chisq(anova(fit_parentinv_svy_EXTnc, fit_timeparentinv_svy_EXTnc), 'CNC_timeparentinv_EXTnc')
Chisq_7a <- extract_chisq(anova(fit_base_svy, fit_compinv_svy_EXT), 'CNC_compinv_EXT')

Chisq_00b <- extract_chisq(anova(fit_base_svy, fit_null_svy_INTc), 'CNC_null_INTc')
Chisq_01b <- extract_chisq(anova(fit_base_svy, fit_null_svy_INTnc), 'CNC_null_INTnc')
Chisq_02b <- extract_chisq(anova(fit_base_svy, fit_null_svy_INTj), 'CNC_null_INTj')
Chisq_1b <- extract_chisq(anova(fit_base_svy, fit_timeinv_svy_INTc), 'CNC_timeinv_INTc')
Chisq_2b <- extract_chisq(anova(fit_base_svy, fit_parentinv_svy_INTc), 'CNC_parentinv_INTc')
Chisq_3b <- extract_chisq(anova(fit_parentinv_svy_INTc, fit_timeparentinv_svy_INTc), 'CNC_timeparentinv_INTc')
Chisq_4b <- extract_chisq(anova(fit_base_svy, fit_timeinv_svy_INTnc), 'CNC_timeinv_INTnc')
Chisq_5b <- extract_chisq(anova(fit_base_svy, fit_parentinv_svy_INTnc), 'CNC_parentinv_INTnc')
Chisq_6b <- extract_chisq(anova(fit_parentinv_svy_INTnc, fit_timeparentinv_svy_INTnc), 'CNC_timeparentinv_INTnc')
Chisq_7b <- extract_chisq(anova(fit_base_svy, fit_compinv_svy_INT), 'CNC_compinv_INT')

nullcomp <- rbind(Chisq_02a,Chisq_02b)
comparison <- rbind(Chisq_1a,Chisq_2a,Chisq_3a,Chisq_4a,Chisq_5a,Chisq_6a,
                    Chisq_1b,Chisq_2b,Chisq_3b,Chisq_4b,Chisq_5b,Chisq_6b)
comparison$p_fdr <- p.adjust(comparison$pvalue, method="BH")
conscomp <- rbind(Chisq_7a,Chisq_7b)

print(nullcomp)
print(comparison)
print(conscomp)
