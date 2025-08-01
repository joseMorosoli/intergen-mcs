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
outcome_vars <- c("BEXT_res","CEXT_res","DEXT_res","FEXT_res",
                  "BINT_res","CINT_res","DINT_res","FINT_res")
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

### Base model (all effects free)
formulas_base <- sapply(outcome_vars, function(o) paste(o, "~", paste(c(pgs11, pgs12, pgs13, pgs21, pgs22, pgs23), collapse=" + ")))

### Null models: Joint, Cog-only, Noncog-only (EXT & INT)
make_null <- function(pattern, zero_cog=FALSE, zero_noncog=FALSE) {
  sapply(outcome_vars, function(o) {
    if (grepl(pattern, o)) {
      paste(o, "~",
            paste(c(if (zero_cog) paste0("0*",pgs11) else pgs11,
                    if (zero_cog) paste0("0*",pgs12) else pgs12,
                    if (zero_cog) paste0("0*",pgs13) else pgs13,
                    if (zero_noncog) paste0("0*",pgs21) else pgs21,
                    if (zero_noncog) paste0("0*",pgs22) else pgs22,
                    if (zero_noncog) paste0("0*",pgs23) else pgs23),
                  collapse=" + "))
    } else {
      paste(o, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + "))
    }
  })
}
formulas_null_EXTj <- make_null("EXT", TRUE, TRUE)
formulas_null_INTj <- make_null("INT", TRUE, TRUE)
formulas_null_EXTc <- make_null("EXT", TRUE, FALSE)
formulas_null_EXTnc <- make_null("EXT", FALSE, TRUE)
formulas_null_INTc <- make_null("INT", TRUE, FALSE)
formulas_null_INTnc <- make_null("INT", FALSE, TRUE)

### Constrained models
make_constrained <- function(pattern, labels) {
  sapply(outcome_vars, function(o) {
    if (grepl(pattern, o)) {
      lab <- substr(o,1,1)
      paste(o, "~", paste(sprintf(labels, lab), collapse=" + "))
    } else {
      paste(o, "~", paste(c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23), collapse=" + "))
    }
  })
}
# Time-invariant
time_labels_cog <- c("b_cog_C_%s*%s","b_cog_M_%s*%s","b_cog_F_%s*%s","b_%snoncog_C_%s*%s","b_%snoncog_M_%s*%s","b_%snoncog_F_%s*%s")
time_labels_noncog <- c("b_%scog_C_%s*%s","b_%scog_M_%s*%s","b_%scog_F_%s*%s","b_noncog_C_%s*%s","b_noncog_M_%s*%s","b_noncog_F_%s*%s")
formulas_timeinv_EXTc <- make_constrained("EXT", sprintf(time_labels_cog,"EXT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_timeinv_EXTnc <- make_constrained("EXT", sprintf(time_labels_noncog,"EXT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_timeinv_INTc <- make_constrained("INT", sprintf(time_labels_cog,"INT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_timeinv_INTnc <- make_constrained("INT", sprintf(time_labels_noncog,"INT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))

# Parent-invariant & Time+Parent-invariant (EXT/INT)
parent_labels <- c("b_%scog_C_%s*%s","b_%scog_P_%s*%s","b_%scog_P_%s*%s","b_%snoncog_C_%s*%s","b_%snoncog_P_%s*%s","b_%snoncog_P_%s*%s")
formulas_parentinv_EXTc <- make_constrained("EXT", sprintf(parent_labels,"EXT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_parentinv_EXTnc <- make_constrained("EXT", sprintf(parent_labels,"EXT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_parentinv_INTc <- make_constrained("INT", sprintf(parent_labels,"INT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_parentinv_INTnc <- make_constrained("INT", sprintf(parent_labels,"INT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_timeparentinv_EXTc <- formulas_parentinv_EXTc
formulas_timeparentinv_EXTnc <- formulas_parentinv_EXTnc
formulas_timeparentinv_INTc <- formulas_parentinv_INTc
formulas_timeparentinv_INTnc <- formulas_parentinv_INTnc

# Component-invariant
comp_labels <- c("b_%s_C_%s*%s","b_%s_M_%s*%s","b_%s_F_%s*%s","b_%s_C_%s*%s","b_%s_M_%s*%s","b_%s_F_%s*%s")
formulas_compinv_EXT <- make_constrained("EXT", sprintf(comp_labels,"EXT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))
formulas_compinv_INT <- make_constrained("INT", sprintf(comp_labels,"INT", c(pgs11,pgs12,pgs13,pgs21,pgs22,pgs23)))

#############################################
# BUILD & FIT MODELS
#############################################
build_and_fit <- function(formulas) fit_model(make_model(formulas))
fit_base_svy <- build_and_fit(formulas_base)

# Nulls
fit_null_svy_EXTc <- build_and_fit(formulas_null_EXTc)
fit_null_svy_EXTnc <- build_and_fit(formulas_null_EXTnc)
fit_null_svy_EXTj <- build_and_fit(formulas_null_EXTj)
fit_null_svy_INTc <- build_and_fit(formulas_null_INTc)
fit_null_svy_INTnc <- build_and_fit(formulas_null_INTnc)
fit_null_svy_INTj <- build_and_fit(formulas_null_INTj)

# Constrained
fit_timeinv_svy_EXTc <- build_and_fit(formulas_timeinv_EXTc)
fit_timeinv_svy_EXTnc <- build_and_fit(formulas_timeinv_EXTnc)
fit_timeinv_svy_INTc <- build_and_fit(formulas_timeinv_INTc)
fit_timeinv_svy_INTnc <- build_and_fit(formulas_timeinv_INTnc)
fit_parentinv_svy_EXTc <- build_and_fit(formulas_parentinv_EXTc)
fit_parentinv_svy_EXTnc <- build_and_fit(formulas_parentinv_EXTnc)
fit_parentinv_svy_INTc <- build_and_fit(formulas_parentinv_INTc)
fit_parentinv_svy_INTnc <- build_and_fit(formulas_parentinv_INTnc)
fit_timeparentinv_svy_EXTc <- build_and_fit(formulas_timeparentinv_EXTc)
fit_timeparentinv_svy_EXTnc <- build_and_fit(formulas_timeparentinv_EXTnc)
fit_timeparentinv_svy_INTc <- build_and_fit(formulas_timeparentinv_INTc)
fit_timeparentinv_svy_INTnc <- build_and_fit(formulas_timeparentinv_INTnc)
fit_compinv_svy_EXT <- build_and_fit(formulas_compinv_EXT)
fit_compinv_svy_INT <- build_and_fit(formulas_compinv_INT)

#############################################
# MODEL COMPARISONS
#############################################
# Externalising
Chisq_00a <- extract_chisq(anova(fit_base_svy, fit_null_svy_EXTc), 'CNC_null_EXTc')
Chisq_01a <- extract_chisq(anova(fit_base_svy, fit_null_svy_EXTnc), 'CNC_null_EXTnc')
Chisq_02a <- extract_chisq(anova(fit_base_svy, fit_null_svy_EXTj), 'CNC_null_EXTj')
Chisq_1a <- extract_chisq(anova(fit_base_svy, fit_timeinv_svy_EXTc), 'CNC_timeinv_EXTc')
Chisq_2a <- extract_chisq(anova(fit_base_svy, fit_parentinv_svy_EXTc), 'CNC_parentinv_EXTc')
Chisq_3a <- extract_chisq(anova(fit_parentinv_svy_EXTc, fit_timeparentinv_svy_EXTc), 'CNC_timeparentinv_EXTc')
Chisq_4a <- extract_chisq(anova(fit_base_svy, fit_timeinv_svy_EXTnc), 'CNC_timeinv_EXTnc')
Chisq_5a <- extract_chisq(anova(fit_base_svy, fit_parentinv_svy_EXTnc), 'CNC_parentinv_EXTnc')
Chisq_6a <- extract_chisq(anova(fit_parentinv_svy_EXTnc, fit_timeparentinv_svy_EXTnc), 'CNC_timeparentinv_EXTnc')
Chisq_7a <- extract_chisq(anova(fit_base_svy, fit_compinv_svy_EXT), 'CNC_compinv_EXTc')

# Internalising
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

# Combine results
nullcomp <- rbind(Chisq_02a, Chisq_02b)
comparison <- rbind(Chisq_1a,Chisq_2a,Chisq_3a,Chisq_4a,Chisq_5a,Chisq_6a,
                    Chisq_1b,Chisq_2b,Chisq_3b,Chisq_4b,Chisq_5b,Chisq_6b)
comparison$p_fdr <- p.adjust(comparison$pvalue, method="BH")
conscomp <- rbind(Chisq_7a, Chisq_7b)

print(nullcomp); print(comparison); print(conscomp)
