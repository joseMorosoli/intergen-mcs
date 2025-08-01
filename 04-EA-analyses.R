#########################################
# Trio SEM Models: Educational Attainment
# Author: Jose J. Morosoli
# Date: 31-07-2025
#########################################

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
aframe_svy_na <- aframe_svy[!is.na(aframe_svy$FOVWT2), ]

### --- Define Complex Survey Design ---
svy_design <- svydesign(
  ids = ~SPTN00,
  strata = ~PTTYPE2,
  weights = ~FOVWT2,
  data = aframe_svy_na,
  nest = TRUE
)

### --- Define Outcomes & Predictors ---
outcomes <- c("BEXT_res", "CEXT_res", "DEXT_res", "FEXT_res",
              "BINT_res", "CINT_res", "DINT_res", "FINT_res")
trait <- "lee"
modelname <- paste0(trait, "_EXTINT")
pgs_child <- paste0("EA_", trait, "_C_res")
pgs_mother <- paste0("EA_", trait, "_M_res")
pgs_father <- paste0("EA_", trait, "_F_res")

### --- Quick Checks ---
describe(aframe_svy_na[, c(33:52, 89:132, 146:154, 157:164)])
table(aframe_svy_na$sex)

#########################################
# MODELS
#########################################

### --- 0. Null Models (Constrained to 0) ---
# Externalising-only
myformula_nullext <- sapply(outcomes[1:4], function(outcome) paste0(outcome, "~ 0*", pgs_child, "+ 0*", pgs_mother, "+ 0*", pgs_father))
myformula_nullext <- c(
  myformula_nullext,
  sapply(outcomes[5:8], function(outcome) paste0(outcome, "~ bi11a*", pgs_child, "+ bi12a*", pgs_mother, "+ bi13a*", pgs_father))
)
model_nullext <- paste(myformula_nullext,
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
                       "EA_lee_M_res ~~ EA_lee_F_res", sep = "\n")
model_nullextTidy <- noquote(strsplit(model_nullext, "\n")[[1]])
nullext_fit <- sem(model_nullextTidy, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
nullext_fit_svy <- lavaan.survey(lavaan.fit = nullext_fit, survey.design = svy_design, estimator = "MLM")

# (Repeat same structure for nullint_ea and null_ea – left unchanged)

### --- Base Model (Fully Free) ---
myformula_base <- sapply(outcomes, function(outcome) paste0(outcome, "~ ", pgs_child, "+", pgs_mother, "+", pgs_father))
model_base <- paste(myformula_base,
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
                    "EA_lee_M_res ~~ EA_lee_F_res", sep = "\n")
model_baseTidy <- noquote(strsplit(model_base, "\n")[[1]])
base_fit <- sem(model_baseTidy, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
base_fit_svy <- lavaan.survey(lavaan.fit = base_fit, survey.design = svy_design, estimator = "MLM")

### --- Other Models ---
# (ea1_fit_svy, ea2_fit_svy, ea3_fit_svy, ea5_fit_svy, ea6_fit_svy, ea7_fit_svy, ea8_fit_svy, ea9_fit_svy, ea10_fit_svy)
# These follow the same pattern as above and are kept unchanged for reproducibility.

#########################################
# MODEL COMPARISONS
#########################################

### --- Cluster 1: Null vs Base ---
model <- anova(nullext_ea_fit_svy, base_ea_fit_svy)
Chisq_0aa <- extract_chisq(model, modelname = "null_ea_ext")

model <- anova(nullint_ea_fit_svy, base_ea_fit_svy)
Chisq_0ab <- extract_chisq(model, modelname = "null_ea_int")

model <- anova(null_ea_fit_svy, base_ea_fit_svy)
Chisq_0ac <- extract_chisq(model, modelname = "null_ea_both")


### --- Cluster 2: Source / Time Constraints ---
model <- anova(base_ea_fit_svy, ea2_fit_svy)
Chisq_2 <- extract_chisq(model, modelname = "EA_parentEq_EXT")

model <- anova(base_ea_fit_svy, ea3_fit_svy)
Chisq_3 <- extract_chisq(model, modelname = "EA_parentEq_INT")

model <- anova(base_ea_fit_svy, ea6_fit_svy)
Chisq_4 <- extract_chisq(model, modelname = "EA_timeinv_EXT")

model <- anova(base_ea_fit_svy, ea7_fit_svy)
Chisq_5 <- extract_chisq(model, modelname = "EA_timeinv_INT")

model <- anova(ea2_fit_svy, ea9_fit_svy)
Chisq_6 <- extract_chisq(model, modelname = "EA_timeparentEq_EXT")

model <- anova(ea3_fit_svy, ea10_fit_svy)
Chisq_7 <- extract_chisq(model, modelname = "EA_timeparentEq_INT")
