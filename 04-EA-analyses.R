#########################################
# Trio SEM Models: Educational Attainment
# Author: Jose J. Morosoli
# Date: 2025-08-01
# Purpose: Fit full set of Trio SEM models 
# (null, base, and constrained) with 
# complex survey adjustment, and run model 
# comparisons with FDR correction.
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
pgs_child <- paste0("EA_", trait, "_C_res")
pgs_mother <- paste0("EA_", trait, "_M_res")
pgs_father <- paste0("EA_", trait, "_F_res")

### --- Define Shared Covariances ---
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
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep = "\n"
)

### --- Helper: Build Model String ---
make_model <- function(formulas) {
  paste(paste(formulas, collapse = "\n"), covariances, sep = "\n")
}

#########################################
# MODELS
#########################################

### --- 0. Null Models ---
# A. Externalising only
formulas_nullext <- c(
  sapply(outcomes[1:4], function(outcome) paste0(outcome, "~ 0*", pgs_child, "+0*", pgs_mother, "+0*", pgs_father)),
  sapply(outcomes[5:8], function(outcome) paste0(outcome, "~ bi11a*", pgs_child, "+bi12a*", pgs_mother, "+bi13a*", pgs_father))
)
model_nullext <- make_model(formulas_nullext)
nullext_fit <- sem(model_nullext, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
nullext_fit_svy <- lavaan.survey(lavaan.fit = nullext_fit, survey.design = svy_design, estimator = "MLM")

# B. Internalising only
formulas_nullint <- c(
  sapply(outcomes[1:4], function(outcome) paste0(outcome, "~ be11a*", pgs_child, "+be12a*", pgs_mother, "+be13a*", pgs_father)),
  sapply(outcomes[5:8], function(outcome) paste0(outcome, "~ 0*", pgs_child, "+0*", pgs_mother, "+0*", pgs_father))
)
model_nullint <- make_model(formulas_nullint)
nullint_fit <- sem(model_nullint, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
nullint_fit_svy <- lavaan.survey(lavaan.fit = nullint_fit, survey.design = svy_design, estimator = "MLM")

# C. Both Externalising & Internalising
formulas_nullboth <- sapply(outcomes, function(outcome) paste0(outcome, "~ 0*", pgs_child, "+0*", pgs_mother, "+0*", pgs_father))
model_nullboth <- make_model(formulas_nullboth)
nullboth_fit <- sem(model_nullboth, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
nullboth_fit_svy <- lavaan.survey(lavaan.fit = nullboth_fit, survey.design = svy_design, estimator = "MLM")

### --- Base Model (Fully Free) ---
formulas_base <- sapply(outcomes, function(outcome) paste0(outcome, "~ ", pgs_child, "+", pgs_mother, "+", pgs_father))
model_base <- make_model(formulas_base)
base_fit <- sem(model_base, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
base_fit_svy <- lavaan.survey(lavaan.fit = base_fit, survey.design = svy_design, estimator = "MLM")

### --- Constrained Models ---
# ea1: Domain-invariant
formulas_ea1 <- sapply(outcomes, function(outcome) paste0(outcome, "~ b11*", pgs_child, "+b12*", pgs_mother, "+b13*", pgs_father))
model_ea1 <- make_model(formulas_ea1)
ea1_fit <- sem(model_ea1, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea1_fit_svy <- lavaan.survey(lavaan.fit = ea1_fit, survey.design = svy_design, estimator = "MLM")

# ea2: Source-invariant (Externalising)
formulas_ea2 <- c(
  sapply(outcomes[1:4], function(outcome) paste0(outcome, "~ be11*", pgs_child, "+be12*", pgs_mother, "+be12*", pgs_father)),
  sapply(outcomes[5:8], function(outcome) paste0(outcome, "~ bi11*", pgs_child, "+bi12*", pgs_mother, "+bi13*", pgs_father))
)
model_ea2 <- make_model(formulas_ea2)
ea2_fit <- sem(model_ea2, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea2_fit_svy <- lavaan.survey(lavaan.fit = ea2_fit, survey.design = svy_design, estimator = "MLM")

# ea3: Source-invariant (Internalising)
formulas_ea3 <- c(
  sapply(outcomes[1:4], function(outcome) paste0(outcome, "~ be11*", pgs_child, "+be12*", pgs_mother, "+be13*", pgs_father)),
  sapply(outcomes[5:8], function(outcome) paste0(outcome, "~ bi11*", pgs_child, "+bi12*", pgs_mother, "+bi12*", pgs_father))
)
model_ea3 <- make_model(formulas_ea3)
ea3_fit <- sem(model_ea3, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea3_fit_svy <- lavaan.survey(lavaan.fit = ea3_fit, survey.design = svy_design, estimator = "MLM")

# ea6: Time-invariant (Externalising)
formulas_ea6 <- c(
  sapply(outcomes[1:4], function(outcome) paste0(outcome, "~ be11*", pgs_child, "+be12*", pgs_mother, "+be13*", pgs_father)),
  sapply(outcomes[5:8], function(outcome) paste0(outcome, "~ bi11*", pgs_child, "+bi12*", pgs_mother, "+bi13*", pgs_father))
)
model_ea6 <- make_model(formulas_ea6)
ea6_fit <- sem(model_ea6, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea6_fit_svy <- lavaan.survey(lavaan.fit = ea6_fit, survey.design = svy_design, estimator = "MLM")

# ea7: Time-invariant (Internalising)
formulas_ea7 <- formulas_ea6 # same pattern as ea6
model_ea7 <- make_model(formulas_ea7)
ea7_fit <- sem(model_ea7, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea7_fit_svy <- lavaan.survey(lavaan.fit = ea7_fit, survey.design = svy_design, estimator = "MLM")

# ea9: Time + Source-invariant (Externalising)
formulas_ea9 <- formulas_ea2 # reuse pattern
model_ea9 <- make_model(formulas_ea9)
ea9_fit <- sem(model_ea9, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea9_fit_svy <- lavaan.survey(lavaan.fit = ea9_fit, survey.design = svy_design, estimator = "MLM")

# ea10: Time + Source-invariant (Internalising)
formulas_ea10 <- formulas_ea3 # reuse pattern
model_ea10 <- make_model(formulas_ea10)
ea10_fit <- sem(model_ea10, data = aframe_svy_na, estimator = "MLR", missing = "FIML")
ea10_fit_svy <- lavaan.survey(lavaan.fit = ea10_fit, survey.design = svy_design, estimator = "MLM")

#########################################
# MODEL COMPARISONS
#########################################

### --- Cluster 1: Null vs Base ---
Chisq_0aa <- extract_chisq(anova(nullext_fit_svy, base_fit_svy), "null_ea_ext")
Chisq_0ab <- extract_chisq(anova(nullint_fit_svy, base_fit_svy), "null_ea_int")
comparison0 <- rbind(Chisq_0aa, Chisq_0ab)

### --- Cluster 2: Constrained Models vs Base ---
Chisq_2 <- extract_chisq(anova(base_fit_svy, ea2_fit_svy), "EA_parentEq_EXT")
Chisq_3 <- extract_chisq(anova(base_fit_svy, ea3_fit_svy), "EA_parentEq_INT")
Chisq_4 <- extract_chisq(anova(base_fit_svy, ea6_fit_svy), "EA_timeinv_EXT")
Chisq_5 <- extract_chisq(anova(base_fit_svy, ea7_fit_svy), "EA_timeinv_INT")
Chisq_6 <- extract_chisq(anova(ea2_fit_svy, ea9_fit_svy), "EA_timeparentEq_EXT")
Chisq_7 <- extract_chisq(anova(ea3_fit_svy, ea10_fit_svy), "EA_timeparentEq_INT")

