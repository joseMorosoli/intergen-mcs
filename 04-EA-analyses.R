#########################################
# Trio SEM Models: Educational Attainment
# Author: Jose J. Morosoli
# Date: 31-07-2025
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
options(stringsAsFactors = FALSE)

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

### --- Load Data ---
load(file = "../NATCOMMS_R1.RData")
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
outcome1 <- "BEXT_res"
outcome2 <- "CEXT_res"
outcome3 <- "DEXT_res" 
outcome4 <- "FEXT_res" 
outcome5 <- "BINT_res" 
outcome6 <- "CINT_res" 
outcome7 <- "DINT_res" 
outcome8 <- "FINT_res"
# What predictor? 
trait <- 'lee'
# Label model and specify trio PGSs
modelname <- paste(trait,"_","EXTINT",sep="") #   EXT INT
ind <- "C"
ind2 <- "M"
ind3 <- "F"
# Educational attainment
pgs1 <- paste0("EA_",trait, "_", ind, "_res")
pgs2 <- paste0("EA_",trait, "_", ind2, "_res")
pgs3 <- paste0("EA_",trait, "_", ind3, "_res")

#-----------------------------------------------
# 0. NULL: All Parameters Constrained To Zero Model 
#-----------------------------------------------

# A. Externalising only
#formulas where trio PGS effects are set to ZERO for externalising
myformula_nullext_ea1 <- paste(outcome1, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullext_ea2 <- paste(outcome2, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullext_ea3 <- paste(outcome3, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullext_ea4 <- paste(outcome4, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullext_ea5 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi13a*",pgs3, sep="")
myformula_nullext_ea6 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi13b*",pgs3, sep="")
myformula_nullext_ea7 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi13c*",pgs3, sep="")
myformula_nullext_ea8 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi13d*",pgs3, sep="")

#specify the free model using these formulas
model_nullext_ea <- paste(# regressions 
  myformula_nullext_ea1,
  myformula_nullext_ea2,
  myformula_nullext_ea3,
  myformula_nullext_ea4,
  myformula_nullext_ea5,
  myformula_nullext_ea6,
  myformula_nullext_ea7,
  myformula_nullext_ea8,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_nullext_eaTidy <- noquote(strsplit(model_nullext_ea, "\n")[[1]])
# Fit the free model
nullext_ea_fit <- sem(model_nullext_eaTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
nullext_ea_fit_svy <- lavaan.survey(lavaan.fit = nullext_ea_fit, survey.design = svy_design, estimator = "MLM")


# B. Internalising only

#formulas where trio PGS effects are set to ZERO for internalising
myformula_nullint_ea1 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be13a*",pgs3, sep="")
myformula_nullint_ea2 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be13b*",pgs3, sep="")
myformula_nullint_ea3 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be13c*",pgs3, sep="")
myformula_nullint_ea4 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be13d*",pgs3, sep="")
myformula_nullint_ea5 <- paste(outcome5, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullint_ea6 <- paste(outcome6, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullint_ea7 <- paste(outcome7, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_nullint_ea8 <- paste(outcome8, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")

#specify the free model using these formulas
model_nullint_ea <- paste(# regressions 
  myformula_nullint_ea1,
  myformula_nullint_ea2,
  myformula_nullint_ea3,
  myformula_nullint_ea4,
  myformula_nullint_ea5,
  myformula_nullint_ea6,
  myformula_nullint_ea7,
  myformula_nullint_ea8,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_nullint_eaTidy <- noquote(strsplit(model_nullint_ea, "\n")[[1]])
# Fit the free model
nullint_ea_fit <- sem(model_nullint_eaTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
nullint_ea_fit_svy <- lavaan.survey(lavaan.fit = nullint_ea_fit, survey.design = svy_design, estimator = "MLM")

# C. Both Externalising and Internalising

#formulas where trio PGS effects are set to ZERO for both phenotypes
myformula_null_ea1 <- paste(outcome1, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea2 <- paste(outcome2, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea3 <- paste(outcome3, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea4 <- paste(outcome4, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea5 <- paste(outcome5, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea6 <- paste(outcome6, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea7 <- paste(outcome7, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")
myformula_null_ea8 <- paste(outcome8, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3, sep="")

#specify the free model using these formulas
model_null_ea <- paste(# regressions 
  myformula_null_ea1,
  myformula_null_ea2,
  myformula_null_ea3,
  myformula_null_ea4,
  myformula_null_ea5,
  myformula_null_ea6,
  myformula_null_ea7,
  myformula_null_ea8,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_null_eaTidy <- noquote(strsplit(model_null_ea, "\n")[[1]])
# Fit the free model
null_ea_fit <- sem(model_null_eaTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
null_ea_fit_svy <- lavaan.survey(lavaan.fit = null_ea_fit, survey.design = svy_design, estimator = "MLM")



#-----------------------------------------------
# BASE: Fully Free Model 
#-----------------------------------------------

#formulas where trio PGS effects are FREE to vary over time
myformula_base_ea1 <- paste(outcome1, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea2 <- paste(outcome2, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea3 <- paste(outcome3, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea4 <- paste(outcome4, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea5 <- paste(outcome5, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea6 <- paste(outcome6, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea7 <- paste(outcome7, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
myformula_base_ea8 <- paste(outcome8, "~",pgs1,"+",pgs2,"+",pgs3, sep="")
#specify the free model using these formulas
model_base_ea <- paste(# regressions 
  myformula_base_ea1,
  myformula_base_ea2,
  myformula_base_ea3,
  myformula_base_ea4,
  myformula_base_ea5,
  myformula_base_ea6,
  myformula_base_ea7,
  myformula_base_ea8,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_base_eaTidy <- noquote(strsplit(model_base_ea, "\n")[[1]])
# Fit the free model
base_ea_fit <- sem(model_base_eaTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
base_ea_fit_svy <- lavaan.survey(lavaan.fit = base_ea_fit, survey.design = svy_design, estimator = "MLM")


#-----------------------------------------------
# 1. DOMAIN-INVARIANT (EXT = INT)
#-----------------------------------------------

myformula_ea11 <- paste(outcome1, "~","b11a*",pgs1,"+","b12a*",pgs2,"+","b13a*",pgs3, sep="")
myformula_ea12 <- paste(outcome2, "~","b11b*",pgs1,"+","b12b*",pgs2,"+","b13b*",pgs3, sep="")
myformula_ea13 <- paste(outcome3, "~","b11c*",pgs1,"+","b12c*",pgs2,"+","b13c*",pgs3, sep="")
myformula_ea14 <- paste(outcome4, "~","b11d*",pgs1,"+","b12d*",pgs2,"+","b13d*",pgs3, sep="")
myformula_ea15 <- paste(outcome5, "~","b11a*",pgs1,"+","b12a*",pgs2,"+","b13a*",pgs3, sep="")
myformula_ea16 <- paste(outcome6, "~","b11b*",pgs1,"+","b12b*",pgs2,"+","b13b*",pgs3, sep="")
myformula_ea17 <- paste(outcome7, "~","b11c*",pgs1,"+","b12c*",pgs2,"+","b13c*",pgs3, sep="")
myformula_ea18 <- paste(outcome8, "~","b11d*",pgs1,"+","b12d*",pgs2,"+","b13d*",pgs3, sep="")

#specify the free model using these formulas
model_ea1 <- paste(# regressions 
  myformula_ea11,
  myformula_ea12,
  myformula_ea13,
  myformula_ea14,
  myformula_ea15,
  myformula_ea16,
  myformula_ea17,
  myformula_ea18,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea1Tidy <- noquote(strsplit(model_ea1, "\n")[[1]])

# Fit the free model
ea1_fit <- sem(model_ea1Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea1_fit_svy <- lavaan.survey(lavaan.fit = ea1_fit, survey.design = svy_design, estimator = "MLM")



#-----------------------------------------------
# 2. SOURCE-INVARIANT (Maternal = Paternal)
#-----------------------------------------------

# A. Externalising only (cog + noncog source-invariant)

myformula_ea21 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be12a*",pgs3, sep="")
myformula_ea22 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be12b*",pgs3, sep="")
myformula_ea23 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be12c*",pgs3, sep="")
myformula_ea24 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be12d*",pgs3, sep="")
myformula_ea25 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi13a*",pgs3, sep="")
myformula_ea26 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi13b*",pgs3, sep="")
myformula_ea27 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi13c*",pgs3, sep="")
myformula_ea28 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi13d*",pgs3, sep="")

#specify the free model using these formulas
model_ea2 <- paste(# regressions 
  myformula_ea21,
  myformula_ea22,
  myformula_ea23,
  myformula_ea24,
  myformula_ea25,
  myformula_ea26,
  myformula_ea27,
  myformula_ea28,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea2Tidy <- noquote(strsplit(model_ea2, "\n")[[1]])

# Fit the free model
ea2_fit <- sem(model_ea2Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea2_fit_svy <- lavaan.survey(lavaan.fit = ea2_fit, survey.design = svy_design, estimator = "MLM")


# B. Internalising only (cog + noncog source-invariant)

myformula_ea31 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be13a*",pgs3, sep="")
myformula_ea32 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be13b*",pgs3, sep="")
myformula_ea33 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be13c*",pgs3, sep="")
myformula_ea34 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be13d*",pgs3, sep="")
myformula_ea35 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi12a*",pgs3, sep="")
myformula_ea36 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi12b*",pgs3, sep="")
myformula_ea37 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi12c*",pgs3, sep="")
myformula_ea38 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi12d*",pgs3, sep="")

#specify the free model using these formulas
model_ea3 <- paste(# regressions 
  myformula_ea31,
  myformula_ea32,
  myformula_ea33,
  myformula_ea34,
  myformula_ea35,
  myformula_ea36,
  myformula_ea37,
  myformula_ea38,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea3Tidy <- noquote(strsplit(model_ea3, "\n")[[1]])
model_ea3Tidy

# Fit the free model
ea3_fit <- sem(model_ea3Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea3_fit_svy <- lavaan.survey(lavaan.fit = ea3_fit, survey.design = svy_design, estimator = "MLM")


# C. Both Externalising and Internalising (cog + noncog source-invariant)

myformula_ea51 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be12a*",pgs3, sep="")
myformula_ea52 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be12b*",pgs3, sep="")
myformula_ea53 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be12c*",pgs3, sep="")
myformula_ea54 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be12d*",pgs3, sep="")
myformula_ea55 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi12a*",pgs3, sep="")
myformula_ea56 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi12b*",pgs3, sep="")
myformula_ea57 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi12c*",pgs3, sep="")
myformula_ea58 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi12d*",pgs3, sep="")

#specify the free model using these formulas
model_ea5 <- paste(# regressions 
  myformula_ea51,
  myformula_ea52,
  myformula_ea53,
  myformula_ea54,
  myformula_ea55,
  myformula_ea56,
  myformula_ea57,
  myformula_ea58,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea5Tidy <- noquote(strsplit(model_ea5, "\n")[[1]])
model_ea5Tidy

# Fit the free model
ea5_fit <- sem(model_ea5Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea5_fit_svy <- lavaan.survey(lavaan.fit = ea5_fit, survey.design = svy_design, estimator = "MLM")

#-----------------------------------------------
# 3. TIME-INVARIANT MODELS
#-----------------------------------------------

# A. Externalising only (cog + noncog time-invariant)
#formulas where trio PGS effects are constrained across time
myformula_ea61 <- paste(outcome1, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea62 <- paste(outcome2, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea63 <- paste(outcome3, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea64 <- paste(outcome4, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea65 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi13a*",pgs3, sep="")
myformula_ea66 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi13b*",pgs3, sep="")
myformula_ea67 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi13c*",pgs3, sep="")
myformula_ea68 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi13d*",pgs3, sep="")

#specify the free model using these formulas
model_ea6 <- paste(# regressions 
  myformula_ea61,
  myformula_ea62,
  myformula_ea63,
  myformula_ea64,
  myformula_ea65,
  myformula_ea66,
  myformula_ea67,
  myformula_ea68,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea6Tidy <- noquote(strsplit(model_ea6, "\n")[[1]])

# Fit the free model
ea6_fit <- sem(model_ea6Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea6_fit_svy <- lavaan.survey(lavaan.fit = ea6_fit, survey.design = svy_design, estimator = "MLM")


# B. Internalising only (cog + noncog time-invariant)

#formulas where trio PGS effects are constrained across time
myformula_ea71 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be13a*",pgs3, sep="")
myformula_ea72 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be13b*",pgs3, sep="")
myformula_ea73 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be13c*",pgs3, sep="")
myformula_ea74 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be13d*",pgs3, sep="")
myformula_ea75 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
myformula_ea76 <- paste(outcome6, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
myformula_ea77 <- paste(outcome7, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
myformula_ea78 <- paste(outcome8, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
#specify the free model using these formulas
model_ea7 <- paste(# regressions 
  myformula_ea71,
  myformula_ea72,
  myformula_ea73,
  myformula_ea74,
  myformula_ea75,
  myformula_ea76,
  myformula_ea77,
  myformula_ea78,
  ## correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea7Tidy <- noquote(strsplit(model_ea7, "\n")[[1]])

# Fit the free model
ea7_fit <- sem(model_ea7Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea7_fit_svy <- lavaan.survey(lavaan.fit = ea7_fit, survey.design = svy_design, estimator = "MLM")

# C. Both Externalising and Internalising (cog + noncog time-invariant))
#formulas where trio PGS effects are constrained across time
myformula_ea81 <- paste(outcome1, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea82 <- paste(outcome2, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea83 <- paste(outcome3, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea84 <- paste(outcome4, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3, sep="")
myformula_ea85 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
myformula_ea86 <- paste(outcome6, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
myformula_ea87 <- paste(outcome7, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
myformula_ea88 <- paste(outcome8, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3, sep="")
#specify the free model using these formulas
model_ea8 <- paste(# regressions 
  myformula_ea81,
  myformula_ea82,
  myformula_ea83,
  myformula_ea84,
  myformula_ea85,
  myformula_ea86,
  myformula_ea87,
  myformula_ea88,
  ## correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea8Tidy <- noquote(strsplit(model_ea8, "\n")[[1]])

# Fit the free model
ea8_fit <- sem(model_ea8Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea8_fit_svy <- lavaan.survey(lavaan.fit = ea8_fit, survey.design = svy_design, estimator = "MLM")

#-----------------------------------------------
# 4. TIME- AND SOURCE-INVARIANT 
#-----------------------------------------------

# A. Externalising only
#formulas where trio PGS effects are constrained across time and parent
myformula_ea91 <- paste(outcome1, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea92 <- paste(outcome2, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea93 <- paste(outcome3, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea94 <- paste(outcome4, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea95 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi13a*",pgs3, sep="")
myformula_ea96 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi13b*",pgs3, sep="")
myformula_ea97 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi13c*",pgs3, sep="")
myformula_ea98 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi13d*",pgs3, sep="")

#specify the free model using these formulas
model_ea9 <- paste(# regressions 
  myformula_ea91,
  myformula_ea92,
  myformula_ea93,
  myformula_ea94,
  myformula_ea95,
  myformula_ea96,
  myformula_ea97,
  myformula_ea98,
  # correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea9Tidy <- noquote(strsplit(model_ea9, "\n")[[1]])

# Fit the free model
ea9_fit <- sem(model_ea9Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea9_fit_svy <- lavaan.survey(lavaan.fit = ea9_fit, survey.design = svy_design, estimator = "MLM")


# B. Internalising only
#formulas where trio PGS effects are constrained across time and parent
myformula_ea101 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be13a*",pgs3, sep="")
myformula_ea102 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be13b*",pgs3, sep="")
myformula_ea103 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be13c*",pgs3, sep="")
myformula_ea104 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be13d*",pgs3, sep="")
myformula_ea105 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea106 <- paste(outcome6, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea107 <- paste(outcome7, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea108 <- paste(outcome8, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
#specify the free model using these formulas
model_ea10 <- paste(# regressions 
  myformula_ea101,
  myformula_ea102,
  myformula_ea103,
  myformula_ea104,
  myformula_ea105,
  myformula_ea106,
  myformula_ea107,
  myformula_ea108,
  ## correlated residuals
  "BINT_res ~~ CINT_res + DINT_res + FINT_res", 
  "CINT_res ~~ DINT_res + FINT_res",
  "DINT_res ~~ FINT_res",
  "BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res", 
  "CEXT_res ~~ DEXT_res + FEXT_res",
  "DEXT_res ~~ FEXT_res",
  # correlated residuals between int and ext
  "BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  "CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res",
  "FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res", 
  # Correlated parental phenotypes
  "EA_lee_M_res ~~ EA_lee_F_res",
  sep="\n")

#remove quotation marks and separate formulae
model_ea10Tidy <- noquote(strsplit(model_ea10, "\n")[[1]])

# Fit the free model
ea10_fit <- sem(model_ea10Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
ea10_fit_svy <- lavaan.survey(lavaan.fit = ea10_fit, survey.design = svy_design, estimator = "MLM")

#########################################
# MODEL COMPARISONS
#########################################

### --- Cluster 1: Null vs Base ---
model <- anova(nullext_ea_fit_svy,base_ea_fit_svy)
(Chisq_0aa <- extract_chisq(model, modelname = 'null_ea_ext'))
model <- anova(nullint_ea_fit_svy,base_ea_fit_svy)
(Chisq_0ab <- extract_chisq(model, modelname = 'null_ea_int'))
model <- anova(null_ea_fit_svy,base_ea_fit_svy)
(Chisq_0ac <- extract_chisq(model, modelname = 'null_ea_both'))

### --- Cluster 2: Constrained Models vs Base ---
#~~~~~~~~~ FIRST STEP: PARENTS ~~~~~~~#
model <- anova(base_ea_fit_svy,ea2_fit_svy)
(Chisq_2 <- extract_chisq(model, modelname = 'EA_parentEq_EXT'))
model <- anova(base_ea_fit_svy,ea3_fit_svy)
(Chisq_3 <- extract_chisq(model, modelname = 'EA_parentEq_INT'))
#~~~~~~~~~ SECOND STEP: TIME ~~~~~~~#
model <- anova(base_ea_fit_svy,ea6_fit_svy)
(Chisq_4 <- extract_chisq(model, modelname = 'EA_timeinv_EXT'))
model <- anova(base_ea_fit_svy,ea7_fit_svy)
(Chisq_5 <- extract_chisq(model, modelname = 'EA_timeinv_INT'))
#~~~~~~~~~ THIRD STEP: PARENTS + TIME ~~~~~~~#
model <- anova(ea2_fit_svy,ea9_fit_svy)
(Chisq_6 <- extract_chisq(model, modelname = 'EA_timeparentEq_EXT'))
model <- anova(ea3_fit_svy,ea10_fit_svy)
(Chisq_7 <- extract_chisq(model, modelname = 'EA_timeparentEq_INT'))
