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

outcome1 <- "BEXT_res"
outcome2 <- "CEXT_res"
outcome3 <- "DEXT_res" 
outcome4 <- "FEXT_res" 
outcome5 <- "BINT_res" 
outcome6 <- "CINT_res" 
outcome7 <- "DINT_res" 
outcome8 <- "FINT_res"

trait1 <- "cog" 
trait2 <- "noncog" 

ind <- "C"
ind2 <- "M"
ind3 <- "F"
pgs11 <- paste0("EA_",trait1,"_",ind,"_res")
pgs12 <- paste0("EA_",trait1,"_",ind2,"_res")
pgs13 <- paste0("EA_",trait1,"_",ind3,"_res")
pgs21 <- paste0("EA_",trait2,"_",ind,"_res")
pgs22 <- paste0("EA_",trait2,"_",ind2,"_res")
pgs23 <- paste0("EA_",trait2,"_",ind3,"_res")

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
  sep="\n")

# ----------------------------------------------------------
# HELPERS
# ----------------------------------------------------------

create_model <- function(formulas) {
  paste(
    paste(formulas, collapse = "\n"),
    covariances,
    sep = "\n"
  )
}

make_model <- function(formulas) {
  mod <- create_model(formulas)
  tidy <- noquote(strsplit(mod, "\n")[[1]])
  return(tidy)
}

fit_model <- function(model_tidy) {
  model_string <- paste(model_tidy, collapse="\n")
  fit <- sem(model = model_string, 
             data = aframe_svy_na, 
             estimator = "MLR", 
             missing = "FIML")
  fit_svy <- lavaan.survey(lavaan.fit = fit, 
                           survey.design = svy_design, 
                           estimator = "MLM")
  return(fit_svy)
}

# ----------------------------------------------------------
# BASE MODEL (Fully Free)
# ----------------------------------------------------------

formulas_base <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    paste(
      outcome, "~",
      paste(
        c(
          pgs11,
          pgs12,
          pgs13,
          pgs21,
          pgs22,
          pgs23
        ),
        collapse = " + "
      )
    )
  }
)

# ----------------------------------------------------------
# NULL MODELS
# ----------------------------------------------------------

#### JOINT NULL MODELS ####

# Null model (separate for EXT)
formulas_null_EXTj <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      paste(
        outcome, "~",
        paste(
          c(
            paste0("0*", pgs11),
            paste0("0*", pgs12),
            paste0("0*", pgs13),
            paste0("0*", pgs21),
            paste0("0*", pgs22),
            paste0("0*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


# Null model (separate for INT)
formulas_null_INTj <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      paste(
        outcome, "~",
        paste(
          c(
            paste0("0*", pgs11),
            paste0("0*", pgs12),
            paste0("0*", pgs13),
            paste0("0*", pgs21),
            paste0("0*", pgs22),
            paste0("0*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


# Null model (separate for EXT)
# Cognitive
formulas_null_EXTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      paste(
        outcome, "~",
        paste(
          c(
            paste0("0*", pgs11),
            paste0("0*", pgs12),
            paste0("0*", pgs13),
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)
# Non-cognitive
formulas_null_EXTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            paste0("0*", pgs21),
            paste0("0*", pgs22),
            paste0("0*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Null model (separate for INT)
# Cognitive
formulas_null_INTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      paste(
        outcome, "~",
        paste(
          c(
            paste0("0*", pgs11),
            paste0("0*", pgs12),
            paste0("0*", pgs13),
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Non-cognitive
formulas_null_INTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            paste0("0*", pgs21),
            paste0("0*", pgs22),
            paste0("0*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


# ----------------------------------------------------------
# CONSTRAINED MODELS
# ----------------------------------------------------------

# Time invariance (separate for EXT)
#### Cognitive component ####
formulas_timeinv_EXTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_cog_C_", domain, "*", pgs11),
            paste0("b_cog_M_", domain, "*", pgs12),
            paste0("b_cog_F_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_M_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

#### Non-cognitive component ####
formulas_timeinv_EXTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_M_", domain, "*", pgs12),
            paste0("b_",lab,"cog_F_", domain, "*", pgs13),
            paste0("b_noncog_C_", domain, "*", pgs21),
            paste0("b_noncog_M_", domain, "*", pgs22),
            paste0("b_noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

#### Time invariance (separate for INT) ####
# Cognitive component
formulas_timeinv_INTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_cog_C_", domain, "*", pgs11),
            paste0("b_cog_M_", domain, "*", pgs12),
            paste0("b_cog_F_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_M_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Time invariance (separate for INT)
# Non-cognitive component
formulas_timeinv_INTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_M_", domain, "*", pgs12),
            paste0("b_",lab,"cog_F_", domain, "*", pgs13),
            paste0("b_noncog_C_", domain, "*", pgs21),
            paste0("b_noncog_M_", domain, "*", pgs22),
            paste0("b_noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


# Parent invariance (separate for EXT)
# Cognitive component
formulas_parentinv_EXTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_P_", domain, "*", pgs12),
            paste0("b_",lab,"cog_P_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_M_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Parent invariance (separate for EXT)
# Non-cognitive component
formulas_parentinv_EXTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_M_", domain, "*", pgs12),
            paste0("b_",lab,"cog_F_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_P_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Parent invariance (separate for INT)
# Cognitive component
formulas_parentinv_INTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_P_", domain, "*", pgs12),
            paste0("b_",lab,"cog_P_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_M_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Parent invariance (separate for INT)
# Non-cognitive component
formulas_parentinv_INTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_M_", domain, "*", pgs12),
            paste0("b_",lab,"cog_F_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_P_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


#### Time + Parent invariance (separate for EXT) ####
# Cognitive component
formulas_timeparentinv_EXTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_cog_C_", domain, "*", pgs11),
            paste0("b_cog_P_", domain, "*", pgs12),
            paste0("b_cog_P_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_M_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Non-Cognitive component
formulas_timeparentinv_EXTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_M_", domain, "*", pgs12),
            paste0("b_",lab,"cog_F_", domain, "*", pgs13),
            paste0("b_noncog_C_", domain, "*", pgs21),
            paste0("b_noncog_P_", domain, "*", pgs22),
            paste0("b_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


#### Time + Parent invariance (separate for INT) ####
# Cognitive component
formulas_timeparentinv_INTc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_cog_C_", domain, "*", pgs11),
            paste0("b_cog_P_", domain, "*", pgs12),
            paste0("b_cog_P_", domain, "*", pgs13),
            paste0("b_",lab,"noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"noncog_M_", domain, "*", pgs22),
            paste0("b_",lab,"noncog_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Non-Cognitive component
formulas_timeparentinv_INTnc <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"cog_M_", domain, "*", pgs12),
            paste0("b_",lab,"cog_F_", domain, "*", pgs13),
            paste0("b_noncog_C_", domain, "*", pgs21),
            paste0("b_noncog_P_", domain, "*", pgs22),
            paste0("b_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


# Component invariance (separate for EXT)
formulas_compinv_EXT <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("EXT", outcome)) {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"_C_", domain, "*", pgs11),
            paste0("b_",lab,"_M_", domain, "*", pgs12),
            paste0("b_",lab,"_F_", domain, "*", pgs13),
            paste0("b_",lab,"_C_", domain, "*", pgs21),
            paste0("b_",lab,"_M_", domain, "*", pgs22),
            paste0("b_",lab,"_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)

# Component invariance (separate for INT)
formulas_compinv_INT <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"_C_", domain, "*", pgs11),
            paste0("b_",lab,"_M_", domain, "*", pgs12),
            paste0("b_",lab,"_F_", domain, "*", pgs13),
            paste0("b_",lab,"_C_", domain, "*", pgs21),
            paste0("b_",lab,"_M_", domain, "*", pgs22),
            paste0("b_",lab,"_F_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      paste(
        outcome, "~",
        paste(
          c(
            pgs11,
            pgs12,
            pgs13,
            pgs21,
            pgs22,
            pgs23
          ),
          collapse = " + "
        )
      )
    }
  }
)


# ----------------------------------------------------------
# BUILD MODELS
# ----------------------------------------------------------

# Null models
model_nullTidy_EXTc <- make_model(formulas_null_EXTc)
model_nullTidy_EXTnc <- make_model(formulas_null_EXTnc)
model_nullTidy_INTc <- make_model(formulas_null_INTc)
model_nullTidy_INTnc <- make_model(formulas_null_INTnc)
model_nullTidy_EXTj <- make_model(formulas_null_EXTj)
model_nullTidy_INTj <- make_model(formulas_null_INTj)
# Main comparison
model_baseTidy <- make_model(formulas_base)
# Externalising
model_timeinvTidy_EXTc <- make_model(formulas_timeinv_EXTc)
model_timeinvTidy_EXTnc <- make_model(formulas_timeinv_EXTnc)
model_parentinvTidy_EXTc <- make_model(formulas_parentinv_EXTc)
model_parentinvTidy_EXTnc <- make_model(formulas_parentinv_EXTnc)
model_timeparentinvTidy_EXTc <- make_model(formulas_timeparentinv_EXTc)
model_timeparentinvTidy_EXTnc <- make_model(formulas_timeparentinv_EXTnc)
# Optional: constrain across cog-noncog
model_compinvTidy_EXT <- make_model(formulas_compinv_EXT)

# Internalising
model_timeinvTidy_INTc <- make_model(formulas_timeinv_INTc)
model_timeinvTidy_INTnc <- make_model(formulas_timeinv_INTnc)
model_parentinvTidy_INTc <- make_model(formulas_parentinv_INTc)
model_parentinvTidy_INTnc <- make_model(formulas_parentinv_INTnc)
model_timeparentinvTidy_INTc <- make_model(formulas_timeparentinv_INTc)
model_timeparentinvTidy_INTnc <- make_model(formulas_timeparentinv_INTnc)
# Optional: constrain across cog-noncog
model_compinvTidy_INT <- make_model(formulas_compinv_INT)

# ----------------------------------------------------------
# FIT MODELS
# ----------------------------------------------------------

# 0
fit_null_EXTc <- sem(model_nullTidy_EXTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_null_svy_EXTc <- lavaan.survey(lavaan.fit = fit_null_EXTc, survey.design = svy_design, estimator = "MLM")
fit_null_EXTnc <- sem(model_nullTidy_EXTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_null_svy_EXTnc <- lavaan.survey(lavaan.fit = fit_null_EXTnc, survey.design = svy_design, estimator = "MLM")
fit_null_EXTj <- sem(model_nullTidy_EXTj, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_null_svy_EXTj <- lavaan.survey(lavaan.fit = fit_null_EXTj, survey.design = svy_design, estimator = "MLM")


fit_null_INTc <- sem(model_nullTidy_INTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_null_svy_INTc <- lavaan.survey(lavaan.fit = fit_null_INTc, survey.design = svy_design, estimator = "MLM")
fit_null_INTnc <- sem(model_nullTidy_INTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_null_svy_INTnc <- lavaan.survey(lavaan.fit = fit_null_INTnc, survey.design = svy_design, estimator = "MLM")
fit_null_INTj <- sem(model_nullTidy_INTj, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_null_svy_INTj <- lavaan.survey(lavaan.fit = fit_null_INTj, survey.design = svy_design, estimator = "MLM")


# 1
fit_base <- sem(model_baseTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_base_svy <- lavaan.survey(lavaan.fit = fit_base, survey.design = svy_design, estimator = "MLM")

# 2
fit_timeinv_EXTc <- sem(model_timeinvTidy_EXTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeinv_svy_EXTc <- lavaan.survey(lavaan.fit = fit_timeinv_EXTc, survey.design = svy_design, estimator = "MLM")
fit_timeinv_EXTnc <- sem(model_timeinvTidy_EXTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeinv_svy_EXTnc <- lavaan.survey(lavaan.fit = fit_timeinv_EXTnc, survey.design = svy_design, estimator = "MLM")
# 2
fit_timeinv_INTc <- sem(model_timeinvTidy_INTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeinv_svy_INTc <- lavaan.survey(lavaan.fit = fit_timeinv_INTc, survey.design = svy_design, estimator = "MLM")
fit_timeinv_INTnc <- sem(model_timeinvTidy_INTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeinv_svy_INTnc <- lavaan.survey(lavaan.fit = fit_timeinv_INTnc, survey.design = svy_design, estimator = "MLM")

# 3
fit_parentinv_EXTc  <- sem(model_parentinvTidy_EXTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_parentinv_svy_EXTc <- lavaan.survey(lavaan.fit = fit_parentinv_EXTc , survey.design = svy_design, estimator = "MLM")
fit_parentinv_EXTnc  <- sem(model_parentinvTidy_EXTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_parentinv_svy_EXTnc <- lavaan.survey(lavaan.fit = fit_parentinv_EXTnc , survey.design = svy_design, estimator = "MLM")
# 3
fit_parentinv_INTc  <- sem(model_parentinvTidy_INTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_parentinv_svy_INTc <- lavaan.survey(lavaan.fit = fit_parentinv_INTc , survey.design = svy_design, estimator = "MLM")
fit_parentinv_INTnc  <- sem(model_parentinvTidy_INTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_parentinv_svy_INTnc <- lavaan.survey(lavaan.fit = fit_parentinv_INTnc , survey.design = svy_design, estimator = "MLM")

# 4
fit_timeparentinv_EXTc <- sem(model_timeparentinvTidy_EXTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeparentinv_svy_EXTc <- lavaan.survey(lavaan.fit = fit_timeparentinv_EXTc, survey.design = svy_design, estimator = "MLM")
fit_timeparentinv_EXTnc <- sem(model_timeparentinvTidy_EXTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeparentinv_svy_EXTnc <- lavaan.survey(lavaan.fit = fit_timeparentinv_EXTnc, survey.design = svy_design, estimator = "MLM")
# 4
fit_timeparentinv_INTc <- sem(model_timeparentinvTidy_INTc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeparentinv_svy_INTc <- lavaan.survey(lavaan.fit = fit_timeparentinv_INTc, survey.design = svy_design, estimator = "MLM")
fit_timeparentinv_INTnc <- sem(model_timeparentinvTidy_INTnc, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeparentinv_svy_INTnc <- lavaan.survey(lavaan.fit = fit_timeparentinv_INTnc, survey.design = svy_design, estimator = "MLM")

# 5
fit_compinv_EXT <- sem(model_compinvTidy_EXT, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_compinv_svy_EXT <- lavaan.survey(lavaan.fit = fit_compinv_EXT, survey.design = svy_design, estimator = "MLM")
fit_compinv_INT <- sem(model_compinvTidy_INT, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_compinv_svy_INT <- lavaan.survey(lavaan.fit = fit_compinv_INT, survey.design = svy_design, estimator = "MLM")

#########################################
# MODEL COMPARISONS
#########################################

# Externalising
# Null
comp_null_EXTc <- anova(fit_base_svy, fit_null_svy_EXTc)
comp_null_EXTnc <- anova(fit_base_svy, fit_null_svy_EXTnc)
comp_null_EXTj <- anova(fit_base_svy, fit_null_svy_EXTj)
# Main comparisons
# Cognitive
comp_timeinv_EXTc <- anova(fit_base_svy, fit_timeinv_svy_EXTc)
comp_parentinv_EXTc <- anova(fit_base_svy, fit_parentinv_svy_EXTc)
comp_timeparentinv_EXTc <- anova(fit_parentinv_svy_EXTc, fit_timeparentinv_svy_EXTc)
# Non-cognitive
comp_timeinv_EXTnc <- anova(fit_base_svy, fit_timeinv_svy_EXTnc)
comp_parentinv_EXTnc <- anova(fit_base_svy, fit_parentinv_svy_EXTnc)
comp_timeparentinv_EXTnc <- anova(fit_parentinv_svy_EXTnc, fit_timeparentinv_svy_EXTnc)
# Component constrain
comp_compinv_EXT <- anova(fit_base_svy, fit_compinv_svy_EXT)

# Null
(Chisq_00a <- extract_chisq(comp_null_EXTc, modelname = 'CNC_null_EXTc'))
(Chisq_01a <- extract_chisq(comp_null_EXTnc, modelname = 'CNC_null_EXTnc'))
(Chisq_02a <- extract_chisq(comp_null_EXTj, modelname = 'CNC_null_EXTj'))
# Main comparisons
# Cognitive
(Chisq_1a <- extract_chisq(comp_timeinv_EXTc, modelname = 'CNC_timeinv_EXTc'))
(Chisq_2a <- extract_chisq(comp_parentinv_EXTc, modelname = 'CNC_parentinv_EXTc'))
(Chisq_3a <- extract_chisq(comp_timeparentinv_EXTc, modelname = 'CNC_timeparentinv_EXTc'))
# Non-cognitive
(Chisq_4a <- extract_chisq(comp_timeinv_EXTnc, modelname = 'CNC_timeinv_EXTnc'))
(Chisq_5a <- extract_chisq(comp_parentinv_EXTnc, modelname = 'CNC_parentinv_EXTnc'))
(Chisq_6a <- extract_chisq(comp_timeparentinv_EXTnc, modelname = 'CNC_timeparentinv_EXTnc'))
# Component constrain
(Chisq_7a <- extract_chisq(comp_compinv_EXT, modelname = 'CNC_compinv_EXTc'))

# Internalising
# Null models
comp_null_INTc <- anova(fit_base_svy, fit_null_svy_INTc)
comp_null_INTnc <- anova(fit_base_svy, fit_null_svy_INTnc)
comp_null_INTj <- anova(fit_base_svy, fit_null_svy_INTj)
# Main comparisons
# Cognitive
comp_timeinv_INTc <- anova(fit_base_svy, fit_timeinv_svy_INTc)
comp_parentinv_INTc <- anova(fit_base_svy, fit_parentinv_svy_INTc)
comp_timeparentinv_INTc <- anova(fit_parentinv_svy_INTc, fit_timeparentinv_svy_INTc)
# Non-cognitive
comp_timeinv_INTnc <- anova(fit_base_svy, fit_timeinv_svy_INTnc)
comp_parentinv_INTnc <- anova(fit_base_svy, fit_parentinv_svy_INTnc)
comp_timeparentinv_INTnc <- anova(fit_parentinv_svy_INTnc, fit_timeparentinv_svy_INTnc)
# Component constrain
comp_compinv_INT <- anova(fit_base_svy, fit_compinv_svy_INT)

#### Run comparisons ####
# Null models
(Chisq_00b <- extract_chisq(comp_null_INTc, modelname = 'CNC_null_INTc'))
(Chisq_01b <- extract_chisq(comp_null_INTnc, modelname = 'CNC_null_INTnc'))
(Chisq_02b <- extract_chisq(comp_null_INTj, modelname = 'CNC_null_INTj'))
# Main comparisons
# Cognitive
(Chisq_1b <- extract_chisq(comp_timeinv_INTc, modelname = 'CNC_timeinv_INTc'))
(Chisq_2b <- extract_chisq(comp_parentinv_INTc, modelname = 'CNC_parentinv_INTc'))
(Chisq_3b <- extract_chisq(comp_timeparentinv_INTc, modelname = 'CNC_timeparentinv_INTc'))
# Non-cognitive
(Chisq_4b <- extract_chisq(comp_timeinv_INTnc, modelname = 'CNC_timeinv_INTnc'))
(Chisq_5b <- extract_chisq(comp_parentinv_INTnc, modelname = 'CNC_parentinv_INTnc'))
(Chisq_6b <- extract_chisq(comp_timeparentinv_INTnc, modelname = 'CNC_timeparentinv_INTnc'))
# Component constrain
(Chisq_7b <- extract_chisq(comp_compinv_INT, modelname = 'CNC_compinv_INT'))