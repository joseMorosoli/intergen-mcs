#########################################
# Trio SEM Models: Cognitive and Non-cognitive PGS
# Author: Jose J. Morosoli
# Date: 26-01-2026
# Purpose: Extract parameters estimates,
# obtain CIs and correct CIs using Bonferroni.
# Note: Run after 04-EA-analyses.R
#########################################

# Load required packages
library(writexl)

#### STEP 1 ####

# ESTIMATE TIME-SPECIFIC MODEL #
formulas_timeCNC <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"_cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"_cog_P_", domain, "*", pgs12),
            paste0("b_",lab,"_cog_P_", domain, "*", pgs13),
            paste0("b_",lab,"_noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"_noncog_P_", domain, "*", pgs22),
            paste0("b_",lab,"_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_",lab,"_cog_C_", domain, "*", pgs11),
            paste0("b_",lab,"_cog_P_", domain, "*", pgs12),
            paste0("b_",lab,"_cog_P_", domain, "*", pgs13),
            paste0("b_",lab,"_noncog_C_", domain, "*", pgs21),
            paste0("b_",lab,"_noncog_P_", domain, "*", pgs22),
            paste0("b_",lab,"_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    }
  }
)


# ----------------------------------------------------------
# BUILD MODEL
# ----------------------------------------------------------

model_timeCNCTidy <- make_model(formulas_timeCNC)

# ----------------------------------------------------------
# FIT MODEL
# ----------------------------------------------------------

fit_timeCNC <- sem(model_timeCNCTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_timeCNC_svy <- lavaan.survey(lavaan.fit = fit_timeCNC, survey.design = svy_design, estimator = "MLM")

# Get standardized estimates with confidence intervals
estimates_cistd <- standardizedSolution(fit_base) # fit_base fit_timeCNC_svy
estimates <- parameterEstimates(fit_timeCNC_svy , standardized = TRUE, ci = TRUE, level = 0.95)
summary(fit_timeCNC_svy, standardize = T, rsquare = T)

# Select and rename key columns
table_data <- estimates[, c("lhs","op","rhs", "std.all", "ci.lower", "ci.upper","est","se")]
#names(table_data) <- c("outcome", "operator","predictor", "beta.std", "ci.lower","ci.upper","est","se")
table_data <- estimates_cistd[, c("lhs","op","rhs", "est.std", "ci.lower", "ci.upper","se")]
names(table_data) <- c("outcome", "operator","predictor", "beta.std", "ci.lower","ci.upper","se")

# Load required packages
library(writexl)

# Export to Excel
write_xlsx(table_data, path = "time_CNC_svy.xlsx")



# ESTIMATE time FITTING MODEL #
model_COGNON_time <- '
# Same coefficients for parents for EXT and INT
# Same coefficients for time for INT
BEXT_res ~ b11c*EA_cog_C_res+b11p*EA_cog_M_res+b11p*EA_cog_F_res+b21c*EA_noncog_C_res+b21p*EA_noncog_M_res+b21p*EA_noncog_F_res
CEXT_res ~ b12c*EA_cog_C_res+b12p*EA_cog_M_res+b12p*EA_cog_F_res+b22c*EA_noncog_C_res+b22p*EA_noncog_M_res+b22p*EA_noncog_F_res
DEXT_res ~ b13c*EA_cog_C_res+b13p*EA_cog_M_res+b13p*EA_cog_F_res+b23c*EA_noncog_C_res+b23p*EA_noncog_M_res+b23p*EA_noncog_F_res
FEXT_res ~ b14c*EA_cog_C_res+b14p*EA_cog_M_res+b14p*EA_cog_F_res+b24c*EA_noncog_C_res+b24p*EA_noncog_M_res+b24p*EA_noncog_F_res

BINT_res ~ b31c*EA_cog_C_res+b31p*EA_cog_M_res+b31p*EA_cog_F_res+b41c*EA_noncog_C_res+b41p*EA_noncog_M_res+b41p*EA_noncog_F_res
CINT_res ~ b32c*EA_cog_C_res+b32p*EA_cog_M_res+b31p*EA_cog_F_res+b42c*EA_noncog_C_res+b42p*EA_noncog_M_res+b42p*EA_noncog_F_res
DINT_res ~ b33c*EA_cog_C_res+b33p*EA_cog_M_res+b31p*EA_cog_F_res+b43c*EA_noncog_C_res+b43p*EA_noncog_M_res+b43p*EA_noncog_F_res
FINT_res ~ b34c*EA_cog_C_res+b34p*EA_cog_M_res+b31p*EA_cog_F_res+b44c*EA_noncog_C_res+b44p*EA_noncog_M_res+b44p*EA_noncog_F_res
              
# Correlated residuals among outcomes
BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res
CEXT_res ~~ DEXT_res + FEXT_res
DEXT_res ~~ FEXT_res

# Cross-domain covariances
BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res
CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res
DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res
FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res

# Internalizing covariances (if needed)
BINT_res ~~ CINT_res + DINT_res + FINT_res
CINT_res ~~ DINT_res + FINT_res
DINT_res ~~ FINT_res

# Parental genotype correlations
EA_cog_M_res ~~ EA_cog_F_res
EA_noncog_M_res ~~ EA_noncog_F_res
'

fit_COGNON_time <- sem(model_COGNON_time, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_COGNON_time_svy <- lavaan.survey(fit_COGNON_time, survey.design=svy_design, estimator="MLM")




# ESTIMATE BEST-FITTING MODEL #
formulas_bestCNC <- sapply(
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
            paste0("b_noncog_C_", domain, "*", pgs21),
            paste0("b_noncog_P_", domain, "*", pgs22),
            paste0("b_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      domain <- "EXT"
      lab = substr(outcome,1,1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_cog_C_", domain, "*", pgs11),
            paste0("b_cog_P_", domain, "*", pgs12),
            paste0("b_cog_P_", domain, "*", pgs13),
            paste0("b_noncog_C_", domain, "*", pgs21),
            paste0("b_noncog_P_", domain, "*", pgs22),
            paste0("b_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    }
  }
)


# ----------------------------------------------------------
# BUILD MODEL
# ----------------------------------------------------------

model_bestCNCTidy <- make_model(formulas_bestCNC)

# ----------------------------------------------------------
# FIT MODEL
# ----------------------------------------------------------

fit_bestCNC <- sem(model_bestCNCTidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_bestCNC_svy <- lavaan.survey(lavaan.fit = fit_bestCNC, survey.design = svy_design, estimator = "MLM")



# ESTIMATE BEST FITTING MODEL #
model_COGNON_best <- '
# Same coefficients for parents for EXT and INT
# Same coefficients for time for INT
BEXT_res ~ b11c*EA_cog_C_res+b11p*EA_cog_M_res+b11p*EA_cog_F_res+b21c*EA_noncog_C_res+b21p*EA_noncog_M_res+b21p*EA_noncog_F_res
CEXT_res ~ b11c*EA_cog_C_res+b11p*EA_cog_M_res+b11p*EA_cog_F_res+b21c*EA_noncog_C_res+b21p*EA_noncog_M_res+b21p*EA_noncog_F_res
DEXT_res ~ b11c*EA_cog_C_res+b11p*EA_cog_M_res+b11p*EA_cog_F_res+b21c*EA_noncog_C_res+b21p*EA_noncog_M_res+b21p*EA_noncog_F_res
FEXT_res ~ b11c*EA_cog_C_res+b11p*EA_cog_M_res+b11p*EA_cog_F_res+b21c*EA_noncog_C_res+b21p*EA_noncog_M_res+b21p*EA_noncog_F_res

BINT_res ~ b31c*EA_cog_C_res+b31p*EA_cog_M_res+b31p*EA_cog_F_res+b41c*EA_noncog_C_res+b41p*EA_noncog_M_res+b41p*EA_noncog_F_res
CINT_res ~ b31c*EA_cog_C_res+b31p*EA_cog_M_res+b31p*EA_cog_F_res+b41c*EA_noncog_C_res+b41p*EA_noncog_M_res+b41p*EA_noncog_F_res
DINT_res ~ b31c*EA_cog_C_res+b31p*EA_cog_M_res+b31p*EA_cog_F_res+b41c*EA_noncog_C_res+b41p*EA_noncog_M_res+b41p*EA_noncog_F_res
FINT_res ~ b31c*EA_cog_C_res+b31p*EA_cog_M_res+b31p*EA_cog_F_res+b41c*EA_noncog_C_res+b41p*EA_noncog_M_res+b41p*EA_noncog_F_res
              
# Correlated residuals among outcomes
BEXT_res ~~ CEXT_res + DEXT_res + FEXT_res
CEXT_res ~~ DEXT_res + FEXT_res
DEXT_res ~~ FEXT_res

# Cross-domain covariances
BINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res
CINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res
DINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res
FINT_res ~~ BEXT_res + CEXT_res + DEXT_res + FEXT_res

# Internalizing covariances (if needed)
BINT_res ~~ CINT_res + DINT_res + FINT_res
CINT_res ~~ DINT_res + FINT_res
DINT_res ~~ FINT_res

# Parental genotype correlations
EA_cog_M_res ~~ EA_cog_F_res
EA_noncog_M_res ~~ EA_noncog_F_res
'

fit_COGNON_best <- sem(model_COGNON_best, data=aframe_svy_na, estimator = "MLR", missing = "FIML")
fit_COGNON_best_svy <- lavaan.survey(fit_COGNON_best, survey.design=svy_design, estimator="MLM")



#### STEP 2 ####
# Get standardized estimates with confidence intervals
estimates_cistd_base <- standardizedSolution(fit_base_svy )
estimates_cistd_time <- standardizedSolution(fit_COGNON_time_svy )
estimates_cistd_best <- standardizedSolution(fit_COGNON_best_svy )
# Get unstandardized estimates with confidence intervals
estimates_unstd_base <- parameterEstimates(fit_base_svy, standardized = TRUE, ci = TRUE)
estimates_unstd_time <- parameterEstimates(fit_COGNON_time_svy, standardized = TRUE, ci = TRUE)
estimates_unstd_best <- parameterEstimates(fit_COGNON_best_svy, standardized = TRUE, ci = TRUE)
# Select and rename key columns
table_data_base <- cbind(estimates_cistd_base[, c("lhs","op","rhs")],
                         estimates_unstd_base[, c("est","se")],
                         estimates_cistd_base[, c("est.std", "ci.lower", "ci.upper","se","pvalue")])
table_data_time <- cbind(estimates_cistd_time[, c("lhs","op","rhs")],
                         estimates_unstd_time[, c("est","se")],
                         estimates_cistd_time[, c("est.std", "ci.lower", "ci.upper","se","pvalue")])
table_data_best <- cbind(estimates_cistd_best[, c("lhs","op","rhs")],
                                              estimates_unstd_best[, c("est","se")],
                                              estimates_cistd_best[, c("est.std", "ci.lower", "ci.upper","se","pvalue")])
table_data_base$p_fdr <- p.adjust(table_data_base$pvalue, method="fdr")
table_data_time$p_fdr <- p.adjust(table_data_time$pvalue, method="fdr")
table_data_best$p_fdr <- p.adjust(table_data_best$pvalue, method="fdr")
names(table_data_base) <- c("outcome", "operator","predictor", 'beta.obs','se.obs',"beta.std", "ci.lower","ci.upper","se","pvalue","qvalue")
names(table_data_time) <- c("outcome", "operator","predictor", 'beta.obs','se.obs',"beta.std", "ci.lower","ci.upper","se","pvalue","qvalue")
names(table_data_best) <- c("outcome", "operator","predictor", 'beta.obs','se.obs',"beta.std", "ci.lower","ci.upper","se","pvalue","qvalue")

#### STEP 3 ####
# Get multiple-testing-adjusted CIs (Bonferroni)
alpha_family <- 0.05
m_tests <- 32                       # based on our Methods
alpha_adj <- alpha_family / m_tests
ci_level_adj <- 1 - alpha_adj
z_adj <- qnorm(1 - alpha_adj/2)     # two-sided
# Base model
table_data_base$ci.lower.adj = table_data_base$beta.std - z_adj * table_data_base$se
table_data_base$ci.upper.adj = table_data_base$beta.std + z_adj * table_data_base$se
table_data_base$ci_level_adj = ci_level_adj
table_data_base$alpha_adj = alpha_adj
table_data_base$m_tests = m_tests
# Time-specific model
table_data_time$ci.lower.adj = table_data_time$beta.std - z_adj * table_data_time$se
table_data_time$ci.upper.adj = table_data_time$beta.std + z_adj * table_data_time$se
table_data_time$ci_level_adj = ci_level_adj
table_data_time$alpha_adj = alpha_adj
table_data_time$m_tests = m_tests
# Fully-constrained model
table_data_best$ci.lower.adj = table_data_best$beta.std - z_adj * table_data_best$se
table_data_best$ci.upper.adj = table_data_best$beta.std + z_adj * table_data_best$se
table_data_best$ci_level_adj = ci_level_adj
table_data_best$alpha_adj = alpha_adj
table_data_best$m_tests = m_tests

#### STEP 4 ####
# Export to Excel
write_xlsx(table_data_base, path = "fit_cnc_base_svy_R2.xlsx") # base
write_xlsx(table_data_time, path = "fit_cnc_time_svy_R2.xlsx") # time-specific only
write_xlsx(table_data_best, path = "fit_cnc_best_svy_R2.xlsx") # best-fitting model
