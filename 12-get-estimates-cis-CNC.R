#########################################
# Trio SEM Models: Cognitive and Non-cognitive PGS
# Author: Jose J. Morosoli
# Date: 26-01-2026
# Purpose: Extract parameter estimates,
# obtain CIs and multiplicity-adjusted CIs.
# Note: Run after 04-EA-analyses.R
#########################################

# Load required packages
library(writexl)

#### STEP 1 ####

# ----------------------------------------------------------
# PARENT-INVARIANT, TIME-SPECIFIC MODEL
# Same parental coefficients within each age.
# Coefficients remain free to vary across ages for EXT and INT.
# ----------------------------------------------------------

formulas_timeCNC <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
      lab <- substr(outcome, 1, 1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_", lab, "_cog_C_", domain, "*", pgs11),
            paste0("b_", lab, "_cog_P_", domain, "*", pgs12),
            paste0("b_", lab, "_cog_P_", domain, "*", pgs13),
            paste0("b_", lab, "_noncog_C_", domain, "*", pgs21),
            paste0("b_", lab, "_noncog_P_", domain, "*", pgs22),
            paste0("b_", lab, "_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    } else {
      domain <- "EXT"
      lab <- substr(outcome, 1, 1)
      paste(
        outcome, "~",
        paste(
          c(
            paste0("b_", lab, "_cog_C_", domain, "*", pgs11),
            paste0("b_", lab, "_cog_P_", domain, "*", pgs12),
            paste0("b_", lab, "_cog_P_", domain, "*", pgs13),
            paste0("b_", lab, "_noncog_C_", domain, "*", pgs21),
            paste0("b_", lab, "_noncog_P_", domain, "*", pgs22),
            paste0("b_", lab, "_noncog_P_", domain, "*", pgs23)
          ),
          collapse = " + "
        )
      )
    }
  }
)

model_timeCNCTidy <- make_model(formulas_timeCNC)

fit_timeCNC <- sem(
  model_timeCNCTidy,
  data = aframe_svy_na,
  estimator = "MLR",
  missing = "FIML"
)

fit_timeCNC_svy <- lavaan.survey(
  lavaan.fit = fit_timeCNC,
  survey.design = svy_design,
  estimator = "MLM"
)


# ----------------------------------------------------------
# BEST-FITTING FULLY CONSTRAINED MODEL
# Same parental coefficients across parents and ages.
# Child coefficients are also constrained equal across ages.
# Cog and NonCog remain separate, as do EXT and INT.
# ----------------------------------------------------------

formulas_bestCNC <- sapply(
  list(outcome1, outcome2, outcome3, outcome4, outcome5, outcome6, outcome7, outcome8),
  function(outcome) {
    if (grepl("INT", outcome)) {
      domain <- "INT"
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

model_bestCNCTidy <- make_model(formulas_bestCNC)

fit_bestCNC <- sem(
  model_bestCNCTidy,
  data = aframe_svy_na,
  estimator = "MLR",
  missing = "FIML"
)

fit_bestCNC_svy <- lavaan.survey(
  lavaan.fit = fit_bestCNC,
  survey.design = svy_design,
  estimator = "MLM"
)


#### STEP 2 ####
# Extract estimates for:
# 1. Base model: fit_base_svy (created upstream)
# 2. Parent-invariant, time-specific model: fit_timeCNC_svy
# 3. Best-fitting fully constrained model: fit_bestCNC_svy

# Standardized estimates with confidence intervals
estimates_cistd_base <- standardizedSolution(fit_base_svy)
estimates_cistd_time <- standardizedSolution(fit_timeCNC_svy)
estimates_cistd_best <- standardizedSolution(fit_bestCNC_svy)

# Unstandardized estimates
estimates_unstd_base <- parameterEstimates(
  fit_base_svy,
  standardized = TRUE,
  ci = TRUE
)

estimates_unstd_time <- parameterEstimates(
  fit_timeCNC_svy,
  standardized = TRUE,
  ci = TRUE
)

estimates_unstd_best <- parameterEstimates(
  fit_bestCNC_svy,
  standardized = TRUE,
  ci = TRUE
)

# Combine key columns
table_data_base <- cbind(
  estimates_cistd_base[, c("lhs", "op", "rhs")],
  estimates_unstd_base[, c("est", "se")],
  estimates_cistd_base[, c("est.std", "ci.lower", "ci.upper", "se", "pvalue")]
)

table_data_time <- cbind(
  estimates_cistd_time[, c("lhs", "op", "rhs")],
  estimates_unstd_time[, c("est", "se")],
  estimates_cistd_time[, c("est.std", "ci.lower", "ci.upper", "se", "pvalue")]
)

table_data_best <- cbind(
  estimates_cistd_best[, c("lhs", "op", "rhs")],
  estimates_unstd_best[, c("est", "se")],
  estimates_cistd_best[, c("est.std", "ci.lower", "ci.upper", "se", "pvalue")]
)

# FDR-adjusted parameter P values retained for Supplementary Data output
table_data_base$p_fdr <- p.adjust(table_data_base$pvalue, method = "fdr")
table_data_time$p_fdr <- p.adjust(table_data_time$pvalue, method = "fdr")
table_data_best$p_fdr <- p.adjust(table_data_best$pvalue, method = "fdr")

names(table_data_base) <- c(
  "outcome", "operator", "predictor", "beta.obs", "se.obs",
  "beta.std", "ci.lower", "ci.upper", "se", "pvalue", "qvalue"
)

names(table_data_time) <- c(
  "outcome", "operator", "predictor", "beta.obs", "se.obs",
  "beta.std", "ci.lower", "ci.upper", "se", "pvalue", "qvalue"
)

names(table_data_best) <- c(
  "outcome", "operator", "predictor", "beta.obs", "se.obs",
  "beta.std", "ci.lower", "ci.upper", "se", "pvalue", "qvalue"
)


#### STEP 3 ####
# Get multiple-testing-adjusted CIs (Bonferroni)

alpha_family <- 0.05
m_tests <- 32
alpha_adj <- alpha_family / m_tests
ci_level_adj <- 1 - alpha_adj
z_adj <- qnorm(1 - alpha_adj / 2)

# Base model
table_data_base$ci.lower.adj <- table_data_base$beta.std - z_adj * table_data_base$se
table_data_base$ci.upper.adj <- table_data_base$beta.std + z_adj * table_data_base$se
table_data_base$ci_level_adj <- ci_level_adj
table_data_base$alpha_adj <- alpha_adj
table_data_base$m_tests <- m_tests

# Time-specific model
table_data_time$ci.lower.adj <- table_data_time$beta.std - z_adj * table_data_time$se
table_data_time$ci.upper.adj <- table_data_time$beta.std + z_adj * table_data_time$se
table_data_time$ci_level_adj <- ci_level_adj
table_data_time$alpha_adj <- alpha_adj
table_data_time$m_tests <- m_tests

# Best-fitting model
table_data_best$ci.lower.adj <- table_data_best$beta.std - z_adj * table_data_best$se
table_data_best$ci.upper.adj <- table_data_best$beta.std + z_adj * table_data_best$se
table_data_best$ci_level_adj <- ci_level_adj
table_data_best$alpha_adj <- alpha_adj
table_data_best$m_tests <- m_tests


#### STEP 4 ####
# Export to Excel

write_xlsx(
  table_data_base,
  path = "fit_cnc_base_svy_R3.xlsx"
)

write_xlsx(
  table_data_time,
  path = "fit_cnc_time_svy_R3.xlsx"
)

write_xlsx(
  table_data_best,
  path = "fit_cnc_best_svy_R3.xlsx"
)
