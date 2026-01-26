#########################################
# Trio SEM Models: Educational Attainment
# Author: Jose J. Morosoli
# Date: 26-01-2026
# Purpose: Extract parameters estimates,
# obtain CIs and correct CIs using Bonferroni.
# Note: Run after 04-EA-analyses.R
#########################################


# Load required packages
library(writexl)

#### STEP 1 ####
# DEFINE SHARED OBJECTS
outcome1 <- "BEXT_res"
outcome2 <- "CEXT_res"
outcome3 <- "DEXT_res" 
outcome4 <- "FEXT_res" 
outcome5 <- "BINT_res" 
outcome6 <- "CINT_res" 
outcome7 <- "DINT_res" 
outcome8 <- "FINT_res"

## Modelling Educational Attainment
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

# ESTIMATE TIME-SPECIFIC MODEL #
# Same coefficients for parents for EXT and INT
# Same coefficients for time for INT
myformula_ea_time_1 <- paste(outcome1, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea_time_2 <- paste(outcome2, "~","be21*",pgs1,"+","be22*",pgs2,"+","be22*",pgs3, sep="")
myformula_ea_time_3 <- paste(outcome3, "~","be31*",pgs1,"+","be32*",pgs2,"+","be32*",pgs3, sep="")
myformula_ea_time_4 <- paste(outcome4, "~","be41*",pgs1,"+","be42*",pgs2,"+","be42*",pgs3, sep="")
myformula_ea_time_5 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea_time_6 <- paste(outcome6, "~","bi21*",pgs1,"+","bi22*",pgs2,"+","bi22*",pgs3, sep="")
myformula_ea_time_7 <- paste(outcome7, "~","bi31*",pgs1,"+","bi32*",pgs2,"+","bi32*",pgs3, sep="")
myformula_ea_time_8 <- paste(outcome8, "~","bi41*",pgs1,"+","bi42*",pgs2,"+","bi42*",pgs3, sep="")
#specify the free model using these formulas
model_ea_time_ <- paste(# regressions 
  myformula_ea_time_1,
  myformula_ea_time_2,
  myformula_ea_time_3,
  myformula_ea_time_4,
  myformula_ea_time_5,
  myformula_ea_time_6,
  myformula_ea_time_7,
  myformula_ea_time_8,
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
model_ea_time_Tidy <- noquote(strsplit(model_ea_time_, "\n")[[1]])

# Fit the free model
fit_ea_time <- sem(model_ea_time_Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
fit_ea_time_svy <- lavaan.survey(lavaan.fit = fit_ea_time, survey.design = svy_design, estimator = "MLM")


# ESTIMATE FULLY CONSTRAINED MODEL
# ESTIMATE TIME-SPECIFIC MODEL #
# Same coefficients for parents for EXT and INT
# Same coefficients for time for INT
myformula_ea_best_1 <- paste(outcome1, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea_best_2 <- paste(outcome2, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea_best_3 <- paste(outcome3, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea_best_4 <- paste(outcome4, "~","be11*",pgs1,"+","be12*",pgs2,"+","be12*",pgs3, sep="")
myformula_ea_best_5 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea_best_6 <- paste(outcome6, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea_best_7 <- paste(outcome7, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
myformula_ea_best_8 <- paste(outcome8, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3, sep="")
#specify the free model using these formulas
model_ea_best_ <- paste(# regressions 
  myformula_ea_best_1,
  myformula_ea_best_2,
  myformula_ea_best_3,
  myformula_ea_best_4,
  myformula_ea_best_5,
  myformula_ea_best_6,
  myformula_ea_best_7,
  myformula_ea_best_8,
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
model_ea_best_Tidy <- noquote(strsplit(model_ea_best_, "\n")[[1]])

# Fit the free model
fit_ea_best <- sem(model_ea_best_Tidy, data=aframe_svy_na, estimator = "MLR", missing = "FIML")  
# Apply the survey design to the lavaan model 
fit_ea_best_svy <- lavaan.survey(lavaan.fit = fit_ea_best, survey.design = svy_design, estimator = "MLM")


#### STEP 2 ####
# Get standardized estimates with confidence intervals
estimates_cistd_base <- standardizedSolution(fit_ea_best_svy)
estimates_cistd_time <- standardizedSolution(fit_ea_time_svy ) 
estimates_cistd_best <- standardizedSolution(fit_ea_best_svy ) 

# Select and rename key columns
table_data_base <- estimates_cistd_base[, c("lhs","op","rhs", "est.std", "ci.lower", "ci.upper","se","pvalue")]
table_data_time <- estimates_cistd_time[, c("lhs","op","rhs", "est.std", "ci.lower", "ci.upper","se","pvalue")]
table_data_best <- estimates_cistd_best[, c("lhs","op","rhs", "est.std", "ci.lower", "ci.upper","se","pvalue")]
table_data_base$p_fdr <- p.adjust(table_data_base$pvalue, method="fdr")
table_data_time$p_fdr <- p.adjust(table_data_time$pvalue, method="fdr")
table_data_best$p_fdr <- p.adjust(table_data_best$pvalue, method="fdr")
names(table_data_base) <- c("outcome", "operator","predictor", "beta.std", "ci.lower","ci.upper","se","pvalue","qvalue")
names(table_data_time) <- c("outcome", "operator","predictor", "beta.std", "ci.lower","ci.upper","se","pvalue","qvalue")
names(table_data_best) <- c("outcome", "operator","predictor", "beta.std", "ci.lower","ci.upper","se","pvalue","qvalue")

#### STEP 3 ####
# Get multiple-testing-adjusted CIs (Bonferroni)
alpha_family <- 0.05
m_tests <- 16                       # based on our Methods
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
write_xlsx(table_data_base, path = "fit_ea_base_svy_R2.xlsx") # base
write_xlsx(table_data_time, path = "fit_ea_time_svy_R2.xlsx") # time-specific only
write_xlsx(table_data_best, path = "fit_ea_best_svy_R2.xlsx") # best-fitting model
