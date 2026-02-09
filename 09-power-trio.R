#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   Power Calculations for SEM Models
# Purpose:  Perform Monte Carlo simulations to estimate 
#           statistical power for detecting direct and 
#           indirect genetic effects in trio-based 
#           structural equation models.
#
# Author:   Jose J. Morosoli
# Date:     09-02-2026
#############################################################

### --- Load Required Packages ---
library(tidyverse)
library(lavaan)
library(foreach)
library(parallel)
library(doParallel)
library(glue)

# ---- SESOI (standardized beta) ----
SESOI_beta <- 0.05  # half a "small effect" benchmark

### --- Define Population Model Generator ---
population.model <- function(ec, em, ef) {
  
  # implied covariances given:
  # Xc = 0.5*Xm + 0.5*Xf + e, Var(e)=0.5 => Var(Xc)=1
  cov_c_m <- 0.5
  cov_c_f <- 0.5
  cov_m_f <- 0.0
  
  # variance explained in Y by linear predictor b' Sigma b
  var_lp <- (ec^2)*1 + (em^2)*1 + (ef^2)*1 +
    2*ec*em*cov_c_m + 2*ec*ef*cov_c_f + 2*em*ef*cov_m_f
  
  resid_var <- 1 - var_lp
  if (resid_var <= 1e-6) resid_var <- 1e-6  # guard against non-positive variance
  
  glue('
    Xc ~ 0.5*Xm + 0.5*Xf
    Xc ~~ 0.5*Xc
    Xm ~~ 1*Xm
    Xf ~~ 1*Xf
    Xm ~~ 0*Xf
    Yc ~ {ec}*Xc + {em}*Xm + {ef}*Xf
    Yc ~~ {resid_var}*Yc
    gt := 0.5*{ec}
  ')
}

### --- Define Sample Model to Fit ---
samp.model <- "
  Yc ~ d*Xc + Xm + Xf
  gt := 0.5*d
"

### --- Define Effect Sizes (include SESOI = 0.10) ---
direct_effects   <- c(SESOI_beta, 0.08, 0.12, 0.16)
indirect_effects <- c(0.02, 0.04, SESOI_beta, 0.08)
param_grid <- expand.grid(ec = direct_effects, em = indirect_effects, ef = indirect_effects)

### --- Simulation Settings ---
N <- 3223
iterations <- 1000
numCores <- detectCores()
cl <- makeCluster(numCores - 1)
registerDoParallel(cl)
clusterEvalQ(cl, { 
  library(lavaan)
  library(tidyverse)
  library(glue)
})

### --- Run Simulations ---
results <- foreach(i = 1:nrow(param_grid), .combine = rbind) %dopar% {
  
  ec <- param_grid$ec[i]
  em <- param_grid$em[i]
  ef <- param_grid$ef[i]
  mod_text <- population.model(ec, em, ef)
  
  # Minimal, robust change: avoid nested foreach/%do% inside workers
  out_list <- vector("list", iterations)
  
  for (k in seq_len(iterations)) {
    dat <- simulateData(mod_text, model.type = "sem", sample.nobs = N)
    fit <- sem(samp.model, data = dat)
    est <- parameterEstimates(fit, standardized = TRUE)
    est$ec <- ec
    est$em <- em
    est$ef <- ef
    out_list[[k]] <- est
  }
  
  bind_rows(out_list)
}

stopCluster(cl)

### --- Summarize Power ---
power_summary <- results %>%
  filter(op == "~", lhs == "Yc", rhs %in% c("Xc", "Xm", "Xf")) %>%  # (optional) focus on key paths
  group_by(lhs, rhs, ec, em, ef) %>%
  summarise(power = mean(pvalue < 0.05), .groups = "drop")

### --- Print Power Table ---
print(power_summary)

#------------------------------------------------------------
#  Lakens (2022) Sensitivity analysis: smallest detectable effect
#  Given N, alpha, target power -> minimum detectable f^2
#------------------------------------------------------------

# If needed:
library(pwr)

alpha <- 0.05
target_power <- 0.80

# In the fitted model: Yc ~ Xc + Xm + Xf  (3 predictors)
p_predictors <- 3
u <- 1
v <- N - p_predictors - 1  # denominator df for testing one predictor

sens <- pwr.f2.test(u = u, v = v, sig.level = alpha, power = target_power)

f2_mde <- sens$f2
partialR2_mde <- f2_mde / (1 + f2_mde)

# Approximate mapping to standardized beta (benchmark)
beta_mde_approx <- sqrt(partialR2_mde)

cat("\n--- Sensitivity analysis (Lakens, 2022) ---\n")
cat("N =", N, " | alpha =", alpha, " | target power =", target_power, "\n")
cat("Minimum detectable f^2 =", round(f2_mde, 4), "\n")
cat("Minimum detectable partial R^2 =", round(partialR2_mde, 4), "\n")
cat("Approx. minimum detectable standardized beta =", round(beta_mde_approx, 3), "\n")

# --- Export to Excel ---
library(openxlsx)
out_file <- "power_results_trio_SEM.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "power_summary")
writeData(wb, "power_summary", power_summary)
addWorksheet(wb, "param_grid")
writeData(wb, "param_grid", param_grid)
saveWorkbook(wb, out_file, overwrite = TRUE)
message("Saved: ", out_file)
