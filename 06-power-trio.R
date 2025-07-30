# Load required packages
library(tidyverse)
library(lavaan)
library(foreach)
library(parallel)
library(doParallel)
library(glue)

# Define model generator
population.model <- function(ec, em, ef) {
  glue('
    Xc ~ 0.5*Xm + 0.5*Xf
    Xc ~~ 0.5*Xc
    Xm ~~ 1*Xm
    Xf ~~ 1*Xf
    Xm ~~ 0*Xf
    Yc ~ {ec}*Xc + {em}*Xm + {ef}*Xf
    Yc ~~ (1 - ({ec}^2 + {em}^2 + {ef}^2 + {ec}*{em} + {ec}*{ef}))*Yc
    gt := 0.5*{ec}
  ')
}

# Model to fit
samp.model <- "
  Yc ~ d*Xc + Xm + Xf
  gt := 0.5*d
"

# Define effect sizes
direct_effects <- c(0.08, 0.12, 0.16)
indirect_effects <- c(0.02, 0.04, 0.08)
param_grid <- expand.grid(ec = direct_effects, em = indirect_effects, ef = indirect_effects)

# Simulation settings
N <- 3223
iterations <- 1000
numCores <- detectCores()
cl <- makeCluster(numCores - 1)
registerDoParallel(cl)
clusterEvalQ(cl, { library(lavaan); library(tidyverse); library(glue); library(tidyverse); library(foreach)})

# Run simulation
results <- foreach(i = 1:nrow(param_grid), .combine = rbind) %dopar% {
  ec <- param_grid$ec[i]
  em <- param_grid$em[i]
  ef <- param_grid$ef[i]
  mod_text <- population.model(ec, em, ef)
  
  res <- foreach(k = 1:iterations, .combine = rbind) %do% {
    dat <- simulateData(mod_text, model.type = "sem", sample.nobs = N)
    fit <- sem(samp.model, data = dat)
    est <- parameterEstimates(fit, standardized = TRUE)
    est$ec <- ec
    est$em <- em
    est$ef <- ef
    est
  }
  res
}

stopCluster(cl)

# Summarize power
power_summary <- results %>%
  filter(op == "~") %>%
  group_by(lhs, rhs, ec, em, ef) %>%
  summarise(power = mean(pvalue < 0.05), .groups = "drop")

# Print power table
print(power_summary)
