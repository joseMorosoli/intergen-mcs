#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   Plotting PGS Effects (Direct & Indirect) Over Time
# Purpose:  Visualize standardized beta coefficients and 
#           confidence intervals for cognitive, non-cognitive, 
#           and EA polygenic scores across multiple time points 
#           and domains (externalising & internalising).
#
# Approach:
#   - Extract standardized SEM estimates with confidence intervals
#   - Combine PGS components (EA, cognitive, non-cognitive)
#   - Compute variance explained (R²) per predictor
#   - Generate bar plots with error bars, grouped by pathway
#
# Outputs:  
#   - Bar plots of standardized betas with 95% CIs for each domain
#   - Annotated subtitles showing mean R² per predictor/pathway
#
# NOTE: 
#   - This script must be re-run for each outcome/domain (EXT vs INT) 
#     based on the specific model you want to plot.
#   - For the manuscript, plots were generated using the 
#     **parent-invariant** model for each outcome.
#
# Author:   Jose J. Morosoli
# Date:     31-07-2025
#############################################################

### --- Load Required Packages ---
library(dplyr)
library(ggplot2)

#------------------------------------------------------------
# 1. Cognitive & Non-cognitive standardized betas with CI
#------------------------------------------------------------
std_parentinv <- standardizedSolution(fit_parentinv_svy_EXTc, se = TRUE, ci = TRUE) %>%
  filter(op == "~") %>%
  mutate(
    Component = case_when(
      grepl("EA_cog", rhs) ~ "Cognitive",
      grepl("EA_noncog", rhs) ~ "Non-cognitive",
      TRUE ~ "Other"
    ),
    Pathway = case_when(
      grepl("_C_", rhs) ~ "Transmission",
      grepl("_M_", rhs) | grepl("_F_", rhs) ~ "Nurture",
      TRUE ~ "Other"
    ),
    Domain = case_when(
      grepl("EXT", lhs) ~ "Externalising",
      grepl("INT", lhs) ~ "Internalising",
      TRUE ~ "Other"
    ),
    Time = recode(lhs,
                  "BEXT_res" = "Age 3", "BINT_res" = "Age 3",
                  "CEXT_res" = "Age 5", "CINT_res" = "Age 5",
                  "DEXT_res" = "Age 7", "DINT_res" = "Age 7",
                  "FEXT_res" = "Age 14", "FINT_res" = "Age 14"
    ),
    Time = factor(Time, levels = c("Age 3", "Age 5", "Age 7", "Age 14"))
  )

std_parentinv_collapsed <- std_parentinv %>%
  group_by(lhs, Component, Pathway, Domain, Time) %>%
  summarise(
    est.std = first(est.std),
    ci_lower = first(ci.lower),
    ci_upper = first(ci.upper),
    .groups = "drop"
  )

#------------------------------------------------------------
# 2. Educational Attainment (EA) standardized betas with CI
#------------------------------------------------------------
std_EA <- standardizedSolution(ea2_fit_svy, se = TRUE, ci = TRUE) %>%
  filter(op == "~") %>%
  mutate(
    Component = "EA",
    Pathway = case_when(
      grepl("_C_", rhs) ~ "Transmission",
      grepl("_M_", rhs) | grepl("_F_", rhs) ~ "Nurture",
      TRUE ~ "Other"
    ),
    Domain = case_when(
      grepl("EXT", lhs) ~ "Externalising",
      grepl("INT", lhs) ~ "Internalising",
      TRUE ~ "Other"
    ),
    Time = recode(lhs,
                  "BEXT_res" = "Age 3", "BINT_res" = "Age 3",
                  "CEXT_res" = "Age 5", "CINT_res" = "Age 5",
                  "DEXT_res" = "Age 7", "DINT_res" = "Age 7",
                  "FEXT_res" = "Age 14", "FINT_res" = "Age 14"
    ),
    Time = factor(Time, levels = c("Age 3", "Age 5", "Age 7", "Age 14"))
  )

std_EA_collapsed <- std_EA %>%
  group_by(lhs, Component, Pathway, Domain, Time) %>%
  summarise(
    est.std = first(est.std),
    ci_lower = first(ci.lower),
    ci_upper = first(ci.upper),
    .groups = "drop"
  )

#------------------------------------------------------------
# 3. Combine All Components
#------------------------------------------------------------
std_all_combined <- bind_rows(std_parentinv_collapsed, std_EA_collapsed)

#------------------------------------------------------------
# 4. Compute R² per Predictor
#------------------------------------------------------------
r2_per_predictor <- std_all_combined %>%
  mutate(R2 = est.std^2)

#------------------------------------------------------------
# 5. Plotting Function with CI
#------------------------------------------------------------
plot_pgseffects_direct_indirect <- function(data, domain_name, r2_table) {
  
  r2_sub <- r2_table %>%
    filter(Domain == domain_name) %>%
    group_by(Component, Pathway) %>%
    summarise(MeanR2 = mean(R2), .groups = "drop")
  
  subtitle_text <- paste0(
    "Mean R²:",
    paste(
      vapply(
        split(r2_sub, factor(r2_sub$Component, levels = c("EA", "Cognitive", "Non-cognitive"))),
        function(df) {
          paste0(
            unique(df$Component), ": ",
            paste0(
              df$Pathway, " = ", sprintf("%.3f", df$MeanR2),
              collapse = " | "
            )
          )
        },
        character(1)
      ),
      collapse = "\n"
    )
  )
  
  ggplot(
    filter(data, Domain == domain_name),
    aes(x = Time, y = est.std,
        fill = Component,
        group = interaction(Component, Pathway))
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge(0.8),
      width = 0.7,
      color = "black"
    ) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(0.8),
      width = 0.2,
      linewidth = 0.6
    ) +
    facet_wrap(~ Pathway) +
    labs(
      title = paste(domain_name, ": Genetic transmission and Genetic nurture effects"),
      subtitle = subtitle_text,
      x = "Time Point",
      y = "Standardized Beta",
      fill = "Predictor"
    ) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 13, face = "bold"),
      legend.position = "top"
    )
}

#------------------------------------------------------------
# 6. Set Factor Levels for Consistent Bar Order
#------------------------------------------------------------
std_all_combined$Component <- factor(
  std_all_combined$Component,
  levels = c("EA", "Cognitive", "Non-cognitive")
)

#------------------------------------------------------------
# 7. Plot Externalising
#------------------------------------------------------------
plot_pgseffects_direct_indirect(std_all_combined, "Externalising", r2_per_predictor)

#------------------------------------------------------------
# 8. Plot Internalising
#------------------------------------------------------------
plot_pgseffects_direct_indirect(std_all_combined, "Internalising", r2_per_predictor)
