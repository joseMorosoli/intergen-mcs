#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   Descriptive statistics and sample checks
# Purpose:  Reproduce descriptive/sample analyses reported in
#           the manuscript but not currently included in the
#           main GitHub workflow, and provide statistics needed
#           for final editorial reporting.
#
# Inputs:   - ../NATCOMMS_R1.RData
#           - MCS phenotype files already used in 03-data-preparation.R
#
# Outputs:  - descriptive_table1_check.csv
#           - sample_comparison_tests.csv
#           - sample_comparison_wave_tests.csv
#           - ethnicity_comparison.csv
#           - ethnicity_comparison_diagnostic.csv
#           - sex_difference_tests.csv
#           - age_anova_results.csv
#           - age_anova_effect_sizes.csv
#           - SDQ_distribution_summary.csv
#           - Supplementary_Figure_1.png
#
# Notes:
# - Welch independent-samples t-tests are used (R t.test default).
# - Mean differences and their 95% CIs are retained as effect estimates.
# - `ez` is used only to reproduce the Greenhouse-Geisser-corrected
#   repeated-measures ANOVA reflected by the fractional dfs in the manuscript.
# - `effectsize` is used only to obtain partial eta-squared with 95% CIs.
#############################################################

library(foreign)
library(tidyverse)
library(e1071)
library(ez)
library(effectsize)

options(stringsAsFactors = FALSE)

results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE)

# ----------------------------------------------------------
# 1. LOAD DATA
# ----------------------------------------------------------

file_path <- "C:/UCL_DATA/00_LOCAL_MCS/"

load("../NATCOMMS_R1.RData")
aframe_svy_na <- aframe_svy[!is.na(aframe_svy$FOVWT2), ]

parentstructure <- read.spss(
  paste0(file_path, "GENDAC_PINGAULT_mcs_parent_structure_2021_05_12.sav"),
  to.data.frame = TRUE
)
CMstructure <- read.spss(
  paste0(file_path, "GENDAC_PINGAULT_mcs_cm_structure_2021_05_12.sav"),
  to.data.frame = TRUE
)
matches <- read.spss(
  paste0(file_path, "MDAC-2020-0016-05A-PINGAULT_unifed_ega_2021_12_13.sav"),
  to.data.frame = TRUE
)

parentstructure$Pingault_FID <- trimws(parentstructure$Pingault_FID)
CMstructure$Pingault_FID <- trimws(CMstructure$Pingault_FID)
matches$Pingault_FID <- trimws(matches$Pingault_FID)
matches$MFC <- trimws(matches$MFC)

# ----------------------------------------------------------
# 2. BASIC SAMPLE CHECKS + TABLE 1
# ----------------------------------------------------------

outcomes_raw <- c(
  "BEXT", "CEXT", "DEXT", "FEXT",
  "BINT", "CINT", "DINT", "FINT"
)

cat("Number of rows:", nrow(aframe_svy_na), "\n")
cat("Unique families:", length(unique(aframe_svy_na$Pingault_FID)), "\n")
cat("Duplicated family IDs:", sum(duplicated(aframe_svy_na$Pingault_FID)), "\n")
cat("Complete data for all eight outcome variables:",
    sum(complete.cases(aframe_svy_na[, outcomes_raw])), "\n")
print(table(aframe_svy_na$SEX_C, useNA = "ifany"))

age_info <- data.frame(
  Age = c(3, 5, 7, 14),
  EXT = c("BEXT", "CEXT", "DEXT", "FEXT"),
  INT = c("BINT", "CINT", "DINT", "FINT"),
  TOTAL = c("BEBDTOT", "CEBDTOT", "DDDEBDTOT", "FEBDTOT")
)

complete_all <- complete.cases(aframe_svy_na[, outcomes_raw])
sex_levels <- unique(na.omit(aframe_svy_na$SEX_C))

table1_check <- data.frame()

for (i in seq_len(nrow(age_info))) {
  for (s in sex_levels) {
    idx <- aframe_svy_na$SEX_C == s
    ext <- aframe_svy_na[[age_info$EXT[i]]]
    int <- aframe_svy_na[[age_info$INT[i]]]
    total <- aframe_svy_na[[age_info$TOTAL[i]]]
    both <- idx & !is.na(ext) & !is.na(int)
    
    table1_check <- rbind(
      table1_check,
      data.frame(
        Sex = as.character(s),
        Age = age_info$Age[i],
        N_ALL = sum(idx & complete_all),
        N_WAVE = sum(both),
        EXT_M = mean(ext[both]),
        EXT_SD = sd(ext[both]),
        INT_M = mean(int[both]),
        INT_SD = sd(int[both]),
        CCT_percent = mean(total[idx] >= 15, na.rm = TRUE) * 100
      )
    )
  }
}

print(table1_check)
write.csv(table1_check, file.path(results_dir, "descriptive_table1_check.csv"), row.names = FALSE)

# ----------------------------------------------------------
# 3. TRIO SAMPLE VS EXCLUDED MCS SAMPLE
# ----------------------------------------------------------

# Parent age comparison used in the archived analysis
included_age <- merge(
  aframe_svy_na,
  parentstructure[, c("Pingault_FID", "ADDAGB00")],
  by = "Pingault_FID"
)
excluded_parent <- parentstructure[
  !(parentstructure$Pingault_FID %in% aframe_svy_na$Pingault_FID),
]

age_test <- t.test(
  included_age$ADDAGB00,
  excluded_parent$ADDAGB00
)

# Parental education comparison used in the archived analysis
excluded_parent$ADACAQ00[excluded_parent$ADACAQ00 == 96] <- NA
excluded_parent$ADACAQ00[excluded_parent$ADACAQ00 == 95] <- 3

education_test <- t.test(
  c(aframe_svy_na$ADACAQ00_M, aframe_svy_na$ADACAQ00_F),
  excluded_parent$ADACAQ00
)

# Create externalising/internalising scores in the excluded sample
excluded_child <- CMstructure[
  !(CMstructure$Pingault_FID %in% aframe_svy_na$Pingault_FID),
]

vars_numeric <- c(
  "BEMOTION", "BPEER", "CEMOTION", "CPEER",
  "DDEMOTION", "DDPEER", "FEMOTION", "FPEER",
  "BCONDUCT", "BHYPER", "CCONDUCT", "CHYPER",
  "DDCONDUCT", "DDHYPER", "FCONDUCT", "FHYPER"
)
excluded_child[vars_numeric] <- lapply(
  excluded_child[vars_numeric],
  function(x) as.numeric(as.character(x))
)

excluded_child <- excluded_child %>%
  mutate(
    BINT = BEMOTION + BPEER,
    CINT = CEMOTION + CPEER,
    DINT = DDEMOTION + DDPEER,
    FINT = FEMOTION + FPEER,
    BEXT = BCONDUCT + BHYPER,
    CEXT = CCONDUCT + CHYPER,
    DEXT = DDCONDUCT + DDHYPER,
    FEXT = FCONDUCT + FHYPER
  )

# Parent age and parental education comparisons
sample_tests <- list(
  Parent_age = age_test,
  Parental_education = education_test
)

sample_comparison_tests <- do.call(
  rbind,
  lapply(names(sample_tests), function(name) {
    x <- sample_tests[[name]]
    data.frame(
      Comparison = name,
      Mean_included = unname(x$estimate[1]),
      Mean_excluded = unname(x$estimate[2]),
      Mean_difference = unname(x$estimate[1] - x$estimate[2]),
      t = unname(x$statistic),
      df = unname(x$parameter),
      p = x$p.value,
      CI_lower = x$conf.int[1],
      CI_upper = x$conf.int[2]
    )
  })
)

print(sample_comparison_tests)
write.csv(
  sample_comparison_tests,
  file.path(results_dir, "sample_comparison_tests.csv"),
  row.names = FALSE
)

# Reproduce the archived wave-specific child outcome comparisons
sample_comparison_wave_tests <- data.frame()

for (v in outcomes_raw) {
  tt <- t.test(
    aframe_svy_na[[v]],
    excluded_child[[v]]
  )
  
  sample_comparison_wave_tests <- rbind(
    sample_comparison_wave_tests,
    data.frame(
      Outcome = v,
      Mean_included = unname(tt$estimate[1]),
      Mean_excluded = unname(tt$estimate[2]),
      Mean_difference = unname(tt$estimate[1] - tt$estimate[2]),
      t = unname(tt$statistic),
      df = unname(tt$parameter),
      p = tt$p.value,
      CI_lower = tt$conf.int[1],
      CI_upper = tt$conf.int[2]
    )
  )
}

print(sample_comparison_wave_tests)
write.csv(
  sample_comparison_wave_tests,
  file.path(results_dir, "sample_comparison_wave_tests.csv"),
  row.names = FALSE
)

cat(
  "Archived-style average EXT mean difference:",
  mean(sample_comparison_wave_tests$Mean_difference[1:4]),
  "\n"
)
cat(
  "Archived-style average INT mean difference:",
  mean(sample_comparison_wave_tests$Mean_difference[5:8]),
  "\n"
)

# Ethnicity comparison: cohort members in the study vs excluded cohort members
excluded_matches_child <- matches[
  !(matches$Pingault_FID %in% aframe_svy_na$Pingault_FID) &
    matches$MFC == "C",
]

included_white <- trimws(as.character(aframe_svy_na$ETHNICITY_C)) == "White"
excluded_white_child <- trimws(as.character(excluded_matches_child$ETHNICITY)) == "White"

n_white_included <- sum(included_white, na.rm = TRUE)
n_total_included <- sum(!is.na(included_white))
n_white_excluded_child <- sum(excluded_white_child, na.rm = TRUE)
n_total_excluded_child <- sum(!is.na(excluded_white_child))

ethnicity_test <- prop.test(
  x = c(n_white_included, n_white_excluded_child),
  n = c(n_total_included, n_total_excluded_child),
  correct = FALSE
)

ethnicity_comparison <- data.frame(
  White_included = n_white_included,
  N_included = n_total_included,
  White_excluded = n_white_excluded_child,
  N_excluded = n_total_excluded_child,
  Proportion_included = unname(ethnicity_test$estimate[1]),
  Proportion_excluded = unname(ethnicity_test$estimate[2]),
  Difference = unname(ethnicity_test$estimate[1] - ethnicity_test$estimate[2]),
  Chi_square = unname(ethnicity_test$statistic),
  df = unname(ethnicity_test$parameter),
  p = ethnicity_test$p.value,
  CI_lower = ethnicity_test$conf.int[1],
  CI_upper = ethnicity_test$conf.int[2]
)

print(ethnicity_comparison)
write.csv(
  ethnicity_comparison,
  file.path(results_dir, "ethnicity_comparison.csv"),
  row.names = FALSE
)

# Diagnostic only: reproduce the archived comparison, which did not restrict
# the excluded genetic sample to cohort members.
excluded_matches_all <- matches[
  !(matches$Pingault_FID %in% aframe_svy_na$Pingault_FID),
]
excluded_white_all <- trimws(as.character(excluded_matches_all$ETHNICITY)) == "White"

n_white_excluded_all <- sum(excluded_white_all, na.rm = TRUE)
n_total_excluded_all <- sum(!is.na(excluded_white_all))

ethnicity_test_archived <- prop.test(
  x = c(n_white_included, n_white_excluded_all),
  n = c(n_total_included, n_total_excluded_all),
  correct = FALSE
)

ethnicity_comparison_diagnostic <- data.frame(
  Comparison = c("Child_only", "Archived_all_excluded_family_members"),
  White_included = c(n_white_included, n_white_included),
  N_included = c(n_total_included, n_total_included),
  White_excluded = c(n_white_excluded_child, n_white_excluded_all),
  N_excluded = c(n_total_excluded_child, n_total_excluded_all),
  Proportion_included = c(
    unname(ethnicity_test$estimate[1]),
    unname(ethnicity_test_archived$estimate[1])
  ),
  Proportion_excluded = c(
    unname(ethnicity_test$estimate[2]),
    unname(ethnicity_test_archived$estimate[2])
  ),
  Chi_square = c(
    unname(ethnicity_test$statistic),
    unname(ethnicity_test_archived$statistic)
  ),
  df = c(
    unname(ethnicity_test$parameter),
    unname(ethnicity_test_archived$parameter)
  ),
  p = c(
    ethnicity_test$p.value,
    ethnicity_test_archived$p.value
  )
)

print(ethnicity_comparison_diagnostic)
write.csv(
  ethnicity_comparison_diagnostic,
  file.path(results_dir, "ethnicity_comparison_diagnostic.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 4. SEX DIFFERENCES AT EACH AGE
# ----------------------------------------------------------

sex_difference_tests <- data.frame()

for (v in outcomes_raw) {
  tt <- t.test(aframe_svy_na[[v]] ~ aframe_svy_na$SEX_C)
  
  sex_difference_tests <- rbind(
    sex_difference_tests,
    data.frame(
      Outcome = v,
      Group_1 = names(tt$estimate)[1],
      Mean_1 = unname(tt$estimate[1]),
      Group_2 = names(tt$estimate)[2],
      Mean_2 = unname(tt$estimate[2]),
      Mean_difference = unname(tt$estimate[1] - tt$estimate[2]),
      t = unname(tt$statistic),
      df = unname(tt$parameter),
      p = tt$p.value,
      CI_lower = tt$conf.int[1],
      CI_upper = tt$conf.int[2]
    )
  )
}

print(sex_difference_tests)
write.csv(sex_difference_tests, file.path(results_dir, "sex_difference_tests.csv"), row.names = FALSE)

# ----------------------------------------------------------
# 5. DIFFERENCES ACROSS AGE
# ----------------------------------------------------------

# Complete cases are used separately for externalising and internalising,
# matching a conventional within-subjects repeated-measures ANOVA.

ext_long <- aframe_svy_na %>%
  select(Pingault_FID, BEXT, CEXT, DEXT, FEXT) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(
    cols = c(BEXT, CEXT, DEXT, FEXT),
    names_to = "Age",
    values_to = "Score"
  ) %>%
  mutate(
    Pingault_FID = factor(Pingault_FID),
    Age = factor(Age, levels = c("BEXT", "CEXT", "DEXT", "FEXT"))
  )

int_long <- aframe_svy_na %>%
  select(Pingault_FID, BINT, CINT, DINT, FINT) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(
    cols = c(BINT, CINT, DINT, FINT),
    names_to = "Age",
    values_to = "Score"
  ) %>%
  mutate(
    Pingault_FID = factor(Pingault_FID),
    Age = factor(Age, levels = c("BINT", "CINT", "DINT", "FINT"))
  )

ext_anova <- ezANOVA(
  data = ext_long,
  dv = Score,
  wid = Pingault_FID,
  within = Age,
  type = 3,
  detailed = TRUE
)

int_anova <- ezANOVA(
  data = int_long,
  dv = Score,
  wid = Pingault_FID,
  within = Age,
  type = 3,
  detailed = TRUE
)

print(ext_anova)
print(int_anova)

ext_aov <- ext_anova$ANOVA[ext_anova$ANOVA$Effect == "Age", ]
int_aov <- int_anova$ANOVA[int_anova$ANOVA$Effect == "Age", ]

ext_gg <- ext_anova$`Sphericity Corrections`[
  ext_anova$`Sphericity Corrections`$Effect == "Age",
]
int_gg <- int_anova$`Sphericity Corrections`[
  int_anova$`Sphericity Corrections`$Effect == "Age",
]

ext_df1_gg <- ext_aov$DFn * ext_gg$GGe
ext_df2_gg <- ext_aov$DFd * ext_gg$GGe
int_df1_gg <- int_aov$DFn * int_gg$GGe
int_df2_gg <- int_aov$DFd * int_gg$GGe

age_anova_results <- data.frame(
  Outcome = c("Externalising", "Internalising"),
  F = c(ext_aov$F, int_aov$F),
  df1_GG = c(ext_df1_gg, int_df1_gg),
  df2_GG = c(ext_df2_gg, int_df2_gg),
  p_GG = c(ext_gg$`p[GG]`, int_gg$`p[GG]`),
  generalized_eta_squared = c(ext_aov$ges, int_aov$ges)
)

print(age_anova_results)
write.csv(age_anova_results, file.path(results_dir, "age_anova_results.csv"), row.names = FALSE)

ext_eta <- as.data.frame(
  effectsize::F_to_eta2(
    ext_aov$F,
    df = ext_df1_gg,
    df_error = ext_df2_gg,
    ci = 0.95,
    alternative = "two.sided"
  )
)
int_eta <- as.data.frame(
  effectsize::F_to_eta2(
    int_aov$F,
    df = int_df1_gg,
    df_error = int_df2_gg,
    ci = 0.95,
    alternative = "two.sided"
  )
)

age_anova_effect_sizes <- bind_rows(
  cbind(Outcome = "Externalising", ext_eta),
  cbind(Outcome = "Internalising", int_eta)
)

print(age_anova_effect_sizes)
write.csv(
  age_anova_effect_sizes,
  file.path(results_dir, "age_anova_effect_sizes.csv"),
  row.names = FALSE
)

# ----------------------------------------------------------
# 6. SDQ DISTRIBUTION CHECKS FOR SUPPLEMENTARY TABLE/FIGURE
# ----------------------------------------------------------

outcomes_res <- c(
  "BEXT_res", "CEXT_res", "DEXT_res", "FEXT_res",
  "BINT_res", "CINT_res", "DINT_res", "FINT_res"
)

dist_summary <- data.frame()

for (i in seq_along(outcomes_res)) {
  x_res <- aframe_svy_na[[outcomes_res[i]]]
  x_raw <- aframe_svy_na[[outcomes_raw[i]]]
  
  dist_summary <- rbind(
    dist_summary,
    data.frame(
      Variable = outcomes_raw[i],
      Mean_res = mean(x_res, na.rm = TRUE),
      SD_res = sd(x_res, na.rm = TRUE),
      Min_res = min(x_res, na.rm = TRUE),
      Max_res = max(x_res, na.rm = TRUE),
      Skewness_res = e1071::skewness(x_res, na.rm = TRUE),
      Kurtosis_res = e1071::kurtosis(x_res, na.rm = TRUE),
      Proportion_zero_raw = mean(x_raw == 0, na.rm = TRUE)
    )
  )
}

print(dist_summary)
write.csv(
  dist_summary,
  file.path(results_dir, "SDQ_distribution_summary.csv"),
  row.names = FALSE
)

raw_plot <- aframe_svy_na %>%
  select(all_of(outcomes_raw)) %>%
  pivot_longer(
    everything(),
    names_to = "Variable",
    values_to = "Score"
  ) %>%
  mutate(Type = "Raw")

res_plot <- aframe_svy_na %>%
  select(all_of(outcomes_res)) %>%
  pivot_longer(
    everything(),
    names_to = "Variable",
    values_to = "Score"
  ) %>%
  mutate(
    Variable = sub("_res$", "", Variable),
    Type = "Residual"
  )

plot_data <- bind_rows(raw_plot, res_plot) %>%
  mutate(
    Variable = factor(Variable, levels = outcomes_raw),
    Type = factor(Type, levels = c("Raw", "Residual"))
  )

hist_plot <- ggplot(plot_data, aes(x = Score)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30) +
  geom_density() +
  facet_grid(Type ~ Variable, scales = "free_x") +
  labs(
    x = "Score",
    y = "Density"
  ) +
  theme_minimal()

print(hist_plot)

ggsave(
  file.path(results_dir, "Supplementary_Figure_1.png"),
  plot = hist_plot,
  width = 16,
  height = 5,
  dpi = 300
)
