#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   Combine ANOVA model comparisons & apply FDR
# Purpose:  Combine model comparison results across 
#           educational attainment (EA), cognitive (cog), 
#           and non-cognitive (noncog) models. 
#           Apply FDR correction and export results.
#
# Inputs:   - Chi-square difference test results from EA, 
#             cognitive, and non-cognitive SEM models
# Outputs:  - Combined comparison table with FDR-corrected p-values
#           - Excel file containing results
#
# Author:   Jose J. Morosoli
# Date:     31-07-2025
#############################################################

### --- Combine all ANOVA model comparisons ---
comparison <- rbind(
  Chisq_2, Chisq_3, Chisq_4, Chisq_5, Chisq_6, Chisq_7,
  Chisq_1a, Chisq_2a, Chisq_3a, Chisq_4a, Chisq_5a, Chisq_6a,
  Chisq_1b, Chisq_2b, Chisq_3b, Chisq_4b, Chisq_5b, Chisq_6b
)

### --- Apply FDR correction ---
comparison$p_fdr <- p.adjust(comparison3$pvalue_num, method = "fdr")

### --- Display results ---
comparison
round(comparison$p_fdr, 3)
round(comparison$pvalue_num, 3)

### --- Post-hoc comparison between cognitive and non-cognitive models ---
cognoncog_comp <- rbind(Chisq_7a, Chisq_7b)
cognoncog_comp

### --- Export results to Excel ---
library(readxl)
write.xlsx(comparison, file = "comparison-results.xlsx")
