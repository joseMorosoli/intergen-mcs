#############################################################
# Project: Intergenerational PGS Trio Models
# Script:   Subsample by removing genetic ancestry outliers
# Purpose:  Create a new subsample from the full trio sample 
#           by excluding participants who are outliers in 
#           principal component (PC) space based on the 
#           Mahalanobis distance. 
#
# Outputs:  - A new dataset (`aframe_svy_maha_na`) excluding 
#             ancestry outliers.
#           - A new survey design object (`svy_design`) 
#             using the updated sample.
#
# Usage:    Use this script to generate a clean subsample to 
#           replace the sample used in the SEM scripts for 
#           EA, cognitive, and non-cognitive PGS models. 
#           After running this script, you can re-run the 
#           previously cleaned model comparison scripts 
#           without changing any code.
#
# Author:   Jose J. Morosoli
# Date:     31-07-2025
#############################################################

### --- Load Data ---
load(file = "NATCOMMS_R1.RData")

### --- Identify Ancestry Outliers ---
# Extract only the PC columns (adjust if your PCs are named differently)
pc_matrix <- as.matrix(aframe_svy[, 33:52])

# Calculate the Mahalanobis distance from the multivariate mean
mahal_dist <- mahalanobis(
  x = pc_matrix,
  center = colMeans(pc_matrix),
  cov = cov(pc_matrix)
)

# Set significance level (e.g., 0.01) and compute cutoff
alpha <- 0.01
cutoff <- qchisq(1 - alpha, df = 20)  # 20 PCs → 20 degrees of freedom

# Flag individuals with distances beyond the threshold
outliers <- mahal_dist > cutoff

# Visualize Mahalanobis distance distribution
hist(mahal_dist, breaks = 50, main = "Mahalanobis Distances from PC Mean")
abline(v = cutoff, col = "red", lty = 2)

# Subset: exclude outliers
aframe_svy_maha <- aframe_svy[!outliers, ]

### --- Keep Only Data with Weights ---
aframe_svy_maha_na <- aframe_svy_maha[!is.na(aframe_svy_maha$FOVWT2), ]

### --- Define Survey Design ---
svy_design <- svydesign(
  ids = ~SPTN00,
  strata = ~PTTYPE2,
  weights = ~FOVWT2,
  data = aframe_svy_maha_na,
  nest = TRUE
)

### --- Compare Sample Sizes ---
dim(aframe_svy_maha_na)  # After excluding outliers
aframe_svy_na <- aframe_svy[!is.na(aframe_svy$FOVWT2), ]
dim(aframe_svy_na)       # Original sample
