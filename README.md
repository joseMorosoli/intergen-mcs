# Genetic Transmission & Nurture Trio Analyses

This is a repository for analyses described in "Intergenerational influences of the cognitive and non-cognitive genetic components of parental education on externalising and internalising behaviour in childhood and adolescence" by Morosoli et al. (2024).

Link to preprint: https://osf.io/preprints/psyarxiv/hc3b5_v1

# Overview

- 01-genetic-qc.R - quality control and extraction of complete trios from genotype data.
- 02-calculate-pgs.R - implementation of LDPred2 pipeline to calculate PGS.
- 03-data-preparation.R - new variables and dataframes to use in the analysis phase.
- 04-EA-analyses.R - trio models and model comparison for EA PGS.
- 05-cog_noncog-analyses.R - trio models and model comparison for Cog and NonCog PGS.
- 06-model-comparison-fdr.R - combine all model comparisons and apply FDR.
- 07-sensitivity-analyses.R - create subsample to run sensitivity analyses.
- 08-phenotypic-analyses.R - run analyses including observed parental education.
- 09-power-trio.R - power analyses for polygenic score approach using trios.
- 10-plots-with-cis.R - plotting results.
- how-install-lavaansurvey.txt - instructions to install archived R package.

Some of the scripts have been automatically formatted and summarise to display the original logic while preserving readability.

Feel free to contact the authors if you have any comments: j.morosoli@ucl.ac.uk

The original pipeline to compute polygenic scores with LDpred2 can be found at https://github.com/AndreAllegrini/LDpred2

## Data
The data used in this study are from the Millennium Cohort Study (MCS). The MCS data are not publicly available but can be accessed by researchers through the UK Data Service. Researchers seeking access can contact the Centre for Longitudinal Studies at University College London (https://cls.ucl.ac.uk/), where detailed instructions for applying to use the data are provided. Access to the data requires registration and adherence to the terms of use set by the data providers, including agreements on confidentiality and ethical use.

## Software and code
Data analysis was performed using custom R (lavaan, lme4 packages in R v4.4.1) and Unix (PLINK 1.9, KING 2.2.7).
