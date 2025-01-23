# intergen-mcs
This is a repository for analyses described in "Intergenerational influences of the cognitive and non-cognitive components of parental education on externalising and internalising difficulties in childhood and adolescence" by Morosoli et al. (2024).

Link to preprint: https://osf.io/preprints/psyarxiv/hc3b5

# Overview

- 01-genetic-qc.Rmd - quality control and extraction of complete trios from genotype data.
- 02-calculate-pgs.Rmd - implementation of LDPred2 pipeline to calculate PGS.
- 03-data-preparation.Rmd - new variables and dataframes to use in the analysis phase.
- 04-main-analyses.Rmd - trio models and model comparison.
- 05-power-trio.Rmd - power analyses for polygenic score approach using trios.

The original pipeline to compute polygenic scores with LDpred2 can be found at https://github.com/AndreAllegrini/LDpred2

## Data
The data used in this study are from the Millennium Cohort Study (MCS). The MCS data are not publicly available but can be accessed by researchers through the UK Data Service. Researchers seeking access can contact the Centre for Longitudinal Studies at University College London (https://cls.ucl.ac.uk/), where detailed instructions for applying to use the data are provided. Access to the data requires registration and adherence to the terms of use set by the data providers, including agreements on confidentiality and ethical use.

## Software and code
Data analysis was performed using custom R (lavaan, lme4 packages in R v4.4.1) and Unix (PLINK 1.9, KING 2.2.7). We have shared code for main analyses as supplementary material for peer review purposes but we are currently developing a GitHub repository where code and commentary will be made publicly available.
