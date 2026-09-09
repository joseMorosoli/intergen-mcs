# Genetic Transmission & Nurture Trio Analyses

This repository contains the analysis code for the manuscript:

**"Education-related polygenic scores show different associations with youth mental health in parent-offspring trios"**  
Morosoli et al.

Preprint: https://osf.io/preprints/psyarxiv/hc3b5_v2

## Overview

The analyses use mother-father-child trios from the Millennium Cohort Study to examine associations between education-related polygenic scores and externalising and internalising behaviour across childhood and adolescence.

The main analysis scripts are:

- `01-genetic-qc.R` - genetic quality control and extraction of complete parent-offspring trios.
- `02-calculate-pgs.R` - implementation of the LDpred2 pipeline used to calculate polygenic scores.
- `03-data-preparation.R` - preparation of phenotypic, genetic and analysis variables.
- `04-EA-analyses.R` - structural equation models and model comparisons for the educational attainment PGS.
- `05-cog_noncog-analyses.R` - structural equation models and model comparisons for cognitive and non-cognitive PGSs.
- `06-model-comparison-fdr.R` - combines constraint-based model comparisons and applies false discovery rate correction.
- `07-sensitivity-analyses.R` - sensitivity analyses excluding trios with outlying genetic principal component profiles.
- `08-phenotypic-analyses.R` - analyses of observed parental education and adjustment of genetic models for parental education.
- `09-power-trio.R` - Monte Carlo power analyses for the trio-based polygenic score models.
- `10-plots-with-cis.R` - generates Figures 1 and 2 using parent-invariant, age-specific model estimates.
- `11-get-estimates-cis-EA.R` - extracts standardised estimates and confidence intervals for the educational attainment PGS models.
- `12-get-estimates-cis-CNC.R` - extracts standardised estimates and confidence intervals for the cognitive and non-cognitive PGS models, including supplementary Bonferroni-adjusted confidence intervals.
- `13-descriptive-statistics-and-sample-checks.R` - descriptive statistics, sample comparisons and preliminary sex- and age-based analyses.

Additional files used in polygenic score calculation are provided in the `files/` directory:

- `exclHapMap.sh` - variant filtering step used in polygenic score preparation.
- `hm3plus-pos.list` - extended HapMap3 variant list used for LDpred2.
- `ldpred2_auto_inf_qc.R` - LDpred2 scoring and quality-control script.
- `subLDpred2-single.sh` - shell script used to submit LDpred2 analyses.

`how-to-install-lavaan.survey.txt` provides instructions for installing the archived version of `lavaan.survey` used in the analyses.

Generated analysis outputs and figures are stored locally in the `results/` directory but are not version-controlled.

## Data

The data used in this study are from the Millennium Cohort Study (MCS). MCS phenotype data are available to bona fide researchers through the UK Data Service under the relevant access conditions.

Genotype data are available to qualified researchers through a managed access process subject to approval by the cohort's data access committee. Individual-level genotype data and derived individual-level polygenic scores cannot be distributed through this repository.

Further information on MCS data access is available from the Centre for Longitudinal Studies:

https://cls.ucl.ac.uk/data-access-training/

## Software and code

Analyses were conducted primarily in R v4.4.1. Structural equation models were fitted using `lavaan` and complex survey design and sampling weights were incorporated using `lavaan.survey`.

Genetic data processing used PLINK 1.9 and KING 2.2.7. Polygenic scores were calculated using LDpred2.

The original LDpred2 polygenic score pipeline on which the scoring workflow was based is available at:

https://github.com/AndreAllegrini/LDpred2

## Contact

For questions about the analyses or code, please contact:

José J. Morosoli  
j.morosoli@ucl.ac.uk