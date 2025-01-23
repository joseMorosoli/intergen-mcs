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
