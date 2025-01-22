---
title: "Main analyses"
author: "Jose J. Morosoli"
date: "2025-01-21"
output: html_document
---

# Index

The current document describes:

1. Libraries and Import Data
2. Trio Models
3. Model comparison

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Section 1: Libraries and Import Data  

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

# Load required libraries
library(lavaan)
library(lme4)
library(psych)
library(tidyverse)
library(data.table)
library(tidyverse)
library(data.table)
library(reshape)
library(knitr)
library(kableExtra)

# Set global options
options(stringsAsFactors = FALSE)

# Import data
load(file = "data/myData_checked.RData")

```

# Section 2: Trio models
## Preliminary steps

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}


# DEFINE SHARED OBJECTS
outcome1 <- "BEXT_z"
outcome2 <- "CEXT_z"
outcome3 <- "DEXT_z" 
outcome4 <- "FEXT_z" 
outcome5 <- "BINT_z" 
outcome6 <- "CINT_z" 
outcome7 <- "DINT_z" 
outcome8 <- "FINT_z"

# Create object for covariates
PCs1_10 <- paste("PC", 1:10, "_C", sep="")
sex_10PCs <- c("SEX", PCs1_10)

```

## Modelling Educational Attainment

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

# What predictor? 
trait <- 'EA'
# Label model and specify trio PGSs
modelname <- paste(trait,"_","EXTINT",sep="") #   EXT INT
ind <- "C"
ind2 <- "M"
ind3 <- "F"
# Educational attainment
pgs1 <- paste0(trait, "_PGS_", ind)
pgs2 <- paste0(trait, "_PGS_", ind2)
pgs3 <- paste0(trait, "_PGS_", ind3)

```
## Null model

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are set to ZERO
myformula_null_ea1 <- paste(outcome1, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea2 <- paste(outcome2, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea3 <- paste(outcome3, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea4 <- paste(outcome4, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea5 <- paste(outcome5, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea6 <- paste(outcome6, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea7 <- paste(outcome7, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_null_ea8 <- paste(outcome8, "~","0*",pgs1,"+","0*",pgs2,"+","0*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_null_ea <- paste(# regressions 
  myformula_null_ea1,
  myformula_null_ea2,
  myformula_null_ea3,
  myformula_null_ea4,
  myformula_null_ea5,
  myformula_null_ea6,
  myformula_null_ea7,
  myformula_null_ea8,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_null_eaTidy <- noquote(strsplit(model_null_ea, "\n")[[1]])
# Fit the free model
null_ea_fit <- sem(model_null_eaTidy, data=myData_checked, missing = "ML")

```

## Base model

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}


#formulas where trio PGS effects are FREE to vary over time
myformula_base_ea1 <- paste(outcome1, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea2 <- paste(outcome2, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea3 <- paste(outcome3, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea4 <- paste(outcome4, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea5 <- paste(outcome5, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea6 <- paste(outcome6, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea7 <- paste(outcome7, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_base_ea8 <- paste(outcome8, "~",pgs1,"+",pgs2,"+",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
#specify the free model using these formulas
model_base_ea <- paste(# regressions 
  myformula_base_ea1,
  myformula_base_ea2,
  myformula_base_ea3,
  myformula_base_ea4,
  myformula_base_ea5,
  myformula_base_ea6,
  myformula_base_ea7,
  myformula_base_ea8,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_base_eaTidy <- noquote(strsplit(model_base_ea, "\n")[[1]])

# Fit the free model
base_ea_fit <- sem(model_base_eaTidy, data=myData_checked, missing = "ML")

```

## Parent-specific
### Within Externalising Difficulties

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

myformula_ea11 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be12a*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea12 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be12b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea13 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be12c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea14 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be12d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea15 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi13a*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea16 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi13b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea17 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi13c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea18 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi13d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_ae1 <- paste(# regressions 
  myformula_ea11,
  myformula_ea12,
  myformula_ea13,
  myformula_ea14,
  myformula_ea15,
  myformula_ea16,
  myformula_ea17,
  myformula_ea18,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_ae1Tidy <- noquote(strsplit(model_ae1, "\n")[[1]])

# Fit the free model
ae1_fit <- sem(model_ae1Tidy, data=myData_checked, missing = "ML")

```

### Within Internalising Difficulties

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

myformula_ea31 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be13a*",pgs3,"+",
                        
                        paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea32 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be13b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea33 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be13c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea34 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be13d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea35 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi12a*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea36 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi12b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea37 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi12c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea38 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi12d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_ae3 <- paste(# regressions 
  myformula_ea31,
  myformula_ea32,
  myformula_ea33,
  myformula_ea34,
  myformula_ea35,
  myformula_ea36,
  myformula_ea37,
  myformula_ea38,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_ae3Tidy <- noquote(strsplit(model_ae3, "\n")[[1]])
model_ae3Tidy

# Fit the free model
ae3_fit <- sem(model_ae3Tidy, data=myData_checked, missing = "ML")

```

## Time-specific
### Within externalising difficulties

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are FREE to vary over time
myformula_ea51 <- paste(outcome1, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea52 <- paste(outcome2, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea53 <- paste(outcome3, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea54 <- paste(outcome4, "~","be11*",pgs1,"+","be12*",pgs2,"+","be13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea55 <- paste(outcome5, "~","bi11a*",pgs1,"+","bi12a*",pgs2,"+","bi13a*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea56 <- paste(outcome6, "~","bi11b*",pgs1,"+","bi12b*",pgs2,"+","bi13b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea57 <- paste(outcome7, "~","bi11c*",pgs1,"+","bi12c*",pgs2,"+","bi13c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea58 <- paste(outcome8, "~","bi11d*",pgs1,"+","bi12d*",pgs2,"+","bi13d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_ae5 <- paste(# regressions 
  myformula_ea51,
  myformula_ea52,
  myformula_ea53,
  myformula_ea54,
  myformula_ea55,
  myformula_ea56,
  myformula_ea57,
  myformula_ea58,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_ae5Tidy <- noquote(strsplit(model_ae5, "\n")[[1]])

# Fit the free model
ae5_fit <- sem(model_ae5Tidy, data=myData_checked, missing = "ML")

```
### Within internalising difficulties

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are FREE to vary over time
myformula_ea71 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be13a*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea72 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be13b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea73 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be13c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea74 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be13d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea75 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea76 <- paste(outcome6, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea77 <- paste(outcome7, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea78 <- paste(outcome8, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi13*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
#specify the free model using these formulas
model_ae7 <- paste(# regressions 
  myformula_ea71,
  myformula_ea72,
  myformula_ea73,
  myformula_ea74,
  myformula_ea75,
  myformula_ea76,
  myformula_ea77,
  myformula_ea78,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_ae7Tidy <- noquote(strsplit(model_ae7, "\n")[[1]])

# Fit the free model
ae7_fit <- sem(model_ae7Tidy, data=myData_checked, missing = "ML")

```
## Modelling Cognitive and Non-cognitive components

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

# What predictor? 
trait1 <- "CP" 
trait2 <- "NCP" 
# Create object for covariates
PCs1_10 <- paste("PC", 1:10, "_C", sep="")
sex_10PCs <- c("SEX", PCs1_10)

# Label model and specify trio PGSs
modelname <- paste(trait1,"_",trait2,"_","EXTINT",sep="") #   EXT INT
ind <- "C"
ind2 <- "M"
ind3 <- "F"
pgs11 <- paste0(trait1, "_PGS_", ind)
pgs12 <- paste0(trait1, "_PGS_", ind2)
pgs13 <- paste0(trait1, "_PGS_", ind3)
pgs21 <- paste0(trait2, "_PGS_", ind)
pgs22 <- paste0(trait2, "_PGS_", ind2)
pgs23 <- paste0(trait2, "_PGS_", ind3)

```

## Null model

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are set to ZERO
myformula_null1 <- paste(outcome1, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null2 <- paste(outcome2, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null3 <- paste(outcome3, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null4 <- paste(outcome4, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null5 <- paste(outcome5, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null6 <- paste(outcome6, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null7 <- paste(outcome7, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")
myformula_null8 <- paste(outcome8, "~","0*",pgs11,"+","0*",pgs12,"+","0*",pgs13,"+",
                             "0*",pgs21,"+","0*",pgs22,"+","0*",pgs23,"+",
                             paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_null <- paste(# regressions 
  myformula_null1,
  myformula_null2,
  myformula_null3,
  myformula_null4,
  myformula_null5,
  myformula_null6,
  myformula_null7,
  myformula_null8,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_null_Tidy <- noquote(strsplit(model_null, "\n")[[1]])

# Fit the free model
null_fit <- sem(model_null_Tidy, data=myData_checked, missing = "ML")

```

## Base model

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are FREE to vary over time
myformula_base1 <- paste(outcome1, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base2 <- paste(outcome2, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base3 <- paste(outcome3, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base4 <- paste(outcome4, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base5 <- paste(outcome5, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base6 <- paste(outcome6, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base7 <- paste(outcome7, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")
myformula_base8 <- paste(outcome8, "~",pgs11,"+",pgs12,"+",pgs13,"+",
                         pgs21,"+",pgs22,"+",pgs23,"+",
                         paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_base <- paste(# regressions 
  myformula_base1,
  myformula_base2,
  myformula_base3,
  myformula_base4,
  myformula_base5,
  myformula_base6,
  myformula_base7,
  myformula_base8,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_baseTidy <- noquote(strsplit(model_base, "\n")[[1]])

# Fit the free model
base_fit <- sem(model_baseTidy, data=myData_checked, missing = "ML")

```

## Parent-specific
### Within externalising difficulties
- Cognitive component

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}
myformula_m11 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be12a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m12 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be12b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m13 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be12c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m14 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be12d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m15 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi13a*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m16 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi13b*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m17 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi13c*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m18 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi13d*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

# Specify the free model using these formulas
model_m1 <- paste(# regressions 
  myformula_m11,
  myformula_m12,
  myformula_m13,
  myformula_m14,
  myformula_m15,
  myformula_m16,
  myformula_m17,
  myformula_m18,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_m1Tidy <- noquote(strsplit(model_m1, "\n")[[1]])

# Fit the free model
m1_fit <- sem(model_m1Tidy, data=myData_checked, missing = "ML")

```

- Non-cognitive component

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

myformula_m21 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be13a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be22a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m22 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be13b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be22b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m23 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be13c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be22c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m24 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be13d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be22d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m25 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi13a*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m26 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi13b*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m27 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi13c*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m28 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi13d*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_m2 <- paste(# regressions 
  myformula_m21,
  myformula_m22,
  myformula_m23,
  myformula_m24,
  myformula_m25,
  myformula_m26,
  myformula_m27,
  myformula_m28,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_m2Tidy <- noquote(strsplit(model_m2, "\n")[[1]])

# Fit the free model
m2_fit <- sem(model_m2Tidy, data=myData_checked, missing = "ML")
```
### Within internalising difficulties
- Cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

myformula_m31 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be13a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m32 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be13b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m33 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be13c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m34 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be13d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m35 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi12a*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m36 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi12b*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m37 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi12c*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m38 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi12d*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_m3 <- paste(# regressions 
  myformula_m31,
  myformula_m32,
  myformula_m33,
  myformula_m34,
  myformula_m35,
  myformula_m36,
  myformula_m37,
  myformula_m38,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_m3Tidy <- noquote(strsplit(model_m3, "\n")[[1]])

# Fit the free model
m3_fit <- sem(model_m3Tidy, data=myData_checked, missing = "ML")

```
- Non-cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#### b2) NONCOG ####
myformula_m41 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be13a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m42 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be13b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m43 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be13c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m44 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be13d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m45 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi13a*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi22a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m46 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi13b*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi22b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m47 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi13c*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi22c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m48 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi13d*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi22d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_m4 <- paste(# regressions 
  myformula_m41,
  myformula_m42,
  myformula_m43,
  myformula_m44,
  myformula_m45,
  myformula_m46,
  myformula_m47,
  myformula_m48,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_m4Tidy <- noquote(strsplit(model_m4, "\n")[[1]])
model_m4Tidy

# Fit the free model
m4_fit <- sem(model_m4Tidy, data=myData_checked, missing = "ML")
print(summary(m4_fit),digits=2)

```

## Time-specific
### Within Externalising Difficulties
- Cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are FREE to vary over time
myformula_m51 <- paste(outcome1, "~","be11*",pgs11,"+","be12*",pgs12,"+","be13*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m52 <- paste(outcome2, "~","be11*",pgs11,"+","be12*",pgs12,"+","be13*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m53 <- paste(outcome3, "~","be11*",pgs11,"+","be12*",pgs12,"+","be13*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m54 <- paste(outcome4, "~","be11*",pgs11,"+","be12*",pgs12,"+","be13*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m55 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi13a*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m56 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi13b*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m57 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi13c*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m58 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi13d*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_m5 <- paste(# regressions 
  myformula_m51,
  myformula_m52,
  myformula_m53,
  myformula_m54,
  myformula_m55,
  myformula_m56,
  myformula_m57,
  myformula_m58,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_m5Tidy <- noquote(strsplit(model_m5, "\n")[[1]])

# Fit the free model
m5_fit <- sem(model_m5Tidy, data=myData_checked, missing = "ML")
```
- Non-cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}
#### c2) NONCOG ####
myformula_m61 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be13a*",pgs13,"+",
                       "be21*",pgs21,"+","be22*",pgs22,"+","be23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m62 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be13b*",pgs13,"+",
                       "be21*",pgs21,"+","be22*",pgs22,"+","be23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m63 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be13c*",pgs13,"+",
                       "be21*",pgs21,"+","be22*",pgs22,"+","be23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m64 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be13d*",pgs13,"+",
                       "be21*",pgs21,"+","be22*",pgs22,"+","be23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m65 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi13a*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m66 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi13b*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m67 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi13c*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m68 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi13d*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
#specify the free model using these formulas
model_m6 <- paste(# regressions 
  myformula_m61,
  myformula_m62,
  myformula_m63,
  myformula_m64,
  myformula_m65,
  myformula_m66,
  myformula_m67,
  myformula_m68,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")
#remove quotation marks and separate formulae
model_m6Tidy <- noquote(strsplit(model_m6, "\n")[[1]])
# Fit the free model
m6_fit <- sem(model_m6Tidy, data=myData_checked, missing = "ML")
```
### Within Internalising Difficulties
- Cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are FREE to vary over time
myformula_m71 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be13a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m72 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be13b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m73 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be13c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m74 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be13d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m75 <- paste(outcome5, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi13*",pgs13,"+",
                       "bi21a*",pgs21,"+","bi22a*",pgs22,"+","bi23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m76 <- paste(outcome6, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi13*",pgs13,"+",
                       "bi21b*",pgs21,"+","bi22b*",pgs22,"+","bi23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m77 <- paste(outcome7, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi13*",pgs13,"+",
                       "bi21c*",pgs21,"+","bi22c*",pgs22,"+","bi23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m78 <- paste(outcome8, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi13*",pgs13,"+",
                       "bi21d*",pgs21,"+","bi22d*",pgs22,"+","bi23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_m7 <- paste(# regressions 
  myformula_m71,
  myformula_m72,
  myformula_m73,
  myformula_m74,
  myformula_m75,
  myformula_m76,
  myformula_m77,
  myformula_m78,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")

#remove quotation marks and separate formulae
model_m7Tidy <- noquote(strsplit(model_m7, "\n")[[1]])

# Fit the free model
m7_fit <- sem(model_m7Tidy, data=myData_checked, missing = "ML")
```
- Non-cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}
myformula_m81 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be13a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be23a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m82 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be13b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be23b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m83 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be13c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be23c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m84 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be13d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be23d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m85 <- paste(outcome5, "~","bi11a*",pgs11,"+","bi12a*",pgs12,"+","bi13a*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m86 <- paste(outcome6, "~","bi11b*",pgs11,"+","bi12b*",pgs12,"+","bi13b*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m87 <- paste(outcome7, "~","bi11c*",pgs11,"+","bi12c*",pgs12,"+","bi13c*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_m88 <- paste(outcome8, "~","bi11d*",pgs11,"+","bi12d*",pgs12,"+","bi13d*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi23*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
#specify the free model using these formulas
model_m8 <- paste(# regressions 
  myformula_m81,
  myformula_m82,
  myformula_m83,
  myformula_m84,
  myformula_m85,
  myformula_m86,
  myformula_m87,
  myformula_m88,
  # correlated residuals
  "BINT_z ~~ CINT_z + DINT_z + FINT_z", 
  "CINT_z ~~ DINT_z + FINT_z",
  "DINT_z ~~ FINT_z",
  "BEXT_z ~~ CEXT_z + DEXT_z + FEXT_z", 
  "CEXT_z ~~ DEXT_z + FEXT_z",
  "DEXT_z ~~ FEXT_z",
  # correlated residuals between int and ext
  "BINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  "CINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "DINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z",
  "FINT_z ~~ BEXT_z + CEXT_z + DEXT_z + FEXT_z", 
  sep="\n")
#remove quotation marks and separate formulae
model_m8Tidy <- noquote(strsplit(model_m8, "\n")[[1]])
model_m8Tidy
# Fit the free model
m8_fit <- sem(model_m8Tidy, data=myData_checked, missing = "ML")

```
# Section 3: Model comparison

```
{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}

# Create function to extract ANOVA
extract_chisq <- function(myanova, modelname) {
  # extract data and save in dataframe
  chisq.diff <- myanova[2, "Chisq diff"]
  df.diff <- myanova[2, "Df diff"]
  pvalue <- myanova[2, "Pr(>Chisq)"]
  model_summary= as.data.frame(cbind(modelname, chisq.diff, df.diff, pvalue))
  # label significance level
  model_summary$pvalue_num=as.numeric(as.character(model_summary$pvalue))
  model_summary$pvalue_sig[model_summary$pvalue_num <= 0.05 & model_summary$pvalue_num > 0.01] = "*"
  model_summary$pvalue_sig[model_summary$pvalue_num <= 0.01] = "**"
  model_summary$pvalue_sig[model_summary$pvalue_num > 0.05] = "ns"
  return(model_summary)
}

# Can now specify two models to compare their fit, e.g. model <- anova(fit_model_free_c, fit_model_equal_c)

```
## Educational attainment

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#~~~~~~~~~ NULL MODEL ~~~~~~~#
model <- anova(null_ea_fit,base_fit)
(Chisq_0a <- extract_chisq(model, modelname = 'nul_ea'))
#~~~~~~~~~ FIRST STEP: PARENTS ~~~~~~~#
model <- anova(base_ea_fit,ae1_fit)
(Chisq_1 <- extract_chisq(model, modelname = 'p_ea_ext'))
model <- anova(base_ea_fit,ae3_fit)
(Chisq_2 <- extract_chisq(model, modelname = 'p_ea_int'))
#~~~~~~~~~ SECOND STEP: TIME ~~~~~~~#
model <- anova(base_ea_fit,ae5_fit)
(Chisq_3 <- extract_chisq(model, modelname = 't_ea_ext'))
model <- anova(base_ea_fit,ae7_fit)
(Chisq_4 <- extract_chisq(model, modelname = 't_ea_int'))
comparison1 <- rbind(Chisq_0a,Chisq_1,Chisq_2,Chisq_3,Chisq_4)
comparison1$p_fdr <- p.adjust(comparison1$pvalue_num, method="fdr")
comparison1

```

## Cognitive and non-cognitive components

```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#~~~~~~~~~ FIRST STAGE: PARENTS ~~~~~~~#
model <- anova(base_fit,m1_fit)
(Chisq_5 <- extract_chisq(model, modelname = 'p_cp_ext'))
model <- anova(base_fit,m2_fit)
(Chisq_6 <- extract_chisq(model, modelname = 'p_ncp_ext'))
model <- anova(base_fit,m3_fit)
(Chisq_7 <- extract_chisq(model, modelname = 'p_cp_int'))
model <- anova(base_fit,m4_fit)
(Chisq_8 <- extract_chisq(model, modelname = 'p_ncp_int'))
#~~~~~~~~~ SECOND STAGE: TIME ~~~~~~~#
model <- anova(base_fit,m5_fit)
(Chisq_9 <- extract_chisq(model, modelname = 't_cp_ext'))
model <- anova(base_fit,m6_fit)
(Chisq_10 <- extract_chisq(model, modelname = 't_ncp_ext'))
model <- anova(base_fit,m7_fit)
(Chisq_11 <- extract_chisq(model, modelname = 't_cp_int'))
model <- anova(base_fit,m8_fit)
(Chisq_12 <- extract_chisq(model, modelname = 't_ncp_int'))
#~~~~~~~~~ NULL MODEL ~~~~~~~#
model <- anova(null_fit,base_fit)
(Chisq_0b <- extract_chisq(model, modelname = 'null_cp_ncp'))
comparison2 <- rbind(Chisq_0b,Chisq_5,Chisq_6,Chisq_7,Chisq_8,Chisq_9,Chisq_10,Chisq_11,Chisq_12)
comparison2$p_fdr <- p.adjust(comparison2$pvalue_num, method="fdr")
comparison2

```
