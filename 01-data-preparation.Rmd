---
title: "Data preparation"
author: "Jose J. Morosoli"
date: "2025-01-21"
output: html_document
---

# Index

The current document includes step-by-step instructions on how to:

1. Convert PDF to Word.
2. Import Word to R (string).
3. Data cleaning.
4. Spell check text files using OpenAI.
5. Run basic text mining analyses on enquiry documents.

# Preliminary steps

All source documents should be located in a sub-directories named 'data', 'data/pdf' and 'data/word'.


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r echo=T, results='hide', error=F, warning=F, message=F, eval=F}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
