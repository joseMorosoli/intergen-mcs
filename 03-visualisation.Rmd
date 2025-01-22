---
title: "03-visualisation"
author: "Jose J. Morosoli"
date: "2025-01-22"
output: html_document
---

The current document describes:

1. How to re-run best fitting models.
2. Data handling and how to create the plotting function.
3. Creating plots.

# Section 1: How to re-run best fitting models.
## Educational Attainment
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

# formulas where trio PGS effects are FREE to vary over time (ext) but constrained across parents (both)
myformula_ea_best1 <- paste(outcome1, "~","be11a*",pgs1,"+","be12a*",pgs2,"+","be12a*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best2 <- paste(outcome2, "~","be11b*",pgs1,"+","be12b*",pgs2,"+","be12b*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best3 <- paste(outcome3, "~","be11c*",pgs1,"+","be12c*",pgs2,"+","be12c*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best4 <- paste(outcome4, "~","be11d*",pgs1,"+","be12d*",pgs2,"+","be12d*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best5 <- paste(outcome5, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best6 <- paste(outcome6, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best7 <- paste(outcome7, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
myformula_ea_best8 <- paste(outcome8, "~","bi11*",pgs1,"+","bi12*",pgs2,"+","bi12*",pgs3,"+",paste(sex_10PCs, collapse= "+"), sep="")
#specify the free model using these formulas
model_ea_best <- paste(# regressions 
  myformula_ea_best1,
  myformula_ea_best2,
  myformula_ea_best3,
  myformula_ea_best4,
  myformula_ea_best5,
  myformula_ea_best6,
  myformula_ea_best7,
  myformula_ea_best8,
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
model_ea_bestTidy <- noquote(strsplit(model_ea_best, "\n")[[1]])
model_ea_bestTidy
# Fit the free model
model_ea_best_fit <- sem(model_ea_bestTidy, data=myData_checked, missing = "ML")
```
## Cognitive and non-cognitive component
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

#formulas where trio PGS effects are FREE to vary over time (ext) but constrained across parents (both)
myformula_comp_best1 <- paste(outcome1, "~","be11a*",pgs11,"+","be12a*",pgs12,"+","be12a*",pgs13,"+",
                       "be21a*",pgs21,"+","be22a*",pgs22,"+","be22a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best2 <- paste(outcome2, "~","be11b*",pgs11,"+","be12b*",pgs12,"+","be12b*",pgs13,"+",
                       "be21b*",pgs21,"+","be22b*",pgs22,"+","be22b*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best3 <- paste(outcome3, "~","be11c*",pgs11,"+","be12c*",pgs12,"+","be12c*",pgs13,"+",
                       "be21c*",pgs21,"+","be22c*",pgs22,"+","be22c*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best4 <- paste(outcome4, "~","be11d*",pgs11,"+","be12d*",pgs12,"+","be12d*",pgs13,"+",
                       "be21d*",pgs21,"+","be22d*",pgs22,"+","be22d*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best5 <- paste(outcome5, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi12*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22a*",pgs22,"+","bi22a*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best6 <- paste(outcome6, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi12*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi22*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best7 <- paste(outcome7, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi12*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi22*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")
myformula_comp_best8 <- paste(outcome8, "~","bi11*",pgs11,"+","bi12*",pgs12,"+","bi12*",pgs13,"+",
                       "bi21*",pgs21,"+","bi22*",pgs22,"+","bi22*",pgs23,"+",
                       paste(sex_10PCs, collapse= "+"), sep="")

#specify the free model using these formulas
model_comp_best <- paste(# regressions 
  myformula_comp_best1,
  myformula_comp_best2,
  myformula_comp_best3,
  myformula_comp_best4,
  myformula_comp_best5,
  myformula_comp_best6,
  myformula_comp_best7,
  myformula_comp_best8,
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
model_comp_bestTidy <- noquote(strsplit(model_comp_best, "\n")[[1]])

# Fit the free model
comp_best_fit <- sem(model_comp_bestTidy, data=myData_checked, missing = "ML")
```
# Section 2: Data handling and how to create the plotting function
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

# Create dataframe with estimates from best fitting model
myresults <- data.frame(beta = rep(NA, 72), lower.ci = rep(NA, 72), upper.ci = rep(NA, 72),
                        outcome = rep(c(rep('Externalising Score',12),rep('Internalising Score',12)),3), 
                        trait.name = rep(c(rep('Educational Attainment',24),
                                           rep('Cognitive Skills',24),
                                           rep('Non-cognitive Skills',24))),
                        age = rep(rep(c(rep('3',3),
                                        rep('5',3),
                                        rep('7',3),
                                        rep('14',3)), 2),3),
                        fam = rep(c('Child','Mother','Father'),24),
                        model = c(rep('Multivariate',72))
)

# EA model
multi_ae <- parameterEstimates(model_ea_best_fit) # Save output
multi_ae <- multi_ae[grep('PGS',multi_ae$rhs),] # Keep only PGS betas
multi_ae <- multi_ae[grep('_z',multi_ae$lhs),]
# Save into results dataframe
myresults$beta[1:24] <- multi_ae$est
myresults$lower.ci[1:24] <- multi_ae$ci.lower
myresults$upper.ci[1:24] <- multi_ae$ci.upper
# Cog and NonCog model
multi <- parameterEstimates(comp_best_fit) # Save output
multi$rhs <- gsub(pattern = 'NCP', replacement = 'NCS', x = multi$rhs) # Relabel predictors to make extraction easier
multi <- multi[grep('PGS',multi$rhs),] # Keep only PGS betas
multi <- multi[grep('_z',multi$lhs),]
multi_cp <- multi[grep('CP',multi$rhs),] # Subset cog estimates
multi_ncp <- multi[grep('NCS',multi$rhs),] # Subset noncog estimates
# Save into results dataframe
myresults$beta[25:48] <- multi_cp$est
myresults$lower.ci[25:48] <- multi_cp$ci.lower
myresults$upper.ci[25:48] <- multi_cp$ci.upper
myresults$beta[49:72] <- multi_ncp$est
myresults$lower.ci[49:72] <- multi_ncp$ci.lower
myresults$upper.ci[49:72] <- multi_ncp$ci.upper
# Reshape data for ggplot
myresults_multi <- myresults[myresults$model=='Multivariate',]
myresults_multi$id <- paste(substr(myresults_multi$fam, 1, 5),
                            substr(myresults_multi$trait.name, 1, 3),
                            substr(myresults_multi$outcome, 1, 3),
                            myresults_multi$age, sep = "_")
names(myresults_multi)[1:3] <- c('beta_multi','lower.ci_multi','upper.ci_multi')
myresults_merged <- myresults_multi

# Note: Plot only one outcome
# Create a factor to fix the order of the facets in the figure
myresults_merged$age_f <- factor(myresults_merged$age, levels=c("14","7", "5",  "3")) #age
myresults_merged$trait.name_f <- factor(myresults_merged$trait.name, levels=c("Educational Attainment", "Cognitive Skills", "Non-cognitive Skills")) #PGS trait
myresults_merged$outcome_f <- factor(myresults_merged$outcome, levels=c("Externalising Score",'Internalising Score')) #Mental health difficulty

# Rename levels
levels(myresults_merged$trait.name_f ) <- c("EA PGS", "Cog PGS", "NonCog PGS")

# Transform data from long to wide format
#library(data.table)
setDT(myresults_merged)
myresults_lite <- myresults_merged[,c(1:3,7,10:12)]
beta.10PCs.wide <- dcast(myresults_lite, outcome_f+age_f+trait.name_f~fam,
                         value.var = c('beta_multi','lower.ci_multi','upper.ci_multi'))

# Update names based on function
names(beta.10PCs.wide)[4] <- 'beta_model_multi_child'
names(beta.10PCs.wide)[7] <- 'lower_CI_model_child'
names(beta.10PCs.wide)[10] <- 'upper_CI_model_child'
names(beta.10PCs.wide)[5] <- 'beta_model_multi_father'
names(beta.10PCs.wide)[8] <- 'lower_CI_model_father'
names(beta.10PCs.wide)[11] <- 'upper_CI_model_father'
names(beta.10PCs.wide)[6] <- 'beta_model_multi_mother'
names(beta.10PCs.wide)[9] <- 'lower_CI_model_mother'
names(beta.10PCs.wide)[12] <- 'upper_CI_model_mother'
beta.10PCs.wide_onlyext <- beta.10PCs.wide[beta.10PCs.wide$outcome_f=='Externalising Score',]
```
## Create plotting function
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

plot_beta2 = function(model) {  
  
  #Set positions on the graph for each value
  nudge4 <- position_nudge(y = -.2, x = 0)    # child multivariate
  nudge5 <- position_nudge(y = .2, x = 0)   # father multivariate
  nudge6 <- position_nudge(y = 0, x = 0)  # mother multivariate
  
  #Set the values for axes in ggplot and assign colours/shapes for each value
  myplot <- 
    ggplot(beta.10PCs.wide_onlyext, aes(x=beta_model_multi_child, y=reorder(age_f, desc(age_f)))) + 
    scale_color_manual(name="Data", values=c("1" = "blue3", "2" ="mediumvioletred", "4" = "darkgrey", "5" = "slategray1", "6" = "pink"), 
                       labels = c("Child PGS","Parental PGS")) + 
    scale_shape_manual(name="Data", values=c("1" = 15,"2" = 17,"3" = 19), 
                       labels = c("Child PGS","Parental PGS")) +
    
    #Set point and error bar values for bi- and multi-variate models for each PGS
    #Using position, colour and shape values assigned above
    #Educational attainment PGS
    # Child
    geom_point(aes(x=beta_model_multi_child, y=reorder(age_f, desc(age_f)), color="1", shape="1"), size=2, position = nudge4) + 
    geom_errorbar(aes(xmin=lower_CI_model_child, xmax=upper_CI_model_child,color="1"), width=.2, show.legend=FALSE, position = nudge4) +
    # Parental
    geom_point(aes(x=beta_model_multi_father, y=reorder(age_f, desc(age_f)), color="2", shape="2"), size=2, position = nudge5) + 
    geom_errorbar(aes(xmin=lower_CI_model_father, xmax=upper_CI_model_father,color="2"), width=.2, show.legend=FALSE, position = nudge5) +
    #Cognitive PGS
    # Child
    geom_point(aes(x=beta_model_multi_child, y=reorder(age_f, desc(age_f)), color="1", shape="1"), size=2, position = nudge4) + 
    geom_errorbar(aes(xmin=lower_CI_model_child, xmax=upper_CI_model_child,color="1"), width=.2, show.legend=FALSE, position = nudge4) +
    # Parental
    geom_point(aes(x=beta_model_multi_father, y=reorder(age_f, desc(age_f)), color="2", shape="2"), size=2, position = nudge5) + 
    geom_errorbar(aes(xmin=lower_CI_model_father, xmax=upper_CI_model_father,color="2"), width=.2, show.legend=FALSE, position = nudge5) +
    #Non-cognitive PGS
    # Child
    geom_point(aes(x=beta_model_multi_child, y=reorder(age_f, desc(age_f)), color="1", shape="1"), size=2, position = nudge4) + 
    geom_errorbar(aes(xmin=lower_CI_model_child, xmax=upper_CI_model_child,color="1"), width=.2, show.legend=FALSE, position = nudge4) +
    # parental
    geom_point(aes(x=beta_model_multi_father, y=reorder(age_f, desc(age_f)), color="2", shape="2"), size=2, position = nudge5) + 
    geom_errorbar(aes(xmin=lower_CI_model_father, xmax=upper_CI_model_father,color="2"), width=.2, show.legend=FALSE, position = nudge5) +
    #Set plot elements
    xlim(-0.25,0.09) + #x-axis range
    #ylim(-0.3,0.09) +
    geom_vline(xintercept = 0, linetype="dotted") + #zero line on x-axis
    theme_light(base_size=12) + 
    theme(axis.title.y=element_blank(),
          # axis.title.x=element_blank(),
          axis.title=element_text(size=10),
          legend.title=element_blank(),
          panel.grid=element_blank(),
          plot.background=element_blank()) +
    xlab("Regression Coefficients (95% CI)") +
    ylab("Child Age at Outcome (years)") +
    theme(legend.position="right", legend.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20)) +
    coord_flip()
  
  # Set panels with each trait (EA,cog,non-cog) as rows and mental health difficulty (total, or internal vs. external) as columns
  myplot_panel = myplot + facet_grid(rows=vars(outcome_f), cols=vars(trait.name_f), scales = "free", space = "free") + 
    theme(strip.text = element_text(size=14, colour = 'black'),
          strip.background = element_rect(colour="gray", fill="white"))
  
  plot(myplot_panel)
}
```
# Section 3: Run plotting function and store result
```{r echo=F, results='hide', error=F, warning=F, message=F, eval=F}

plot_10PCs <- plot_beta2(beta.10PCs.wide_onlyext)


```
