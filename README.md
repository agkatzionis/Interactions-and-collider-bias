# Interactions-and-collider-bias

This repository contains R code in support of the paper "Relationship between Collider Bias and Interactions on the Log-Additive Scale" (A. Gkatzionis et al. 2023).

Suppose that we wish to investigate the association between an exposure X and an outcome Y, but we can only observe exposure and outcome values conditional on the value of a third variable S. For example, S may represent selection into a study (1: yes, 0: no), in which case we can only observe X and Y conditional on S = 1. In this case, we are interested in the unconditional X-Y association, but based on the available data we can only estimate the conditional X-Y association given S = 1. If X and Y affect S, the conditional and unconditional associations will differ; this phenomenon is known as **collider bias**.

<p align="center">
  <img src="https://github.com/agkatzionis/Interactions-and-collider-bias/assets/46974026/185b1e42-72b5-45db-833a-0f3e3680887e"/>
</p>

We investigate the relationship between the magnitude of collider bias (i.e. the difference between the conditional and unconditional X-Y associations) and interactions between the exposure and outcome in their effects on S. In our manuscript, we show that if the outcome Y follows a logistic regression or Poisson regression model and the collider S follows a log-additive model, the bias is exactly equal to the exposure-outcome interaction. In addition, if the outcome follows a linear regression model and the collider S follows a log-additive model, the bias is a linearfunction of the exposure-outcome interaction, with slope equal to the standard error of the outcome regression.

In this repository we conduct an asymptotic study to explore the relationship between collider bias and exposure-outcome interactions in a log-additive model for S, when the collider S does not follow a log-additive model. We show that when the exposure X is binary, the magnitude of collider bias in the exposure-outcome association is still equal or proportional to the X-Y interaction in a (misspecified) log-additive model for S. On the other hand, when the exposure X is not binary, thesse results do not hold.

The files "Interactions and Collider Bias - Asymptotic Study.R" and "Asymptotic_Results.RData" contain the R code and results used to run the asymptotic study and create teh Figres for our paper.

The files "Interactions and Collider Bias - ALSPAC.R" and "ALSPAC_Results.RData" contain the R code and results for our real-data application using data from the [Avon Longitudinal Study of Parents and Children](http://www.bristol.ac.uk/alspac/) (ALSPAC). Note that reproducing our real-data analysis requires access to the ALSPAC data, which can be requested directly from ALSPAC at: http://www.bristol.ac.uk/alspac/researchers/access/.

