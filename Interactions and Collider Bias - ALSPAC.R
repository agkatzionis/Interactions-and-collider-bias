
##########   INTERACTIONS AND COLLIDER BIAS   ##########

## This file contains the R code used to run the real-data 
## application included in the paper "Relationship Between
## Collider Bias and Interactions on the Log-additive Scale".

## We use data from the Avon Longitudinal Study of Parents
## and Children (ALSPAC). We explore associations of six 
## maternal traits with offspring sex in the full ALSPAC 
## sample, as well as in the TF4 and CCU subsamples, aiming 
## to illustrate the relationship between collider bias and
## exposure-outcome interactions on the log-additive scale.

## Access to the ALSPAC dataset for this project was obtained
## under application B4189. Any requests to access the data 
## should made directly to ALSPAC, through the URL:
## http://www.bristol.ac.uk/alspac/.

## Here, we only provide the R commands we used to analyze 
## the data for the paper's real-data application.

## ---------- LOAD ALSPAC DATA ---------- ##

## Set up R.
options(digits = 4)
library(xtable)
setwd("...")

## Load the ALSPAC data.
alspac <- read.csv("...")
dim(alspac)   ## 15645 rows and 61 columns.

## ---------- EXTRACT VARIABLES ---------- ##

## Our outcome variable is kz021, representing child sex.
sex <- alspac$kz021
sex[sex < 0] <- NA
sum(is.na(sex))   ## 607, not much.
sum(is.na(sex)) / length(sex)
sex.01 <- sex
sex.01[sex == 2] <- 0
## 1: male, 0/2: female.

## We consider six exposure variables.

## Mother's age.
mother.age <- alspac$mz028b
mother.age[mother.age < 0] <- NA
sum(is.na(mother.age))   ## 1585.
sum(is.na(mother.age)) / length(mother.age)

## Mother's BMI.
mother.bmi <- alspac$dw042
mother.bmi[mother.bmi < 0] <- NA
sum(is.na(mother.bmi))   ## 3993.
sum(is.na(mother.bmi)) / length(mother.bmi)

## Mother's education.
mother.educ <- alspac$c645a
mother.educ[mother.educ < 0] <- NA
sum(is.na(mother.educ))   ## 3167.
sum(is.na(mother.educ)) / length(mother.educ)

## Mother's smoking.
mother.smk <- alspac$b663
mother.smk[mother.smk < 0] <- NA
sum(is.na(mother.smk))   ## 2301.
sum(is.na(mother.smk)) / length(mother.smk)

## Mother's depression status.
mother.depr <- alspac$d171
mother.depr[mother.depr < 0] <- NA
sum(is.na(mother.depr))   ## 3257.
sum(is.na(mother.depr)) / length(mother.depr)

## Mother's gestation.
mother.gest <- alspac$bestgest
mother.gest[mother.gest < 0] <- NA
sum(is.na(mother.gest))   ## 1048.
sum(is.na(mother.gest)) / length(mother.gest)

## ---------- TF4/CCU PARTICIPATION ---------- ##

## Missingness in TF4 can be obtained from the variable
## "FJ002a", which represents month/year of the TF4 visit.
miss1 <- which(alspac$FJ002a < 0 | is.na(alspac$FJ002a))
sel1 <- which(!(1:nrow(alspac) %in% miss1))
sel1.01 <- rep(0, nrow(alspac))
sel1.01[sel1] <- 1
sum(sel1.01)

## Missingness in CCU can be extracted from the variable 
## "CCU0006", representing questionnaire return status.
sel2 <- which(alspac$CCU0006 == 1)
miss2 <- which(!(1:nrow(alspac) %in% sel2))
sel2.01 <- rep(0, nrow(alspac))
sel2.01[sel2] <- 1
sum(sel2.01)


## Does TF4/CCU participation differ by sex?

## Proportions of boys and girls in the ALSPAC sample.
mean(sex.01, na.rm = TRUE)   ## 51.1% boys, 48.9% girls.

## Proportions of male participants in each dataset:
sum(sex.01[sel1.01 == 1] == 1, na.rm = TRUE) / (sum(sex.01[sel1.01 == 1] == 1, na.rm = TRUE) + sum(sex.01[sel1.01 == 1] == 0, na.rm = TRUE))   ## 43.6% boys in TF4.
sum(sex.01[sel2.01 == 1] == 1, na.rm = TRUE) / (sum(sex.01[sel2.01 == 1] == 1, na.rm = TRUE) + sum(sex.01[sel2.01 == 1] == 0, na.rm = TRUE))   ## 39.1% boys in CCU.

## Participation rates for boys and girls in each stage.
sum(sel1.01[sex.01 == 1], na.rm = TRUE) / sum(sex.01 == 1, na.rm = TRUE)   ## 29.5% participation rate for boys at TF4.
sum(sel1.01[sex.01 == 0], na.rm = TRUE) / sum(sex.01 == 0, na.rm = TRUE)   ## 40.0% participation rate for girls at TF4.
sum(sel2.01[sex.01 == 1], na.rm = TRUE) / sum(sex.01 == 1, na.rm = TRUE)   ## 22.0% participation rate for boys at CCU
sum(sel2.01[sex.01 == 0], na.rm = TRUE) / sum(sex.01 == 0, na.rm = TRUE)   ## 35.9% participation rate for girls at CCU

## Here, we study collider bias due to TF4/CCU participation. 
## There are further biases in the dataset, e.g. due to maternal 
## exposure data missing for some individuals. We will ignore 
## these biases and only run complete-case analyses. Missingness
## rates for each maternal exposure are calculated later in this
## file and reported in Supplementary Table 2 of the paper.

## ---------- MAIN ANALYSIS - MARGINAL ---------- ##

## We run the main analysis, fitting a separate model for the
## association of each maternal trait with child sex.

## Store results here.
Marg_Results <- matrix(0, 6, 9)
colnames(Marg_Results) <- c("Beta_T", "StdErr_T", "Pval_T", "Beta_TF4", "StdErr_TF4", "Pval_TF4", "Beta_CCU", "StdErr_CCU", "Pval_CCU")
rownames(Marg_Results) <- c("Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")

## Run the analyses.
Marg_Results[1, 1:3] <- summary(glm(sex.01 ~ mother.age, family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[1, 4:6] <- summary(glm(sex.01[sel1] ~ mother.age[sel1], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[1, 7:9] <- summary(glm(sex.01[sel2] ~ mother.age[sel2], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[2, 1:3] <- summary(glm(sex.01 ~ mother.educ, family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[2, 4:6] <- summary(glm(sex.01[sel1] ~ mother.educ[sel1], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[2, 7:9] <- summary(glm(sex.01[sel2] ~ mother.educ[sel2], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[3, 1:3] <- summary(glm(sex.01 ~ mother.bmi, family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[3, 4:6] <- summary(glm(sex.01[sel1] ~ mother.bmi[sel1], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[3, 7:9] <- summary(glm(sex.01[sel2] ~ mother.bmi[sel2], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[4, 1:3] <- summary(glm(sex.01 ~ mother.depr, family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[4, 4:6] <- summary(glm(sex.01[sel1] ~ mother.depr[sel1], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[4, 7:9] <- summary(glm(sex.01[sel2] ~ mother.depr[sel2], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[5, 1:3] <- summary(glm(sex.01 ~ mother.smk, family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[5, 4:6] <- summary(glm(sex.01[sel1] ~ mother.smk[sel1], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[5, 7:9] <- summary(glm(sex.01[sel2] ~ mother.smk[sel2], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[6, 1:3] <- summary(glm(sex.01 ~ mother.gest, family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[6, 4:6] <- summary(glm(sex.01[sel1] ~ mother.gest[sel1], family = binomial))$coefficients[2, c(1, 2, 4)]
Marg_Results[6, 7:9] <- summary(glm(sex.01[sel2] ~ mother.gest[sel2], family = binomial))$coefficients[2, c(1, 2, 4)]

## Compare results - this is Table 1 of the paper.
Marg_Results

## Mother's age at delivery and gestational age were associated
## with child sex in all three samples. Mother's education was not
## associated with child sex in the full ALSPAC sample but was seen
## to associate with child sex in both TF4 and CCU. Maternal smoking
## was associated with child sex in the CCU sample but not in the 
## TF4 sample or in the full ALSPAC sample, while BMI and depression 
## before pregnancy exhibited no association with child sex in any
## of the three regression analyses. These results suggest collider
## bias may be affecting the association of maternal education and
## smoking with child sex.

## ---------- BIAS ASSESSMENT - MARGINAL ---------- ##

## Since we have data for both TF4/CCU participants and for
## non-participants, we can analyze participation in terms 
## of maternal traits and child sex.

## First, fit a log-additive model with interactions.

## Store results here.
Marg_LogAddI_TF4 <- matrix(0, 6, 9)
colnames(Marg_LogAddI_TF4) <- c("Delta1", "StdErr1", "Pval1", "Delta2", "StdErr2", "Pval2", "Delta3", "StdErr3", "Pval3")
rownames(Marg_LogAddI_TF4) <- c("Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")
Marg_LogAddI_CCU <- matrix(0, 6, 9)
colnames(Marg_LogAddI_CCU) <- c("Delta1", "StdErr1", "Pval1", "Delta2", "StdErr2", "Pval2", "Delta3", "StdErr3", "Pval3")
rownames(Marg_LogAddI_CCU) <- c("Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")

## Do it for the TF4 selection mechanism.
Marg_LogAddI_TF4[1, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.age + sex.01 * mother.age, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_TF4[2, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.educ + sex.01 * mother.educ, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_TF4[3, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.bmi + sex.01 * mother.bmi, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_TF4[4, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.depr + sex.01 * mother.depr, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_TF4[5, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.smk + sex.01 * mother.smk, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_TF4[6, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.gest + sex.01 * mother.gest, family = poisson))$coef[2:4, -3]))

## Do it for the CCU selection mechanism.
Marg_LogAddI_CCU[1, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.age + sex.01 * mother.age, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_CCU[2, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.educ + sex.01 * mother.educ, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_CCU[3, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.bmi + sex.01 * mother.bmi, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_CCU[4, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.depr + sex.01 * mother.depr, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_CCU[5, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.smk + sex.01 * mother.smk, family = poisson))$coef[2:4, -3]))
Marg_LogAddI_CCU[6, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.gest + sex.01 * mother.gest, family = poisson))$coef[2:4, -3]))

## Assess results - this is Table 2 of the paper.
Marg_LogAddI_TF4
Marg_LogAddI_CCU

## Compute the bias for comparison.
Bias_TF4 <- Marg_Results[, 4] - Marg_Results[, 1]
Bias_TF4
Bias_CCU <- Marg_Results[, 7] - Marg_Results[, 1]
Bias_CCU

## Bias figures for both sub-samples are similar to the 
## corresponding delta3 estimates in the log-additive model.



## As an additional analysis, we fit logistic models 
## without interactions instead of log-additive ones.

## Store results here.
Marg_LogitN_TF4 <- matrix(0, 6, 6)
colnames(Marg_LogitN_TF4) <- c("Delta1", "StdErr1", "Pval1", "Delta2", "StdErr2", "Pval2")
rownames(Marg_LogitN_TF4) <- c("Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")
Marg_LogitN_CCU <- matrix(0, 6, 6)
colnames(Marg_LogitN_CCU) <- c("Delta1", "StdErr1", "Pval1", "Delta2", "StdErr2", "Pval2")
rownames(Marg_LogitN_CCU) <- c("Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")

## Do it for the TF4 selection mechanism.
sel1.01 <- rep(0, nrow(alspac))
sel1.01[sel1] <- 1
Marg_LogitN_TF4[1, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.age, family = binomial))$coef[2:3, -3]))
Marg_LogitN_TF4[2, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.educ, family = binomial))$coef[2:3, -3]))
Marg_LogitN_TF4[3, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.bmi, family = binomial))$coef[2:3, -3]))
Marg_LogitN_TF4[4, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.depr, family = binomial))$coef[2:3, -3]))
Marg_LogitN_TF4[5, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.smk, family = binomial))$coef[2:3, -3]))
Marg_LogitN_TF4[6, ] <- c(t(summary(glm(sel1.01 ~ sex.01 + mother.gest, family = binomial))$coef[2:3, -3]))

## Do it for the CCU selection mechanism.
sel2.01 <- rep(0, nrow(alspac))
sel2.01[sel2] <- 1
Marg_LogitN_CCU[1, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.age, family = binomial))$coef[2:3, -3]))
Marg_LogitN_CCU[2, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.educ, family = binomial))$coef[2:3, -3]))
Marg_LogitN_CCU[3, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.bmi, family = binomial))$coef[2:3, -3]))
Marg_LogitN_CCU[4, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.depr, family = binomial))$coef[2:3, -3]))
Marg_LogitN_CCU[5, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.smk, family = binomial))$coef[2:3, -3]))
Marg_LogitN_CCU[6, ] <- c(t(summary(glm(sel2.01 ~ sex.01 + mother.gest, family = binomial))$coef[2:3, -3]))

## Assess results - this is Table 3 of the paper.
Marg_LogitN_TF4
Marg_LogitN_CCU

## Little connection between parameter values and bias.
## In fact, this can be misleading in suggesting that
## all six maternal traits may suffer from collider bias.

## ---------- MERGE WITH SUPPL 2 ---------- ##

## At the same time, some of the maternal traits also have 
## missing values. We will restrict the analysis to pregnancies 
## with complete maternal data and child sex.


## Check if TF4/CCU participation differs by sex. If this is 
## the case, participation will be causally downstream of sex,
## which is concerning for collider bias.


## Missingness rates for each maternal exposure in each subsample.

## Sample size, among children with recorded sex.
sum(sel1.01 == 1 & !is.na(sex))
sum(sel2.01 == 1 & !is.na(sex))

## Missingness counts for TF4.
sum(sel1.01 == 1 & is.na(mother.age) & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.bmi) & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.educ) & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.smk) & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.depr) & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.gest) & !is.na(sex))

## Missingness rates for TF4.
sum(sel1.01 == 1 & is.na(mother.age) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.bmi) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.educ) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.smk) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.depr) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex))
sum(sel1.01 == 1 & is.na(mother.gest) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex))

## Missingness counts for CCU.
sum(sel2.01 == 1 & is.na(mother.age) & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.bmi) & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.educ) & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.smk) & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.depr) & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.gest) & !is.na(sex))

## Missingness rates for CCU.
sum(sel2.01 == 1 & is.na(mother.age) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.bmi) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.educ) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.smk) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.depr) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex))
sum(sel2.01 == 1 & is.na(mother.gest) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex))

## ---------- SUPPLEMENTARY TABLE 2 ---------- ##

## Compute missingness rates for each maternal exposure,
## counted only among pregnancies which resulted in birth.

## Create the table.
Miss_table <- matrix(0, 6, 6)
rownames(Miss_table) <- c("Age", "Education", "BMI", "Depression", "Smoking", "Gest Age")
colnames(Miss_table) <- c("All_prop", "All_n", "TF4_prop", "TF4_n", "CCU_prop", "CCU_n")

## Identify miscarriages and early terminations.
miscarriages <- which(is.na(sex))

## All-ALSPAC missingness rates.
Miss_table[, 1] <- c(sum(is.na(mother.age[- miscarriages])) / (nrow(alspac) - length(miscarriages)), 
                     sum(is.na(mother.educ[- miscarriages])) / (nrow(alspac) - length(miscarriages)),
                     sum(is.na(mother.bmi[- miscarriages])) / (nrow(alspac) - length(miscarriages)),
                     sum(is.na(mother.depr[- miscarriages])) / (nrow(alspac) - length(miscarriages)),
                     sum(is.na(mother.smk[- miscarriages])) / (nrow(alspac) - length(miscarriages)),
                     sum(is.na(mother.gest[- miscarriages])) / (nrow(alspac) - length(miscarriages)))

## All-ALSPAC missingness counts.
Miss_table[, 2] <- c(sum(is.na(mother.age[- miscarriages])), 
                     sum(is.na(mother.educ[- miscarriages])),
                     sum(is.na(mother.bmi[- miscarriages])), 
                     sum(is.na(mother.depr[- miscarriages])), 
                     sum(is.na(mother.smk[- miscarriages])), 
                     sum(is.na(mother.gest[- miscarriages])))

## TF4 missingness rates.
Miss_table[, 3] <- c(sum(sel1.01 == 1 & is.na(mother.age) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.educ) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.bmi) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.depr) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.smk) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.gest) & !is.na(sex)) / sum(sel1.01 == 1 & !is.na(sex)))

## TF4 missingness counts.
Miss_table[, 4] <- c(sum(sel1.01 == 1 & is.na(mother.age) & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.educ) & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.bmi) & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.depr) & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.smk) & !is.na(sex)),
                     sum(sel1.01 == 1 & is.na(mother.gest) & !is.na(sex)))

## CCUMissingness rates.
Miss_table[, 5] <- c(sum(sel2.01 == 1 & is.na(mother.age) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.educ) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.bmi) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.depr) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.smk) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.gest) & !is.na(sex)) / sum(sel2.01 == 1 & !is.na(sex)))

## CCU Missingness counts.
Miss_table[, 6] <- c(sum(sel2.01 == 1 & is.na(mother.age) & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.educ) & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.bmi) & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.depr) & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.smk) & !is.na(sex)),
                     sum(sel2.01 == 1 & is.na(mother.gest) & !is.na(sex)))

## Report results.
Miss_table

## ---------- MAIN ANALYSIS - JOINT ---------- ##

## In practice one would not fit six separate regression 
## models but would include all the maternal traits into 
## a single model. We do so here and repeat our analysis 
## for the coefficients of the multiple regression model. 

## Number of individuals to be included in the analyses.
sum(!(is.na(sex.01)) & !(is.na(mother.age)) & !(is.na(mother.educ)) & !(is.na(mother.bmi)) &  !(is.na(mother.depr)) & !(is.na(mother.smk)) & !(is.na(mother.gest)))
sum(!(is.na(sex.01[sel1])) & !(is.na(mother.age[sel1])) & !(is.na(mother.educ[sel1])) & !(is.na(mother.bmi[sel1])) &  !(is.na(mother.depr[sel1])) & !(is.na(mother.smk[sel1])) & !(is.na(mother.gest[sel1])))
sum(!(is.na(sex.01[sel2])) & !(is.na(mother.age[sel2])) & !(is.na(mother.educ[sel2])) & !(is.na(mother.bmi[sel2])) &  !(is.na(mother.depr[sel2])) & !(is.na(mother.smk[sel2])) & !(is.na(mother.gest[sel2])))

## Missingness rates in each dataset.
1 - sum(!(is.na(sex.01)) & !(is.na(mother.age)) & !(is.na(mother.educ)) & !(is.na(mother.bmi)) &  !(is.na(mother.depr)) & !(is.na(mother.smk)) & !(is.na(mother.gest))) / length(sex.01)
1 - sum(!(is.na(sex.01[sel1])) & !(is.na(mother.age[sel1])) & !(is.na(mother.educ[sel1])) & !(is.na(mother.bmi[sel1])) &  !(is.na(mother.depr[sel1])) & !(is.na(mother.smk[sel1])) & !(is.na(mother.gest[sel1]))) / length(sex.01)
1 - sum(!(is.na(sex.01[sel2])) & !(is.na(mother.age[sel2])) & !(is.na(mother.educ[sel2])) & !(is.na(mother.bmi[sel2])) &  !(is.na(mother.depr[sel2])) & !(is.na(mother.smk[sel2])) & !(is.na(mother.gest[sel2]))) / length(sex.01)
1 - sum(!(is.na(sex.01[sel1])) & !(is.na(mother.age[sel1])) & !(is.na(mother.educ[sel1])) & !(is.na(mother.bmi[sel1])) &  !(is.na(mother.depr[sel1])) & !(is.na(mother.smk[sel1])) & !(is.na(mother.gest[sel1]))) / length(sex.01[sel1])
1 - sum(!(is.na(sex.01[sel2])) & !(is.na(mother.age[sel2])) & !(is.na(mother.educ[sel2])) & !(is.na(mother.bmi[sel2])) &  !(is.na(mother.depr[sel2])) & !(is.na(mother.smk[sel2])) & !(is.na(mother.gest[sel2]))) / length(sex.01[sel2])

## Define selection coefficients.
sel.j0 <- which(!(is.na(sex.01)) & !(is.na(mother.age)) & !(is.na(mother.educ)) & !(is.na(mother.bmi)) &  !(is.na(mother.depr)) & !(is.na(mother.smk)) & !(is.na(mother.gest)))
sel.j1 <- which(!(is.na(sex.01)) & !(is.na(mother.age)) & !(is.na(mother.educ)) & !(is.na(mother.bmi)) &  !(is.na(mother.depr)) & !(is.na(mother.smk)) & !(is.na(mother.gest)) & !(1:nrow(alspac) %in% miss1))
sel.j2 <- which(!(is.na(sex.01)) & !(is.na(mother.age)) & !(is.na(mother.educ)) & !(is.na(mother.bmi)) &  !(is.na(mother.depr)) & !(is.na(mother.smk)) & !(is.na(mother.gest)) & alspac$CCU0006 == 1)
c(length(sel.j0), length(sel.j1), length(sel.j2))


## This runs the X-Y regression on individuals with complete 
## maternal data and stores the results.
Joint_results <- Marg_Results
Joint_results[, 1:3] <- summary(glm(sex.01[sel.j0] ~ mother.age[sel.j0] + mother.educ[sel.j0] + mother.bmi[sel.j0] + mother.depr[sel.j0] + mother.smk[sel.j0] + mother.gest[sel.j0], family = binomial))$coefficients[2:7, c(1, 2, 4)]
Joint_results[, 4:6] <- summary(glm(sex.01[sel.j1] ~ mother.age[sel.j1] + mother.educ[sel.j1] + mother.bmi[sel.j1] + mother.depr[sel.j1] + mother.smk[sel.j1] + mother.gest[sel.j1], family = binomial))$coefficients[2:7, c(1, 2, 4)]
Joint_results[, 7:9] <- summary(glm(sex.01[sel.j2] ~ mother.age[sel.j2] + mother.educ[sel.j2] + mother.bmi[sel.j2] + mother.depr[sel.j2] + mother.smk[sel.j2] + mother.gest[sel.j2], family = binomial))$coefficients[2:7, c(1, 2, 4)]

## Report results - this is Supplementary Table 3 of the paper.
Joint_results

## ---------- BIAS ASSESSMENT - JOINT ---------- ##

## Now fit various models for participation and assess
## how well they explain the magnitude of selection bias.

## Explore bias on the log-additive scale with interactions.
Joint_LogAddI_TF4 <- summary(glm(sel1.01[sel.j0] ~ sex.01[sel.j0] + mother.age[sel.j0] + mother.educ[sel.j0] + mother.bmi[sel.j0] + mother.depr[sel.j0] + mother.smk[sel.j0] + mother.gest[sel.j0] + sex.01[sel.j0] * mother.age[sel.j0] + sex.01[sel.j0] * mother.educ[sel.j0] + sex.01[sel.j0] * mother.bmi[sel.j0] + sex.01[sel.j0] * mother.depr[sel.j0] + sex.01[sel.j0] * mother.smk[sel.j0] + sex.01[sel.j0] * mother.gest[sel.j0], family = poisson))$coef[-1, -3]
rownames(Joint_LogAddI_TF4) <- c("Child Sex", "Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age",
                                 "Age x Sex", "Education x Sex", "BMI x Sex", "Depression x Sex", "Smoking x Sex", "Gestation x Sex")
Joint_LogAddI_CCU <- summary(glm(sel2.01[sel.j0] ~ sex.01[sel.j0] + mother.age[sel.j0] + mother.educ[sel.j0] + mother.bmi[sel.j0] + mother.depr[sel.j0] + mother.smk[sel.j0] + mother.gest[sel.j0] + sex.01[sel.j0] * mother.age[sel.j0] + sex.01[sel.j0] * mother.educ[sel.j0] + sex.01[sel.j0] * mother.bmi[sel.j0] + sex.01[sel.j0] * mother.depr[sel.j0] + sex.01[sel.j0] * mother.smk[sel.j0] + sex.01[sel.j0] * mother.gest[sel.j0], family = poisson))$coef[-1, -3]
rownames(Joint_LogAddI_CCU) <- c("Child Sex", "Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age",
                                 "Age x Sex", "Education x Sex", "BMI x Sex", "Depression x Sex", "Smoking x Sex", "Gestation x Sex")

## Report results - this is Supplementary Table 4 of the paper.
Joint_LogAddI_TF4
Joint_LogAddI_CCU

## Explore bias on the logistic scale without interactions.
Joint_LogitN_TF4 <- summary(glm(sel1.01[sel.j0] ~ sex.01[sel.j0] + mother.age[sel.j0] + mother.educ[sel.j0] + mother.bmi[sel.j0] + mother.depr[sel.j0] + mother.smk[sel.j0] + mother.gest[sel.j0], family = binomial(link = "logit")))$coef[-1, -3]
rownames(Joint_LogitN_TF4) <- c("Child Sex", "Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")
Joint_LogitN_CCU <- summary(glm(sel2.01[sel.j0] ~ sex.01[sel.j0] + mother.age[sel.j0] + mother.educ[sel.j0] + mother.bmi[sel.j0] + mother.depr[sel.j0] + mother.smk[sel.j0] + mother.gest[sel.j0], family = binomial(link = "logit")))$coef[-1, -3]
rownames(Joint_LogitN_CCU) <- c("Child Sex", "Mother's Age", "Mother's Education", "Mother's BMI", "Mother's Depression", "Mother's Smoking", "Gestational Age")

## Report results - this is Supplementary Table 5 of the paper.
Joint_LogitN_TF4
Joint_LogitN_CCU

## Similar conclusions as for the marginal analyses.

## ---------- SAVE RESULTS ---------- ##

## Save the results tables.
save(Marg_Results, Marg_LogAddI_TF4, Marg_LogAddI_CCU, Bias_TF4, Bias_CCU, Marg_LogitN_TF4, Marg_LogitN_CCU,
     Miss_table, Joint_results, Joint_LogAddI_TF4, Joint_LogAddI_CCU, Joint_LogitN_TF4, Joint_LogitN_CCU,
     file = "ALSPAC_Results.RData")

##################################################
