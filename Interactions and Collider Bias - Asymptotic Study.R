
##########   INTERACTIONS AND COLLIDER BIAS   ##########

## This file contains R code in support of the paper "Relationship
## Between Collider Bias and Interactions on the Log-additive Scale".

## We first write an R function to perform the asymptotic study
## in the paper and then run the R function for each of the various
## scenarios we want to explore.

## Set up R.
setwd("...")
options(digits = 4)
library(xtable)

##################################################

##########  R FUNCTIONS   ##########

## Standard expit and logit functions.
expit <- function (x) exp(x) / (1 + exp(x))
logit <- function (x) log(x / (1 - x))

## ---------- SIMULATION GENERATOR ---------- ##

## Here is the main simulation function. This can be used 
## to run standard simulation studies, but here we want to 
## study asymptotic properties of collider bias so we will
## run an asymptotic study instead (iter = 1, large n).
interaction.sims <- function (iter, n, x.distr = "binary", y.model = "logistic", beta0, beta1, s.model = "logadd", sel.prob = 0.5, delta1, delta2, inter.range, latent.sigma = 1.6, n.d0 = 1e6, seed = 4189) {
  
  ## Store main results for the X-Y associations here.
  K <- length(inter.range)
  xy.beta <- matrix(0, iter, K); colnames(xy.beta) <- paste("I = ", inter.range, sep = "")
  xy.beta.se <- matrix(0, iter, K); colnames(xy.beta.se) <- paste("I = ", inter.range, sep = "")
  xy.int <- matrix(0, iter, K); colnames(xy.int) <- paste("I = ", inter.range, sep = "")
  xy.int.se <- matrix(0, iter, K); colnames(xy.int.se) <- paste("I = ", inter.range, sep = "")
  
  ## Store results from fitting selection models here.
  m1.delta <- matrix(0, iter, K); colnames(m1.delta) <- paste("I = ", inter.range, sep = "")
  m1.se <- matrix(0, iter, K); colnames(m1.se) <- paste("I = ", inter.range, sep = "")
  m2.delta <- matrix(0, iter, K); colnames(m2.delta) <- paste("I = ", inter.range, sep = "")
  m2.se <- matrix(0, iter, K); colnames(m2.se) <- paste("I = ", inter.range, sep = "")
  m3.delta <- matrix(0, iter, K); colnames(m3.delta) <- paste("I = ", inter.range, sep = "")
  m3.se <- matrix(0, iter, K); colnames(m3.se) <- paste("I = ", inter.range, sep = "")
  
  ## Store additional diagnostics here.
  min.s.pr <- matrix(0, iter, K); colnames(min.s.pr) <- paste("I = ", inter.range, sep = "")
  max.s.pr <- matrix(0, iter, K); colnames(max.s.pr) <- paste("I = ", inter.range, sep = "")
  prop.selected <- matrix(0, iter, K); colnames(prop.selected) <- paste("I = ", inter.range, sep = "")
  mean.y <- matrix(0, iter, K); colnames(mean.y) <- paste("I = ", inter.range, sep = "")
  cutoff <- matrix(0, iter, K); colnames(cutoff) <- paste("I = ", inter.range, sep = "")
  delta0 <- rep(0, K);  names(delta0) <- paste("I = ", inter.range, sep = "")
  if (s.model == "double") {
    rho1 <- rep(0, K);  names(rho1) <- paste("I = ", inter.range, sep = "")
    rho2 <- rep(0, K);  names(rho2) <- paste("I = ", inter.range, sep = "")
    prop.r1 <- matrix(0, iter, K); colnames(prop.r1) <- paste("I = ", inter.range, sep = "")
    prop.r2 <- matrix(0, iter, K); colnames(prop.r2) <- paste("I = ", inter.range, sep = "")
  }
  
  ## Start the loop across delta3 values.
  for (j in 1:K) {
    
    ## Define the seed.
    set.seed(seed + j * 10000)
    
    ## Specify the value of the intercept in the collider model.
    ## This is done by calling a separate function.
    d0 <- compute.d0(x.distr = x.distr, y.model = y.model, beta0 = beta0, 
                     beta1 = beta1, s.model = s.model, sel.prob = sel.prob,
                     delta1 = delta1, delta2 = delta2, delta3 = inter.range[j], 
                     latent.sigma = latent.sigma, n.d0 = n.d0)
    
    ## Store results.
    if (s.model == "double") {
      delta0[j] <- d0[1]
      rho1[j] <- d0[2]
      rho2[j] <- d0[3]
    } else {
      delta0[j] <- d0
    }
    
    ## Run the simulation.
    for (I in 1:iter) {
      
      ## ----- Simulate data. ------ ##
      
      ## Simulate exposure data.
      if (x.distr == "binary") {
        X <- rbinom(n, 1, 0.3) 
      } else if (x.distr == "unif") {
        X <- runif(n, 0, 1)
      } else if (x.distr == "normal") {
        X <- rnorm(n, 0, 1)
      }
      
      ## Simulate outcome data.
      if (y.model == "linear") {
        Y <- beta0 + beta1 * X + rnorm(n, 0, 0.5)
      } else if (y.model == "logistic") {
        Y.pr <- expit(beta0 + beta1 * X)
        Y <- rbinom(n, 1, prob = Y.pr)
      } else if (y.model == "poisson") {
        Y.lambda <- exp(beta0 + beta1 * X)
        Y <- rpois(n, lambda = Y.lambda)
      }
      
      ## Simulate data for the collider.
      if (s.model == "logadd") {
        
        ## First, for a log-additive collider.
        S.pr <- exp(d0 + delta1 * X + delta2 * Y + inter.range[j] * X * Y)
        cutoff[I, j] <- sum(S.pr > 1)
        S.pr[S.pr > 1] <- 1
        S <- rbinom(n, 1, prob = S.pr)
        
      } else if (s.model == "logistic") {
        
        ## Second, for a logistic collider.
        S.pr <- expit(d0 + delta1 * X + delta2 * Y + inter.range[j] * X * Y)
        S <- rbinom(n, 1, prob = S.pr)
        
      } else if (s.model == "probit") {
        
        ## Third, for a probit collider.
        S.pr <- d0 + delta1 * X + delta2 * Y + inter.range[j] * X * Y + rnorm(n, 0, latent.sigma)
        S <- as.numeric(S.pr > 0)
        
      } else if (s.model == "double") {
        
        ## And fourth, for a "double-threshold" collider.
        S.pr <- d0[1] + delta1 * X + delta2 * Y + inter.range[j] * X * Y + rnorm(n, 0, latent.sigma)
        S <- rep(0, n)
        S[S.pr < d0[2]] <- 1
        S[S.pr > d0[3]] <- 1
        
      }
      
      ## ----- Fit the models. ----- ##
      
      ## Fit the X-Y model for inference.
      if (y.model == "linear") {
        fit1 <- summary(lm(Y[S == 1] ~ X[S == 1]))$coefficients
      } else if (y.model == "logistic") {
        fit1 <- summary(glm(Y[S == 1] ~ X[S == 1], family = binomial(link = "logit")))$coefficients
      } else if (y.model == "poisson") {
        fit1 <- summary(glm(Y[S == 1] ~ X[S == 1], family = poisson(link = "log")))$coefficients
      }
      xy.beta[I, j] <- fit1[2, 1]
      xy.beta.se[I, j] <- fit1[2, 2]
      xy.int[I, j] <- fit1[1, 1]
      xy.int.se[I, j] <- fit1[1, 2]
      
      ## Fit a log-probability model for selection.
      fit2 <- summary(glm(S ~ X + Y + X * Y, family = poisson(link = "log")))$coefficients
      m1.delta[I, j] <- fit2[4, 1]
      m1.se[I, j] <- fit2[4, 2]
      
      ## Fit a logistic model for selection.
      fit2 <- summary(glm(S ~ X + Y + X * Y, family = binomial(link = "logit")))$coefficients
      m2.delta[I, j] <- fit2[4, 1]
      m2.se[I, j] <- fit2[4, 2]
      
      ## Fit a probit model for selection.
      fit2 <- summary(glm(S ~ X + Y + X * Y, family = binomial(link = "probit")))$coefficients
      m3.delta[I, j] <- fit2[4, 1]
      m3.se[I, j] <- fit2[4, 2]
      
      ## ----- Additional diagnostics. ----- ##
      
      ## Store additional diagnostics.
      min.s.pr[I, j] <- min(S.pr)
      max.s.pr[I, j] <- max(S.pr)
      prop.selected[I, j] <- mean(S)
      mean.y[I, j] <- mean(Y)
      
      ## Store further diagnostics for "double-threshold" models.
      if (s.model == "double") {
        prop.r1[I, j] <- sum(S.pr < d0[2])
        prop.r2[I, j] <- sum(S.pr > d0[3])
      }
      
    }
    
    ## Progress report.
    print(paste("Interaction value ", j, " done.", sep = ""))
    
  }
  
  ## Put everything together in a list.
  res.list <- list("beta.est" = xy.beta, "beta.se" = xy.beta.se, "int.est" = xy.int, "int.se" = xy.int.se, 
                   "logadd.est" = m1.delta, "logadd.se" = m1.se, "logit.est" = m2.delta, "logit.se" = m2.se, 
                   "probit.est" = m3.delta, "probit.se" = m3.se, "selection.prob" = prop.selected, 
                   "mean.y" = mean.y, "min.sel.probs" = min.s.pr, "max.sel.probs" = max.s.pr, "delta0" = delta0)
  if (s.model == "logadd") res.list$cutoff <- cutoff
  if (s.model == "double") res.list$rho1 <- rho1
  if (s.model == "double") res.list$rho2 <- rho2
  if (s.model == "double") res.list$prop.rho1 <- prop.r1
  if (s.model == "double") res.list$prop.rho2 <- prop.r2
  
  ## Goodbye.
  return(res.list)
  
}

## Auxiliary function that calculates the intercept delta0
## to yield a specified proportion of individuals selected.
## This can be done analytically for the log-additive
## collider model but only numerically for the other models.
compute.d0 <- function(x.distr, y.model, beta0, beta1, s.model, sel.prob, delta1, delta2, delta3, latent.sigma, n.d0) {
  
  ## Simulate exposure data.
  if (x.distr == "binary") {
    X <- rbinom(n.d0, 1, 0.3) 
  } else if (x.distr == "unif") {
    X <- runif(n.d0, 0, 1)
  } else if (x.distr == "normal") {
    X <- rnorm(n.d0, 0, 1)
  }
  
  ## Simulate outcome data.
  if (y.model == "linear") {
    Y <- beta0 + beta1 * X + rnorm(n.d0, 0, 0.5)
  } else if (y.model == "logistic") {
    Y.pr <- expit(beta0 + beta1 * X)
    Y <- rbinom(n.d0, 1, prob = Y.pr)
  } else if (y.model == "poisson") {
    Y.lambda <- exp(beta0 + beta1 * X)
    Y <- rpois(n.d0, lambda = Y.lambda)
  }
  
  ## Compute delta0.
  if (s.model == "logadd") {
    d0 <- log(sel.prob) - log( mean( exp(delta1 * X + delta2 * Y + delta3 * X * Y) ) )
  } else if (s.model == "logistic") {
    #damnit <- function (d0) sum((expit(d0 + delta1 * X + delta2 * Y + delta3 * X * Y) - sel.prob)^2)
    #d0 <- optimize(damnit, interval = c(-10, 10))$minimum
    find.d0 <- function (d0) mean(expit(d0 + delta1 * X + delta2 * Y + delta3 * X * Y)) - sel.prob
    d0 <- uniroot(find.d0, interval = c(-10, 10))$root
  } else if (s.model == "probit") {
    S.pr <- delta1 * X + delta2 * Y + delta3 * X * Y + rnorm(n.d0, 0, latent.sigma)
    d0 <- - sort(S.pr)[floor((1 - sel.prob) * n.d0)]
  } else if (s.model == "double") {
    d0 <- 0
    S.pr <- delta1 * X + delta2 * Y + delta3 * X * Y + rnorm(n.d0, 0, latent.sigma)
    r1 <- sort(S.pr)[floor(sel.prob / 2 * n.d0)]
    r2 <- sort(S.pr)[floor((1 - sel.prob / 2) * n.d0)]
  }
  
  ## Return.
  if (s.model == "double") {
    return (c(d0, r1, r2))
  } else {
    return (d0)
  }
}

## ---------- PLOTTING FUNCTIONS ---------- ##

## Plotting function: plot bias against interactions either
## on the scale of the fitted model or on the log-additive scale.
plot.sim.results <- function (sim, beta1, inter.range, x.range = NULL, y.range, scale = "modelled", main = "", abl = NULL, plot.type = "p", ci = TRUE, logadd.mark = FALSE) {
  
  ## How many runs?
  K <- length(inter.range)
  
  ## Define x-axis labels.
  if (is.null(x.range) & scale == "modelled") x.range <- inter.range
  
  ## Define the x values.
  if (scale == "modelled") {
    x.values <- inter.range
  } else if (scale == "logadd") {
    x.values <- colMeans(sim$logadd.est)
  }
  
  ## Define confidence interval means and limits.
  means <- colMeans(sim$beta.est) - beta1
  ci1 <- colMeans(sim$beta.est) - beta1 - qnorm(0.975) * colMeans(sim$beta.se)
  ci2 <- colMeans(sim$beta.est) - beta1 + qnorm(0.975) * colMeans(sim$beta.se)
  
  ## Set up the plotting colours.
  if (scale == "modelled") {
    cols <- rep("black", K)
    cols[which(inter.range == 0)] <- "red"
  } else if (scale == "logadd") {
    cols <- rep("black", K)
    if (is.numeric(logadd.mark)) cols[logadd.mark] <- "red"
    #cols[which.min(abs(colMeans(sim$logadd.est)))] <- "red"
  }
  
  ## Define the width of interval bars.
  ci.width <- (max(x.values) - min(x.values)) / 100
  
  ## Start plotting.
  plot(x = x.values, y = means, type = plot.type, axes = FALSE, main = main, xlab = "", ylab = "Bias", xlim = c(min(x.range), max(x.range)), ylim = c(min(y.range), max(y.range)), pch = 19, cex.main = 1, col = cols)
  if (ci == TRUE) {
    for (i in 1:K) {
      lines(y = c(ci1[i], ci2[i]), x = c(x.values[i], x.values[i]), col = cols[i])
      lines(y = c(ci1[i], ci1[i]), x = c(x.values[i] - ci.width, x.values[i] + ci.width), col = cols[i])
      lines(y = c(ci2[i], ci2[i]), x = c(x.values[i] - ci.width, x.values[i] + ci.width), col = cols[i])
    }
  }
  title(xlab = "Interaction", line = 2.4)
  
  ## Add the axes.
  axis(side = 1, at = x.range, labels = x.range, las = 1, cex.axis = 0.8, cex.lab = 1)
  axis(side = 2, at = y.range, labels = y.range, las = 1, cex.axis = 0.9, cex.lab = 1)
  abline(h = 0, lty = 2, lwd = 1, col = "brown")
  if (scale == "logadd") abline(v = 0, lty = 2, lwd = 1, col = "grey")
  
  ## Add a line that predicts the bias based on our theory.
  if (!(is.null(abl))) abline(a = abl[1], b = abl[2], lty = 2, col = "blue")
  
}

## Plot the results of multiple simulation runs simultaneously,
## to compare runs with different selection probabilities.
plot.many.results <- function (sim.list, beta1, inter.range, x.range = NULL, y.range, scale = "modelled", main = "", plot.type = "b", cols) {
  
  ## How many runs?
  K <- length(inter.range)
  
  ## How many simulations to compare?
  L <- length(sim.list)
  
  ## Define x-axis labels.
  if (is.null(x.range) & scale == "modelled") x.range <- inter.range
  
  ## Define the x values.
  if (scale == "modelled") {
    all.x.values <- t(matrix(inter.range, K, L))
  } else if (scale == "logadd") {
    all.x.values <- t(sapply(sim.list, function (x) colMeans(x$logadd.est)))
  }
  
  ## Define means per simulation run..
  all.means <- t(sapply(sim.list, function (x) colMeans(x$beta.est) - beta1))
  
  ## Start plotting.
  plot(x = all.x.values[1, ], y = all.means[1, ], type = plot.type, axes = FALSE, main = main, xlab = "", ylab = "Bias", xlim = c(min(x.range), max(x.range)), ylim = c(min(y.range), max(y.range)), pch = 19, cex.main = 1, col = cols[1])
  for (j in 2:L) {
    lines(x = all.x.values[j, ], y = all.means[j, ], type = plot.type, pch = 19, col = cols[j])
  }
  title(xlab = "Interaction", line = 2.4)
  
  ## Add the axes.
  axis(side = 1, at = x.range, labels = x.range, las = 1, cex.axis = 0.8, cex.lab = 1)
  axis(side = 2, at = y.range, labels = y.range, las = 1, cex.axis = 0.9, cex.lab = 1)
  abline(h = 0, lty = 2, lwd = 1, col = "brown")
  if (scale == "logadd") abline(v = 0, lty = 2, lwd = 1, col = "grey")
  
}

##################################################

##########   ASYMPTOTIC STUDY   ##########

## We use a sample size of 1e7 to approximate an infinite
## sample. We run each asymptotic experiment once.

## Each run takes about 20min in a standard desktop.

## ---------- DEFAULT (50% SELECTED) ---------- ##

## Log-additive collider model.
ASY1.0.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "logadd", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.2, 0.2, length = 11), seed = 3320184)
ASY1.0.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "logadd", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.1, sel.prob = 0.5, inter.range = seq(-0.2, 0.2, length = 11), seed = 3420184)
ASY1.0.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "logadd", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.05, sel.prob = 0.5, inter.range = seq(-0.1, 0.1, length = 11), seed = 3520184)

## Logistic collider model.
ASY1.1.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), seed = 3620184)
ASY1.1.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), seed = 3720184)
ASY1.1.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), seed = 3820184)

## Probit collider model.
ASY1.2.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 3920184)
ASY1.2.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 4020184)
ASY1.2.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 4120184)

## Double threshold collider model.
ASY1.3.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 4220184)
ASY1.3.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 4320184)
ASY1.3.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 4420184)

## ---------- 10% SELECTED ---------- ##

## Logistic collider model.
ASY4.1.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), seed = 6320184)
ASY4.1.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), seed = 6420184)
ASY4.1.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), seed = 6520184)

## Probit collider model.
ASY4.2.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 6620184)
ASY4.2.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 6720184)
ASY4.2.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 6820184)

## Double threshold collider model.
ASY4.3.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 6920184)
ASY4.3.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7020184)
ASY4.3.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.1, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7120184)

## ---------- 30% SELECTED ---------- ##

## Logistic collider model.
ASY5.1.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), seed = 7220184)
ASY5.1.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), seed = 7320184)
ASY5.1.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), seed = 7420184)

## Probit collider model.
ASY5.2.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7520184)
ASY5.2.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7620184)
ASY5.2.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7720184)

## Double threshold collider model.
ASY5.3.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7820184)
ASY5.3.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 7920184)
ASY5.3.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.3, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8020184)

## ---------- 70% SELECTED ---------- ##

## Logistic collider model.
ASY6.1.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), seed = 8120184)
ASY6.1.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), seed = 8220184)
ASY6.1.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), seed = 8320184)

## Probit collider model.
ASY6.2.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8420184)
ASY6.2.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8520184)
ASY6.2.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8620184)

## Double threshold collider model.
ASY6.3.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8720184)
ASY6.3.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8820184)
ASY6.3.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.7, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 8920184)

## ---------- 90% SELECTED ---------- ##

## Logistic collider model.
ASY7.1.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), seed = 9020184)
ASY7.1.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), seed = 9120184)
ASY7.1.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), seed = 9220184)

## Probit collider model.
ASY7.2.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 9320184)
ASY7.2.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 9420184)
ASY7.2.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 9520184)

## Double threshold collider model.
ASY7.3.1 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "logistic", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 9620184)
ASY7.3.2 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "linear", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 9720184)
ASY7.3.3 <- interaction.sims(1, 1e7, x.distr = "binary", y.model = "poisson", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.9, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 9820184)

## ---------- X ~ N(0, 1) ---------- ##

## Log-additive collider model.
ASY9.0.1 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "logistic", s.model = "logadd", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.2, 0.2, length = 11), seed = 11120184)
ASY9.0.2 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "linear", s.model = "logadd", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.1, sel.prob = 0.5, inter.range = seq(-0.2, 0.2, length = 11), seed = 11220184)
ASY9.0.3 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "poisson", s.model = "logadd", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.05, sel.prob = 0.5, inter.range = seq(-0.1, 0.1, length = 11), seed = 11320184)

## Logistic collider model.
ASY9.1.1 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "logistic", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), seed = 11420184)
ASY9.1.2 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "linear", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), seed = 11520184)
ASY9.1.3 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "poisson", s.model = "logistic", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), seed = 11620184)

## Probit collider model.
ASY9.2.1 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "logistic", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 11720184)
ASY9.2.2 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "linear", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 11820184)
ASY9.2.3 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "poisson", s.model = "probit", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 11920184)

## Double threshold collider model.
ASY9.3.1 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "logistic", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 12020184)
ASY9.3.2 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "linear", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 12120184)
ASY9.3.3 <- interaction.sims(1, 1e7, x.distr = "normal", y.model = "poisson", s.model = "double", beta0 = 0, beta1 = 0.2, delta1 = 0.3, delta2 = 0.3, sel.prob = 0.5, inter.range = seq(-0.5, 0.5, length = 11), latent.sigma = 1.6, seed = 12220184)

## ---------- SAVE RESULTS ---------- ##

## Save all asymptotic runs along with the functions.
save(ASY1.0.1, ASY1.0.2, ASY1.0.3, ASY1.1.1, ASY1.1.2, ASY1.1.3, 
     ASY1.2.1, ASY1.2.2, ASY1.2.3, ASY1.3.1, ASY1.3.2, ASY1.3.3, 
     ASY4.1.1, ASY4.1.2, ASY4.1.3, ASY4.2.1, ASY4.2.2, ASY4.2.3, 
     ASY4.3.1, ASY4.3.2, ASY4.3.3, ASY5.1.1, ASY5.1.2, ASY5.1.3, 
     ASY5.2.1, ASY5.2.2, ASY5.2.3, ASY5.3.1, ASY5.3.2, ASY5.3.3, 
     ASY6.1.1, ASY6.1.2, ASY6.1.3, ASY6.2.1, ASY6.2.2, ASY6.2.3, 
     ASY6.3.1, ASY6.3.2, ASY6.3.3, ASY7.1.1, ASY7.1.2, ASY7.1.3, 
     ASY7.2.1, ASY7.2.2, ASY7.2.3, ASY7.3.1, ASY7.3.2, ASY7.3.3, 
     ASY9.0.1, ASY9.0.2, ASY9.0.3, ASY9.1.1, ASY9.1.2, ASY9.1.3, 
     ASY9.2.1, ASY9.2.2, ASY9.2.3, ASY9.3.1, ASY9.3.2, ASY9.3.3, 
     expit, logit, compute.d0, interaction.sims, plot.sim.results,
     plot.many.results, file = "Asymptotic_results.RData")
#load("Asymptotic_results.RData")

##################################################

##########   SUMMARIZE RESULTS   ##########

## Plot results of the asymptotic study, create the paper's figures.

## We create 3x3 grids of plots where the columns correspond to a
## logistic, linear and Poisson regression model for the outcome 
## and the columns correspond to a logistic, probit and 
## "double-threshold" model for the collider.

## ---------- MAIN PAPER - FIGURE 2 ---------- ##

## Collider bias plotted against exposure-outcome interactions on
## the scale of the true selection model. 50% selected.

## Do the plot.
png("AsymptoticPlot1m.png", width = 900, height = 900, pointsize = 14, res = 108)
par(mfrow = c(3, 3))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plot.sim.results(ASY1.1.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "modelled", main = "Logistic Regression", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.1.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "modelled", main = "Linear Regression", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.1.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15), scale = "modelled", main = "Poisson Regression", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.2.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.2.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.2.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.3.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, 0, 0.1, 0.2), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.3.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.05, 0, 0.05, 0.1), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY1.3.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, 0, 0.1, 0.2), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
dev.off()

## ---------- MAIN PAPER - FIGURE 3 ---------- ##

## Collider bias plotted against exposure-outcome interactions on
## the log-additive scale. 50% selected.

## Do the plot.
png("AsymptoticPlot1l.png", width = 900, height = 900, pointsize = 14, res = 108)
par(mfrow = c(3, 3))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plot.sim.results(ASY1.1.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "logadd", main = "Logistic Regression", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY1.1.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "logadd", main = "Linear Regression", abl = c(0, 0.5^2), ci = FALSE)
plot.sim.results(ASY1.1.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15), y.range = c(-0.3, -0.15, 0, 0.15), scale = "logadd", main = "Poisson Regression", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY1.2.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY1.2.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "logadd", main = "", abl = c(0, 0.5^2), ci = FALSE)
plot.sim.results(ASY1.2.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15), y.range = c(-0.3, -0.15, 0, 0.15), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY1.3.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.05, 0, 0.05, 0.1, 0.15), y.range = c(-0.1, 0, 0.1, 0.2), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY1.3.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.05, 0, 0.05, 0.1, 0.15), y.range = c(-0.05, 0, 0.05, 0.1), scale = "logadd", main = "", abl = c(0, 0.5^2), ci = FALSE)
plot.sim.results(ASY1.3.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.1, 0, 0.1, 0.2), y.range = c(-0.1, 0, 0.1, 0.2), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
dev.off()

## ---------- MAIN PAPER - FIGURE 4 ---------- ##

## Collider bias plotted against exposure-outcome interactions on
## the scale of the true selection model. Various selection probabilities.

## Set up the plot.
png("AsymptoticPlot2m.png", width = 900, height = 900, pointsize = 14, res = 108)
par(mfrow = c(3, 3))
par(mar = c(4.1, 4.1, 2.1, 2.1))
cols <- c("darkgreen", "blue", "purple", "red", "orange")

## For a logistic collider model.
plot.many.results(sim.list = list(ASY4.1.1, ASY5.1.1, ASY1.1.1, ASY6.1.1, ASY7.1.1), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.5, -0.25, 0, 0.25, 0.5), scale = "modelled", cols = cols, main = "Logistic Regression")
plot.many.results(sim.list = list(ASY4.1.2, ASY5.1.2, ASY1.1.2, ASY6.1.2, ASY7.1.2), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.2, -0.1, 0, 0.1, 0.2), scale = "modelled", cols = cols, main = "Linear Regression")
plot.many.results(sim.list = list(ASY4.1.3, ASY5.1.3, ASY1.1.3, ASY6.1.3, ASY7.1.3), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.6, -0.3, 0, 0.3), scale = "modelled", cols = cols, main = "Poisson Regression")

## For a probit collider model.
plot.many.results(sim.list = list(ASY4.2.1, ASY5.2.1, ASY1.2.1, ASY6.2.1, ASY7.2.1), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.6, -0.3, 0, 0.3, 0.6), scale = "modelled", cols = cols)
plot.many.results(sim.list = list(ASY4.2.2, ASY5.2.2, ASY1.2.2, ASY6.2.2, ASY7.2.2), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.2, -0.1, 0, 0.1, 0.2), scale = "modelled", cols = cols)
plot.many.results(sim.list = list(ASY4.2.3, ASY5.2.3, ASY1.2.3, ASY6.2.3, ASY7.2.3), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.6, -0.3, 0, 0.3), scale = "modelled", cols = cols)

## For a double-threshold collider model.
plot.many.results(sim.list = list(ASY4.3.1, ASY5.3.1, ASY1.3.1, ASY6.3.1, ASY7.3.1), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.2, 0, 0.2, 0.4), scale = "modelled", cols = cols)
plot.many.results(sim.list = list(ASY4.3.2, ASY5.3.2, ASY1.3.2, ASY6.3.2, ASY7.3.2), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, 0, 0.1, 0.2), scale = "modelled", cols = cols)
plot.many.results(sim.list = list(ASY4.3.3, ASY5.3.3, ASY1.3.3, ASY6.3.3, ASY7.3.3), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, 0, 0.3, 0.6), scale = "modelled", cols = cols)

## Finish.
dev.off()

## ---------- MAIN PAPER - FIGURE 5 ---------- ##

## Collider bias plotted against exposure-outcome interactions on
## the log-additive scale. Various selection probabilities.

## Set up the plot.
png("AsymptoticPlot2l.png", width = 900, height = 900, pointsize = 14, res = 108)
par(mfrow = c(3, 3))
par(mar = c(4.1, 4.1, 2.1, 2.1))
cols <- c("darkgreen", "blue", "purple", "red", "orange")

## For a logistic collider model.
plot.many.results(sim.list = list(ASY4.1.1, ASY5.1.1, ASY1.1.1, ASY6.1.1, ASY7.1.1), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.5, -0.25, 0, 0.25, 0.5), y.range = c(-0.5, -0.25, 0, 0.25, 0.5), scale = "logadd", cols = cols, main = "Logistic Regression")
plot.many.results(sim.list = list(ASY4.1.2, ASY5.1.2, ASY1.1.2, ASY6.1.2, ASY7.1.2), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.5, -0.25, 0, 0.25, 0.5), y.range = c(-0.2, -0.1, 0, 0.1, 0.2), scale = "logadd", cols = cols, main = "Linear Regression")
plot.many.results(sim.list = list(ASY4.1.3, ASY5.1.3, ASY1.1.3, ASY6.1.3, ASY7.1.3), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.5, -0.25, 0, 0.25, 0.5), y.range = c(-0.6, -0.3, 0, 0.3), scale = "logadd", cols = cols, main = "Poisson Regression")

## For a probit collider model.
plot.many.results(sim.list = list(ASY4.2.1, ASY5.2.1, ASY1.2.1, ASY6.2.1, ASY7.2.1), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.6, -0.3, 0, 0.3, 0.6), y.range = c(-0.6, -0.3, 0, 0.3, 0.6), scale = "logadd", cols = cols)
plot.many.results(sim.list = list(ASY4.2.2, ASY5.2.2, ASY1.2.2, ASY6.2.2, ASY7.2.2), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.6, -0.3, 0, 0.3, 0.6), y.range = c(-0.2, -0.1, 0, 0.1, 0.2), scale = "logadd", cols = cols)
plot.many.results(sim.list = list(ASY4.2.3, ASY5.2.3, ASY1.2.3, ASY6.2.3, ASY7.2.3), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.6, -0.3, 0, 0.3), y.range = c(-0.6, -0.3, 0, 0.3), scale = "logadd", cols = cols)

## For a double-threshold collider model.
plot.many.results(sim.list = list(ASY4.3.1, ASY5.3.1, ASY1.3.1, ASY6.3.1, ASY7.3.1), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.2, 0, 0.2, 0.4, 0.6), y.range = c(-0.2, 0, 0.2, 0.4, 0.6), scale = "logadd", cols = cols)
plot.many.results(sim.list = list(ASY4.3.2, ASY5.3.2, ASY1.3.2, ASY6.3.2, ASY7.3.2), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.2, 0, 0.2, 0.4), y.range = c(-0.05, 0, 0.05, 0.1, 0.15), scale = "logadd", cols = cols)
plot.many.results(sim.list = list(ASY4.3.3, ASY5.3.3, ASY1.3.3, ASY6.3.3, ASY7.3.3), beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range <- c(-0.2, 0, 0.2, 0.4, 0.6), y.range = c(-0.2, 0, 0.2, 0.4, 0.6), scale = "logadd", cols = cols)

## Finish.
dev.off()

## ---------- SUPPLEMENTARY TABLE 1 ---------- ##

## We report numerical results for our "default" simulation.

## In this table, we report interaction real values, collider 
## bias figures and interaction estimates in the  log-additive
## model for our asymptotic study (with 50% selected).

## Create the table.
ASY_TABLE <- matrix(0, 33, 7)
colnames(ASY_TABLE) <- c("True Delta3", "Logit Bias", "Logit Inter", "Linear Bias", "Linear Inter", "Poisson Bias", "Poisson Inter")

## Set true interaction parameter values.
ASY_TABLE[, 1] <- rep(seq(-0.5, 0.5, length = 11), times = 3)

## Set bias and interaction values for a logistic collider.
ASY_TABLE[1:11, 2] <- ASY1.1.1$beta.est - 0.2
ASY_TABLE[1:11, 3] <- ASY1.1.1$logadd.est
ASY_TABLE[1:11, 4] <- ASY1.1.2$beta.est - 0.2
ASY_TABLE[1:11, 5] <- ASY1.1.2$logadd.est * 0.5^2
ASY_TABLE[1:11, 6] <- ASY1.1.3$beta.est - 0.2
ASY_TABLE[1:11, 7] <- ASY1.1.3$logadd.est

## Set bias and interaction values for a probit collider.
ASY_TABLE[12:22, 2] <- ASY1.2.1$beta.est - 0.2
ASY_TABLE[12:22, 3] <- ASY1.2.1$logadd.est
ASY_TABLE[12:22, 4] <- ASY1.2.2$beta.est - 0.2
ASY_TABLE[12:22, 5] <- ASY1.2.2$logadd.est * 0.5^2
ASY_TABLE[12:22, 6] <- ASY1.2.3$beta.est - 0.2
ASY_TABLE[12:22, 7] <- ASY1.2.3$logadd.est

## Set bias and interaction values for a double-threshold collider.
ASY_TABLE[23:33, 2] <- ASY1.3.1$beta.est - 0.2
ASY_TABLE[23:33, 3] <- ASY1.3.1$logadd.est
ASY_TABLE[23:33, 4] <- ASY1.3.2$beta.est - 0.2
ASY_TABLE[23:33, 5] <- ASY1.3.2$logadd.est * 0.5^2
ASY_TABLE[23:33, 6] <- ASY1.3.3$beta.est - 0.2
ASY_TABLE[23:33, 7] <- ASY1.3.3$logadd.est

## Display the table.
ASY_TABLE
xtable(ASY_TABLE, digits = c(0, 1, rep(3, 6)))

## ---------- SUPPLEMENTARY FIGURE 1 ---------- ##

## Collider bias plotted against exposure-outcome interactions on
## the scale of the true selection model. X ~ N(0, 1). 50% selected.

## Do the plot.
png("AsymptoticPlot4m.png", width = 900, height = 900, pointsize = 14, res = 108)
par(mfrow = c(3, 3))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plot.sim.results(ASY9.1.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "modelled", main = "Logistic Regression", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.1.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "modelled", main = "Linear Regression", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.1.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15), scale = "modelled", main = "Poisson Regression", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.2.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.2.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.2.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.3, -0.15, 0, 0.15), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.3.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.03, 0, 0.03, 0.06), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.3.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.01, 0, 0.01, 0.02), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
plot.sim.results(ASY9.3.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), y.range = c(-0.05, 0, 0.05, 0.1), scale = "modelled", main = "", ci = FALSE, plot.type = "b")
dev.off()

## ---------- SUPPLEMENTARY FIGURE 2 ---------- ##

## Collider bias plotted against exposure-outcome interactions on
## the log-additive scale. X ~ N(0, 1). 50% selected.

## Do the plot.
png("AsymptoticPlot4l.png", width = 900, height = 900, pointsize = 14, res = 108)
par(mfrow = c(3, 3))
par(mar = c(4.1, 4.1, 2.1, 2.1))
plot.sim.results(ASY9.1.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "logadd", main = "Logistic Regression", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY9.1.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "logadd", main = "Linear Regression", abl = c(0, 0.5^2), ci = FALSE)
plot.sim.results(ASY9.1.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15), y.range = c(-0.3, -0.15, 0, 0.15), scale = "logadd", main = "Poisson Regression", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY9.2.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.3, -0.15, 0, 0.15, 0.3), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY9.2.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.15, 0, 0.15, 0.3), y.range = c(-0.1, -0.05, 0, 0.05, 0.1), scale = "logadd", main = "", abl = c(0, 0.5^2), ci = FALSE)
plot.sim.results(ASY9.2.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.3, -0.2, -0.1, 0, 0.1), y.range = c(-0.3, -0.15, 0, 0.15), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY9.3.1, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.04, 0, 0.04, 0.08), y.range = c(-0.03, 0, 0.03, 0.06), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
plot.sim.results(ASY9.3.2, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(0.04, 0.05, 0.06, 0.07, 0.08), y.range = c(-0.01, 0, 0.01, 0.02), scale = "logadd", main = "", abl = c(0, 0.5^2), ci = FALSE)
plot.sim.results(ASY9.3.3, beta1 = 0.2, inter.range = seq(-0.5, 0.5, length = 11), x.range = c(-0.05, 0, 0.05, 0.1), y.range = c(-0.05, 0, 0.05, 0.1), scale = "logadd", main = "", abl = c(0, 1), ci = FALSE)
dev.off()

##################################################
