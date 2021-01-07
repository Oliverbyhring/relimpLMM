
#' A function that calculates the explained variance, R^2, of a random intercept model with one random effect.
#' The method is described an article by Nakagawa and Schielzeth(2013);
#' A general and simple method for obtaining R^2 from generalized linear mixed-effects models.
#' @param LMM The lmer object you want to calculate the explained variance of
#' @param marginal.r2 whether the marginal(only regard fixed effect variance as explained variance)
#' or the conditional(also regards the random effect variance as explained) R^2 should be computed, defaults to TRUE which computes the marginal
#' @return A numeric value. The estimated explained variance of the input model.
#' @keywords explained variance
#' @export

calc.R2 <- function(LMM, marginal.r2 = TRUE){
  random.intercept <- as.data.frame(lme4::VarCorr(LMM)) # extracting the random intercept and residual standard errors
  # Extracting the random intercept standard and residual error
  intercept.SD <- random.intercept$sdcor[1]
  residual.SD <- random.intercept$sdcor[2]
  # Calculating the random intercept and residual variance
  intercept.variance <- intercept.SD^2
  residual.variance <- residual.SD^2
  # Exctracting the beta estimates from the fitted model
  betas <- LMM@pp$delb[-1]
  betas.squared <- betas^2
  # Extracting the covariates used to fit the model
  LMM.data <- data.frame(LMM@pp$X)
  # Creating a list of all the variances of the covariates
  regressor.variances <- sapply(LMM.data[-1], FUN =  stats::var)
  # Creating a list of all the correlations of the covariates
  regressor.correlations <- stats::cor(LMM.data[-1])
  regressor.covariances <- 0
  n <-  length(LMM.data[-1])
  # Calculating the sum of the covariances
  if(n>1){
    for (j in 1:(n-1)){
      for (k in (j+1):n){
        regressor.covariances <- regressor.covariances
        + betas[j]*betas[k]*sqrt(regressor.variances[j] * regressor.variances[k]) * regressor.correlations[j,k]
      }
    }
  }
  if(length(regressor.variances)==0){
    betas.squared <- 0
    regressor.variances <- 0
  }
  # Calculating the variance of the fixed effects using Equation 2 (groemping 2007 Estimators of Relative Importance in Linear Regression Based on
  # Variance Decomposition)
  variance.fixed.effects <- betas.squared%*%regressor.variances + 2 * regressor.covariances
  if (marginal.r2 == T){
    # Only the fixed effect variance is regarded as explained
    explained.variance <- variance.fixed.effects
  }else{
    # Both the fixed effect variance and the random intercept variance is regarded as explained
    explained.variance <- variance.fixed.effects+intercept.variance
  }
  total.variance <- variance.fixed.effects+intercept.variance+residual.variance
  # R2 is the proportion of the variance that is explained
  R2 <- explained.variance/total.variance
  return(R2)
}

#' A function creates a matrix of orderings for which parameters that are to be included in the calculation of the average increase in R^2
#' The rows of the matrix denotes which parameters that are to be included in the model
#' @param n number of parameters that is going to be permuted
#' @return A matrix

create.ordering <- function(n){
  ordering <- data.frame(FrF2::FrF2(nruns = 2^n, nfactors = n, randomize = F))
  ordering <- sapply(ordering, as.numeric)
  ordering <- ordering-1
  return(ordering)
}

#' A function creates a subset of the full df based on a logical list.
#' @param ordering A binary list
#' @param df The data.frame to create a subset of
#' @return A subset data.frame of df

get.subset <- function(ordering, df){
  boolean.list <- as.logical(ordering)
  return(df[boolean.list])
}


#' A function that calculates the LMG-contribution given a subset model. The contribution is calculated as discribed by Byhring(2019);
#' Relative variable importance in linear regression models with random intercept term
#' @param df A data.frame containing the fixed effect variables
#' @param var1 The response variable
#' @param var2 The random intercept variable
#' @param var3 The variable of interest, that is the one we are interested in calculating the importance of.
#' @param var4 The number of fixed effects in the full model.
#' @return A numeric value containing the LMG contribution of a subset model.


calc.LMG.contribution <- function(df, var1, var2, var3, var4){
  model.df <- data.frame(response.variable = unlist(var1), random.intercept.var = unlist(var2), variable.of.interest = unlist(var3), df)
  model1 <- lme4::lmer(response.variable ~ . - random.intercept.var - variable.of.interest + (1|random.intercept.var), data = model.df)
  model2 <- lme4::lmer(response.variable ~ . - random.intercept.var + (1|random.intercept.var), data = model.df)
  LMG.contribution <- calc.R2(model2,T) - calc.R2(model1,T)
  n.S <- length(names(df))
  n.S.factorial <- factorial(n.S) # number of permutations possible for the variables that appears before the variable of interest
  p.minus.n.S.minus.one.factorial <- factorial(var4-n.S-1) # number of permutations possible for the variables that appears after the variable of interest
  n.equal.models <-  n.S.factorial*p.minus.n.S.minus.one.factorial
  tot.models <- factorial(var4)
  avg.contribution.factor <- n.equal.models/tot.models
  return(avg.contribution.factor * LMG.contribution)
}



#' A function that calculate the relative importance of all fixed effects and random effects in random intercept models created with the lme4 package.
#' It is currently limited to one random effect. The method is an extension of one proposed by Lindeman, Marenda and Gold(1980).
#' The method is described and reviewd by Groemping(2007);
#' A general and simple method for obtaining R^2 from generalized linear mixed-effects models
#' @param lmm.obj An lmer object with one random intercept term
#' @param response.name The name of the response variable as a string
#' @return A numeric value containing the LMG contribution of a subset model.
#' @export
#' @keywords variable importance

calc.relimp.lmm <- function(lmm.obj, response.name){
  model.summary <- summary(lmm.obj) # extracts the summary of the lmm object

  # Exctract the response and random intercept
  response.name <-  response.name # extracts the name of the response variable
  response.data <- data.frame(lmm.obj@frame[response.name]) # extracts the data of the response variable
  random.intercept.name <- names(model.summary$ngrps) # extracts the name of the random intercept variable
  random.intercept.data <- data.frame(lmm.obj@frame[random.intercept.name]) # extracts the data of the random intercept variable


  fixef.names <- names(lme4::fixef(lmm.obj)[-1]) # extracts the name of the fixed effects variables
  ordering <- create.ordering(length(fixef.names)-1) # creates an ordering of which variables that should be included in the set S
  ordering.list <- split(ordering, seq(nrow(ordering))) # making the orderings into seperate lists for future operations

  p <-length(fixef.names) # p is the number of fixed effects that is going to be permuted
  p.factorial <- factorial(p) # total number of submodels

  importances <- c() # empty list of importances
  #looping through the fixed effect variables
  for (i in fixef.names){
    variable.of.interest.data <- lmm.obj@frame[i] # extract the column that contain the data of the fixed effect of interest.
    fixef.df <- lmm.obj@frame[fixef.names[fixef.names != i]] # extract a df that contain the data of all the fixed effects except the one of interest
    subsets.of.full.df <- lapply(ordering.list, get.subset, df = fixef.df) # Creates all the subsets S, these are stored in a large list
    contributions <- lapply(X = subsets.of.full.df, FUN = calc.LMG.contribution, var1 = response.data, var2 = random.intercept.data, var3 = variable.of.interest.data, var4 = p)
    importances <- c(importances, sum(unlist(contributions))) # adding the importance of the variable of interest
  }
  ri.importance <- calc.R2(lmm.obj,F)-calc.R2(lmm.obj,T)
  importances <- c(importances,ri.importance)
  importances <- as.data.frame(importances, row.names(c(fixef.names,random.intercept.name)))
  return(importances)
}
