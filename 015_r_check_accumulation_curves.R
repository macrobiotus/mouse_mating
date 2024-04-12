#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis testing"
#' author: "Paul Czechowski `paul.czechowski@helmholtz-munich.de`"
#' date: "`r Sys.Date()`"
#' output:
#'  html_notebook:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: show
#'  pdf_document:
#'     toc: true
#'     number_sections: true
#'  html_document:
#'     toc: true
#'     toc_float: true
#'     number_sections: true
#'     code_folding: show
#' ---


# Prepare environment ---- 


# _1.) Collect garbage ----

rm(list=ls())
gc()

# _2.) Packages ----

# Nonlinear Mixed-Effects Growth Models: A Tutorial Using 'saemix' in R
# Methodology, 2021, Vol. 17(4), 250–270, https://doi.org/10.5964/meth.7061

library("here")
library("saemix")

# _3.) Research questions ----

# RQ1: Does Sex have an association with the total growth, rate of approach to
# the upper asymptotic, or point of inflection (of the curve found optimal in
# answering the first research question?)

# RQ2: Does adding SES to the model in addition to Sex as a predictor of total
# growth, rate of approach to the upper asymptotic, or point of inflection
# improve model fit? If so, how do SES and Sex relate to achievement growth?

# RQ3: Does Sex moderate the association between SES and total growth, rate of
# approach to the upper asymptotic, or point of inflection?

# Setup data and model ----

# _1.) Read in data ----

NLMEGMExData <- read.csv(here("../communication/Boedeker_2021_Nonlinear_mixed_effects_growth_models_SUPPL_Dataset_Code/7061_Dataset_for_replicating_example_analyses.csv"))

# _2.) Define Gompertz function ----

gompertz.model <- function(psi, id, x) { # psi, id, and x components are passed in from data - data frame, subject variable, time variable
  
  t <- x[ , 1] 
  LwrAsy   <- psi[id, 4]   # (d) lower asymptote for subject - where growth begins on y axis
  Apprch   <- psi[id, 2]   # (b) approach rate to upper asymtote - larger values indicate quicker growth - curve steepness - where growth begins on x axis
  Timing   <- psi[id, 3]   # (c) inflection point - time with greatest change - when large accelerated growth toward the upper limit occurs later
  TtlGrwth <- psi[id, 1]   # (a) total change for subject over time - large values indicate greater total growth -  range of y values
  
  ypred <- LwrAsy + TtlGrwth * exp(-exp(-Apprch * (t - Timing)))
  
  return(ypred)
  
}

# _3.) Set modelling options ----

NLMEGM.options <- list(seed = 1234, displayProgress = FALSE)


# RQ1 Which of the Logistic, Gompertz, or Richards curves models best the growth in achievement? ----

# _1.) Define data ----

GompertzData.RQ1 <- saemixData(
  name.data = NLMEGMExData, header = TRUE, name.group = c("ID"), name.predictors = c("time"), name.response = c("Achievement"), name.X = "time"
)

# _2.) Define model object ----

GompertzModel.RQ1 <- saemixModel(model = gompertz.model, 
                                 description = 'Gompertz', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), # starting values for each of the four growth parameters
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE), # which of the four parameter shoul be estimated? All but the last here
                                 transform.par = c(0, 0, 0, 0) # distribution of each of the four parameter
)

# _3.) Fit model ----

GompertzFit.RQ1 <- saemix(GompertzModel.RQ1, GompertzData.RQ1, NLMEGM.options)

# Nonlinear mixed-effects model fit by the SAEM algorithm

# ----         Data and Model          
# Data
# Dataset NLMEGMExData 
# Longitudinal data: Achievement ~ time | ID 
# 
# [...]
#
# ----                  Results                   ---
#
# Fixed effects
# Parameter Estimate   SE     CV(%)  # estimates in units of measurement per time
# TtlGrwth   5.0692  0.06830  1.35 
# Apprch     0.9489  0.04099  4.32 
# Timing     1.8065  0.04571  2.53 
# LwrAsy    -0.0915  0.02691 29.42   # variation is high as value is constrained?
# a.1        0.2352  0.00932  3.96 
# 
# Variance of random effects
# Parameter           Estimate   SE    CV(%) # large veraition indicates additional factors may contribute to this
# omega2.TtlGrwth      0.18546 0.0439   23.7
# omega2.Apprch        0.12201 0.0223   18.3
# omega2.Timing        0.16876 0.0276   16.4
# cov.TtlGrwth.Apprch  0.02935 0.0232   78.9
# cov.TtlGrwth.Timing -0.00556 0.0246 -441.8
# cov.Apprch.Timing   -0.00837 0.0179 -213.8
# 
# Correlation matrix of random effects   # check for insights into the relationships between these paramter
# omega2.TtlGrwth omega2.Apprch omega2.Timing
# omega2.TtlGrwth  1.0000          0.1951       -0.0315      
# omega2.Apprch    0.1951          1.0000       -0.0583      
# omega2.Timing   -0.0315         -0.0583        1.0000      
# 
# Statistical criteria
# Likelihood computed by linearisation
# -2LL= 623.3004 
# AIC= 645.3004 
# BIC= 673.9573 
# Likelihood computed by importance sampling
# -2LL= 622.5147 
# AIC= 644.5147 
# BIC= 673.1715

# RQ2: Does Sex have an association with the total growth, rate of approach to the upper asymptotic, or point of inflection (of the curve found optimal in answering the first research question?) ----

# _1.) Define data ----

GompertzData.RQ2 <- saemixData(
  name.data = NLMEGMExData, header = TRUE, name.group = c("ID"), name.predictors = c("time"), name.response = c("Achievement"), name.X = "time", 
  name.covariates = c("Sex")
)

# _2.) Define model object ----

GompertzModel.RQ2 <- saemixModel(model = gompertz.model, 
                                 description = 'Gompertz', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), # starting values for each of the four growth parameters
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE), # which of the four parameter shoul be estimated? All but the last here
                                 covariate.model = matrix(c(1,1,1,0), ncol = 4, byrow = TRUE), # which parameters are influenced by the covariate?
                                 transform.par = c(0, 0, 0, 0) # distribution of each of the four parameter
)

# _3.) Fit model ----

GompertzFit.RQ2 <- saemix(GompertzModel.RQ2, GompertzData.RQ2, NLMEGM.options)


# RQ3: Does adding SES (socioeconomic status) to the model in addition to Sex as a predictor of total growth, rate of approach to the upper asymptotic, or point of inflection improve model fit? If so, how do SES and Sex relate to achievement growth? ----

# _1.) Define data ----

GompertzData.RQ3 <- saemixData(
  name.data = NLMEGMExData, header = TRUE, name.group = c("ID"), name.predictors = c("time"), name.response = c("Achievement"), name.X = "time", 
  name.covariates = c("Sex", "SES")
)

# _2.) Define model object ----

GompertzModel.RQ3 <- saemixModel(model = gompertz.model, 
                                 description = 'Gompertz', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), # starting values for each of the four growth parameters
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE), # which of the four parameter shoul be estimated? All but the last here
                                 covariate.model = matrix(c(1,1,1,0, 1,1,1,0), ncol = 4, byrow = TRUE), # which parameters are influenced by the covariates? Length n paramters x covariates
                                 transform.par = c(0, 0, 0, 0) # distribution of each of the four parameter
)

# _3.) Fit model ----

GompertzFit.RQ3 <- saemix(GompertzModel.RQ3, GompertzData.RQ3, NLMEGM.options)


# RQ4: Does Sex moderate the association between SES and total growth, rate of approach to the upper asymptotic, or point of inflection? ----

# A test of moderation allows the investigator to determine if the relationship
# between a predictor and outcome is dependent on the value of a third variable.
# Moderation analyses can be used to determine for whom relationships hold or
# treatments are most effective. In growth modeling, moderation can be used to
# determine if, for example, the positive association between SES and the Total
# Growth parameter is equivalent for males and females or if for one sex the
# relationship is stronger.

# _1.) Define data ----

GompertzData.RQ4 <- saemixData(
  name.data = NLMEGMExData, header = TRUE, name.group = c("ID"), name.predictors = c("time"), name.response = c("Achievement"), name.X = "time", 
  name.covariates = c("Sex", "SES", "SexSESmod") # third covariate added - Conditional Growth / Moderation / Interaction
)

# _2.) Define model object ----

GompertzModel.RQ4 <- saemixModel(model = gompertz.model, 
                                 description = 'Gompertz', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), # starting values for each of the four growth parameters
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE), # which of the four parameter shoul be estimated? All but the last here
                                 covariate.model = matrix(c(1,1,1,0, 1,1,1,0, 1,1,1,0), ncol = 4, byrow = TRUE), # which parameters are influenced by the covariates? Length n paramters x covariates
                                 transform.par = c(0, 0, 0, 0) # distribution of each of the four parameter
)

# _3.) Fit model ----

GompertzFit.RQ4 <- saemix(GompertzModel.RQ4, GompertzData.RQ4, NLMEGM.options)

summary(GompertzFit.RQ4)






