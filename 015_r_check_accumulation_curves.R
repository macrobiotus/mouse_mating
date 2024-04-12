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
library("npde")

# _3.) Revise when script is writte, or erase: Research questions ----

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

# Remove this test date in the line below upon script completion
NLMEGMExData <- read.csv(here("../communication/Boedeker_2021_Nonlinear_mixed_effects_growth_models_SUPPL_Dataset_Code/7061_Dataset_for_replicating_example_analyses.csv"))

mice_f0_slct <- readRDS( file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS( file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Define possible function ----


# __a) Gompertz ----

gompertz.model <- function(psi, id, x) { # psi, id, and x components are passed in from data - data frame, subject variable, time variable
  
  # Gompertz parameters - see https://doi.org/10.5964/meth.7061
  t <- x[ , 1] 
  LwrAsy   <- psi[id, 4]   # (d) lower asymptote for subject - where growth begins on y axis
  Apprch   <- psi[id, 2]   # (b) approach rate to upper asymtote - larger values indicate quicker growth - curve steepness - where growth begins on x axis
  Timing   <- psi[id, 3]   # (c) inflection point - time with greatest change - when large accelerated growth toward the upper limit occurs later
  TtlGrwth <- psi[id, 1]   # (a) total change for subject over time - large values indicate greater total growth -  range of y values
  
  # Gompertz curve formula
  ypred <- LwrAsy + TtlGrwth * exp(-exp(-Apprch * (t - Timing)))
  
  return(ypred)
  
}

# __b) Logistic ----

logistic.model <- function(psi, id, x) { # psi, id, and x components are passed in from data - data frame, subject variable, time variable
  
  # Logistic parameters - see https://doi.org/10.5964/meth.7061
  t <- x[ ,1]
  LwrAsy <- psi[id, 4]
  Apprch <- psi[id, 2]
  Timing <- psi[id, 3]
  TtlGrwth <- psi[id, 1]
  
  # Logistic curve formula
  ypred  <- LwrAsy + TtlGrwth/ (1+exp ( -Apprch * (t-Timing)))
  
  return(ypred)
}


# _3.) Set modelling options ----

NLMEGM.options <- list(seed = 1234, displayProgress = FALSE)



# RQ1 Use Gompertz or logistic curve to model weight gain? ----

# _1.) Define data ----

# __a) Gompertz ----

GompertzData.RQ1 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
)


# __b) Logistic ----

LogisticData.RQ1 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )

# _2.) Define model object ----

# __a) Gompertz ----

GompertzModel.RQ1 <- saemixModel(model = gompertz.model, 
                                 description = 'Gompertz', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), # starting values for each of the four growth parameters
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE), # which of the four parameter shoul be estimated? All but the last here
                                 transform.par = c(0, 0, 0, 0) # distribution of each of the four parameter
                                 )

# __b) Logistic ----

LogisticModel.RQ1 <- saemixModel(model = logistic.model,
                                 description = 'Logistic', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), 
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE),
                                 transform.par=c(0,0,0,0)
                                 )


# _3.) Fit model ----

# __a) Gompertz ----

GompertzFit.RQ1 <- saemix(GompertzModel.RQ1, GompertzData.RQ1, NLMEGM.options)

# Parameter Estimate   SE     CV(%)
# TtlGrwth  16.4638  1.16170  7.06  [16.4 g over time] 
# Apprch     0.0575  0.00279  4.86  [steepness reference value - modest]
# Timing    27.6851  1.81702  6.56  [growth inflection time at 27 days ]
# LwrAsy     8.4239  1.05141 12.48  [not estimated 8.4 g birth weight]
# a.1        0.8498  0.04786  5.63 

TtlGrwth <- psi[id, 1]   # (a) total change for subject over time - large values indicate greater total growth -  range of y values
Apprch   <- psi[id, 2]   # (b) approach rate to upper asymtote - larger values indicate quicker growth - curve steepness - where growth begins on x axis
Timing   <- psi[id, 3]   # (c) inflection point - time with greatest change - when large accelerated growth toward the upper limit occurs later
LwrAsy   <- psi[id, 4]   # (d) lower asymptote for subject - where growth begins on y axis





# __b) Logistic ----

LogisticFIT.RQ1 <- saemix(LogisticModel.RQ1, LogisticData.RQ1, NLMEGM.options)


# _4.) Plot model

# Compute normalised prediction distribution errors (npde) and normalized prediction discrepancies (npd).

# __a) Gompertz ----

plot(GompertzFit.RQ1, plot.type="observations.vs.predictions" )
plot(GompertzFit.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.GompertzFit.RQ1 <- npdeSaemix(GompertzFit.RQ1) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# __b) Logistic ----

plot(LogisticFIT.RQ1, plot.type = "observations.vs.predictions")
plot(LogisticFIT.RQ1, plot.type="observations.vs.predictions" )
plot(LogisticFIT.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.LogisticFIT.RQ1 <- npdeSaemix(LogisticFIT.RQ1)

# Continue here - define effcet on sex

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






