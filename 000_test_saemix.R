
# 2024-4-11 - Paul Czechowski

# Nonlinear Mixed-Effects Growth Models: A Tutorial Using 'saemix' in R
# Methodology, 2021, Vol. 17(4), 250–270, https://doi.org/10.5964/meth.7061

library("here")
library("saemix")

# Research questions ----

# RQ1: Which of the Logistic, Gompertz, or Richards curves models best the
# growth in achievement?

# RQ2: Does Sex have an association with the total growth, rate of approach to
# the upper asymptotic, or point of inflection (of the curve found optimal in
# answering the first research question?)

# RQ3: Does adding SES to the model in addition to Sex as a predictor of total
# growth, rate of approach to the upper asymptotic, or point of inflection
# improve model fit? If so, how do SES and Sex relate to achievement growth?

# RQ4: Does Sex moderate the association between SES and total growth, rate of
# approach to the upper asymptotic, or point of inflection?

# Get example data ----

NLMEGMExData <- read.csv(here("../communication/Boedeker_2021_Nonlinear_mixed_effects_growth_models_SUPPL_Dataset_Code/7061_Dataset_for_replicating_example_analyses.csv"))

# Define Gompertz function ---

gompertz.model <- function(psi, id, x) { # psi, id, and x components are passed in from data - data frame, subject variable, time variable
  
  t <- x[ , 1] 
  LwrAsy   <- psi[id, 4]   # (d) lower asymptote for subject - where growth begins on y axis
  Apprch   <- psi[id, 2]   # (b) approach rate to upper asymtote - larger values indicate quicker growth - curve steepness - where growth begins on x axis
  Timing   <- psi[id, 3]   # (c) inflection point - time with greatest change - when large accelerated growth toward the upper limit occurs later
  TtlGrwth <- psi[id, 1]   # (a) total change for subject over time - large values indicate greater total growth -  range of y values
  
  ypred <- LwrAsy + TtlGrwth * exp(-exp(-Apprch * (t - Timing)))
  
  return(ypred)
  
}

# Set modlling options ----

NLMEGM.options <- list(seed = 1234, displayProgress = FALSE)


# RQ1 Which of the Logistic, Gompertz, or Richards curves models best the growth in achievement? ----

# Defeine data
GompertzData.RQ1 <- saemixData(
  name.data = NLMEGMExData, header = TRUE, name.group = c("ID"), name.predictors = c("time"), name.response = c("Achievement"), name.X = "time"
  )

