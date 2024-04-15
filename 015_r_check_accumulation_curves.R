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

# Setup data and model ----

# _1.) Read in data ----

mice_f0_slct <- readRDS( file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS( file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Check data ----

mice_f1_slct %>% dplyr::select(MeasurementDay, BodyWeightG, AnimalId, AnimalSex) %>% arrange(AnimalSex, AnimalId, MeasurementDay) %>% print(n = Inf)

ggplot(data = mice_f1_slct, aes(x = "Week", y="BodyWeightG",  group = "AnimalId")) +
  geom_line( mapping = aes( x = Week, y = BodyWeightG, group = AnimalId, color = AnimalId)) + 
  facet_grid(. ~ AnimalSex) +
  theme_bw()

# _3.) Define possibly applicable model functions ----

# __a) Gompertz

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

# __b) Logistic

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

# __c) Exponential decay

decay.model <- function(psi, id, xidep) { 
  
  # input: 
  #    psi : matrix of parameters (3 columns, ka, V, CL) 
  #     id : vector of indices 
  #  xidep : dependent variables (same nb of rows as length of id) 
  # returns: 
  #   a vector of predictions of length equal to length of id 
  x <- xidep[ , 1]  
  a <- psi[id, 1]   # upper asymptote
  b <- psi[id, 2]   # starting value
  k <- psi[id, 3]   # approach steepness, a will be reached later
  
  # Exponential Decay (increasing form) - 
  # https://people.richland.edu/james/lecture/m116/logs/models.html
  
  ypred <- a * (1 - b * exp( -k * x))
  
  return(ypred) 
  }

# testing function
# a = 10
# b = 0.5
# k = 0.04
# curve(a * (1 - b * exp( -k * x)), from = 0, to = 150, xlab="x", ylab="y", add = TRUE)

# _3.) Set modelling options ----

saemix.options <- list(algorithms = c(1,1,1), nbiter.saemix = c(200,100), nb.chains=1, save=FALSE, save.graphs=FALSE, seed = 1234, displayProgress = FALSE)

# RQ1 Use Gompertz, logistic, or exponential decay curve to model weight gain? ----

# _1.) Define data ----

ModelData.RQ1 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )

# _2.) Define model objects ----

# __a) Gompertz ----

GompertzModel.RQ1 <- saemixModel(model = gompertz.model, 
                                 description = 'Gompertz', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), # starting values for each of the four growth parameters
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE), # which of the four parameter shoul be estimated? All but the last here
                                 transform.par = c(0, 0, 0, 0) # distribution of each of the four parameter
                                 )

# __b) Logistic ----

LogisticModel.RQ1 <- saemixModel(model = logistic.model,
                                 description = 'Logistic growth', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), 
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE),
                                 transform.par=c(0,0,0,0)
                                 )


# __c) Exponential decay  ----

DecayModel.RQ1 <- saemixModel(model = decay.model,
                              description= "Exponential decay", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,0,0,0,1,0,0,0,1), ncol=3, byrow = TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")


# _3.) Fit model ----

# __a) Gompertz ----

GompertzFit.RQ1 <- saemix(GompertzModel.RQ1, ModelData.RQ1, saemix.options)

# Likelihood computed by linearisation
# -2LL= 944.1051 
# AIC = 966.1051 
# BIC = 987.1373 
# 
# Likelihood computed by importance sampling
# -2LL= 941.8332 
# AIC = 963.8332 
# BIC = 984.8654 

# __b) Logistic ----

LogisticFIT.RQ1 <- saemix(LogisticModel.RQ1, ModelData.RQ1, saemix.options)

# Likelihood computed by linearisation
# -2LL= 946.9261 
# AIC = 968.9261 
# BIC = 989.9584 
# 
# Likelihood computed by importance sampling
# -2LL= 946.2228 
# AIC = 968.2228 
# BIC = 989.2551 <- all values higher then above - so worse

# __c) Decay ----

DecayFIT.RQ1 <- saemix(DecayModel.RQ1, ModelData.RQ1, saemix.options)

# Likelihood computed by linearisation
# -2LL= 906.7097 
# AIC = 920.7097 
# BIC = 934.0938 
# 
# Likelihood computed by importance sampling
# -2LL= 906.5984 
# AIC = 920.5984 
# BIC = 933.9826 <- all values lower then both above - so best of the three

# _4.) Plot models ----

# Compute normalised prediction distribution errors (npde) and normalised prediction discrepancies (npd).

# __a) Gompertz ----

plot(GompertzFit.RQ1, plot.type="observations.vs.predictions" )
plot(GompertzFit.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.GompertzFit.RQ1 <- npdeSaemix(GompertzFit.RQ1) # skewness and kurtosis of normalised prediction discrepancies lower then in log model - but no normal distribution ?

# __b) Logistic ----

plot(LogisticFIT.RQ1, plot.type = "observations.vs.predictions")
plot(LogisticFIT.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.LogisticFIT.RQ1 <- npdeSaemix(LogisticFIT.RQ1) # normality of residual within range

# __b) Decay ----

plot(DecayFIT.RQ1, plot.type = "observations.vs.predictions")
plot(DecayFIT.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.DecayFIT.RQ1 <- npdeSaemix(DecayFIT.RQ1) # skew, curtosis of residuals best of all three - use this model



# RQ2: Does Sex have an association with the total growth, rate of approach to the upper asymptotic, or point of inflection (of the curve found optimal in answering the first research question?) ----

# _1.) Define data with sex covariate ----

ModelData.RQ2 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("AnimalSex")
  )

# _2.) Define model object ----

DecayModel.RQ2 <- saemixModel(model = decay.model,
                              description= "Exponential decay", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")




# _3.) Fit model ----

DecayFIT.RQ2 <- saemix(DecayModel.RQ2, ModelData.RQ2, saemix.options)

# _4.) Plot model ----

plot(DecayFIT.RQ2, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ2, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ2, plot.type="parameters.vs.covariates")
npde.DecayFIT.RQ2 <- npdeSaemix(DecayFIT.RQ2) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _5.) Compare models ----

DecayFIT.RQ1

# Fixed effects
# Parameter Estimate   SE     CV(%)
# A         25.3005  0.39486 1.56  
# B          1.2644  0.05365 4.24  
# k          0.0392  0.00154 3.93  
# a.1        0.7390  0.03702 5.01  
# 
# Variance of random effects
# Parameter Estimate   SE     CV(%)
# omega2.A  0.010457 0.00215  20.6 
# omega2.B  0.005064 0.00438  86.4 
# omega2.k  0.000583 0.00364 625.7 

DecayFIT.RQ2

# Fixed effects
# Parameter         Estimate   SE     CV(%) p-value
#  A                 23.1599  0.38612   1.67 -      
#  beta_AnimalSex(A)  0.1746  0.02097  12.01 0.000   # upper asymptote higher
#  B                  1.1081  0.08347   7.53 -      
#  beta_AnimalSex(B)  0.0962  0.09323  96.89 0.302   # starting value higher - insignificant
#  k                  0.0349  0.00267   7.65 -      
#  beta_AnimalSex(k)  0.0438  0.09197 209.79 0.634   # approach steeper  - insignificant
#  a.1                0.7285  0.03997   5.49 -      
#   
# 
# Correlation matrix of random effects
#         omega2.A omega2.B omega2.k
# omega2.A  1.000   -0.078   -0.176  # upper asymptote negatively correlated with starting value   
# omega2.B -0.078    1.000    0.952   
# omega2.k -0.176    0.952    1.000  

# second model is better in log-likelihood ratio test - see Likelihood Ratio Tests
teststatRQ12 <- -2 * (as.numeric(logLik(DecayFIT.RQ1)) - as.numeric(logLik(DecayFIT.RQ2)))
p.val <- pchisq(teststatRQ12, df = 3, lower.tail = FALSE)
p.val

# **Continue here** ----

# RQ3: What are the effects of diet within each sex? ----

# _1.) Define data with sex covariate ----

ModelData.RQ3 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("AnimalSex", "FatherDiet", "MotherDiet")
)

# _2.) Define model object ----

DecayModel.RQ3 <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")




# _3.) Fit model ----

DecayFIT.RQ3 <- saemix(DecayModel.RQ3, ModelData.RQ3, saemix.options)

# Fixed effects
# Parameter          Estimate   SE     CV(%) p-value
#   A                  24.0758  0.55485   2.30 -      
#   beta_AnimalSex(A)   0.1744  0.01899  10.89 0.00000
#   beta_FatherDiet(A) -0.0498  0.01907  38.27 0.00898 # regardless of sex, fathers HFD lowered body weight by 4%
#   beta_MotherDiet(A) -0.0302  0.01829  60.65 0.09916

# Likelihood computed by linearisation
# -2LL= 832.928 
# AIC= 870.928 
# BIC= 907.2564 
# Likelihood computed by importance sampling
# -2LL= 831.6653 
# AIC= 869.6653 
# BIC= 905.9937

# _4.) Plot model ----

plot(DecayFIT.RQ3, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ3, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ3, plot.type="parameters.vs.covariates", ask=TRUE)
npde.DecayFIT.RQ3 <- npdeSaemix(DecayFIT.RQ2) # skewness and kurtosis of normalised prediction discrepancies lower then in log model


# _5.) Compare models ----


# adding diets resulted in Likelihood increases


teststatRQ12 <- -2 * (as.numeric(logLik(DecayFIT.RQ2)) - as.numeric(logLik(DecayFIT.RQ3)))
p.val <- pchisq(teststatRQ12, df = 3, lower.tail = FALSE)
p.val # significant














# **Sex split below** ----

# RQ4: What are the effects within each sex? ----

# _1.) Define model (Logistic) ----

LogisticModel.RQ3 <- saemixModel(model = logistic.model,
                                 description = 'Logistic', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), 
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE),
                                 transform.par=c(0,0,0,0)
                                 )

# _2.) Define data (Logistic) ----

LogisticData.RQ3.m <- saemixData(
  name.data = (mice_f1_slct %>% dplyr::filter(AnimalSex == "m")), header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )

LogisticData.RQ3.f <- saemixData(
  name.data = (mice_f1_slct %>% dplyr::filter(AnimalSex == "f")), header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )


# _3.) Fit model ----

LogisticModel.RQ3.m <- saemix(LogisticModel.RQ1, LogisticData.RQ3.m, NLMEGM.options)

# Fixed effects
# Parameter Estimate   SE     CV(%)
# TtlGrwth  18.6298  3.04572 16.35  [total weight gain over time (pop ref NA g)] 
# Apprch     0.0714  0.00861 12.05  [steepness (pop ref )]
# Timing    32.4976  4.41948 13.60  [growth inflection (pop ref 27.6 d) ]
# LwrAsy     7.8448  2.90706 37.06  [not estimated  - birth weight ]
# a.1        1.0240  0.07689  7.51 

LogisticModel.RQ3.f <- saemix(LogisticModel.RQ1, LogisticData.RQ3.f, NLMEGM.options)

# Fixed effects
# Parameter Estimate   SE     CV(%)
# TtlGrwth  14.2765  2.19220 15.36  [less total gain then in males, less variable]
# Apprch     0.0747  0.00864 11.57  [less steeper approach, less variable]
# Timing    32.6176  4.06872 12.47  [growth inflection slightly later then in males]
# LwrAsy     7.7248  2.10042 27.19  [birth weight slightly lower then in males, less variable]
# a.1        0.6232  0.05537  8.88 


# _3.) Plot model ----

plot(LogisticModel.RQ3.m, plot.type="observations.vs.predictions" )
plot(LogisticModel.RQ3.f, plot.type="observations.vs.predictions" )

plot(LogisticModel.RQ3.m, plot.type = "both.fit",  ilist = 10:21, smooth = TRUE)
plot(LogisticModel.RQ3.f, plot.type = "both.fit",  ilist = 10:21, smooth = TRUE)

npde.LogisticModel.RQ3.m <- npdeSaemix(LogisticModel.RQ3.m) # skewness and kurtosis of normalised prediction discrepancies lower then in log model
npde.LogisticModel.RQ3.f <- npdeSaemix(LogisticModel.RQ3.f) # skewness and kurtosis of normalised prediction discrepancies lower then in log model


# RQ4: What are the effects of diet within each sex? ----

# _1.) Define data (Logistic) ----

LogisticData.RQ4.m <- saemixData(
  name.data = (mice_f1_slct %>% dplyr::filter(AnimalSex == "m")), header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("MotherDiet", "FatherDiet")
)

LogisticData.RQ4.f <- saemixData(
  name.data = (mice_f1_slct %>% dplyr::filter(AnimalSex == "f")), header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("MotherDiet", "FatherDiet")
)



# _2.) Define model (Logistic) ----

LogisticModel.RQ4 <- saemixModel(model = logistic.model,
                                 description = 'Logistic', 
                                 psi0 = c(TtlGrwth = 0, Apprch = 0, Timing = 0, LwrAsy = 0), 
                                 covariance.model = matrix( c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0), ncol = 4, byrow = TRUE),
                                 covariate.model = matrix(c(1, 0, 0, 0, 1, 0, 0, 0), ncol = 4, byrow = TRUE),
                                 transform.par=c(0,0,0,0)
                                 )

# _2.) Fit models (Logistic) ----

LogisticModel.RQ4.m <- saemix(LogisticModel.RQ4, LogisticData.RQ4.m, NLMEGM.options)
LogisticModel.RQ4.f <- saemix(LogisticModel.RQ4, LogisticData.RQ4.f, NLMEGM.options)

