#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis testing - using {saemix}"
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

library("here")    # handle path names
library("dplyr")   # handle data more easily

library("lattice") # create trellis graphs
library("ggpubr")  # save trellis graphs
library("ggplot2") # save trellis graphs

library("saemix")  # model non-linear mixed-effect dependencies
library("npde")    # model non-linear mixed-effect dependencies

# Get and check data ----

# _1.) Read in data ----

mice_f0_slct <- readRDS( file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS( file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Add litter size to f1 ----

mice_f1_slct <- left_join(mice_f1_slct, {mice_f0_slct %>% dplyr::select(AnimalId, MatingWith, LitterSize) %>% distinct}, by = c("MotherId" = "AnimalId", "FatherId" = "MatingWith"))

# _3.) Check data ----

mice_f1_slct %>% dplyr::select(MeasurementDay, BodyWeightG, AnimalId, AnimalSex) %>% arrange(AnimalSex, AnimalId, MeasurementDay) %>% print(n = Inf)

plot_data_check <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, groups = AnimalSex, data = mice_f1_slct, xlab = "day [d]", ylab = "body weight [g]", auto.key = list(title = "sex"))
plot_data_check

ggsave("015_r_use_saemix__data_check.pdf", plot = ggarrange(plot_data_check), path = here("../manuscript/display_items"),
       scale = 1, width = 9, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# RQ1 Which function is best suited to model weight gain (in a null model)? Gompertz, logistic, or exponential approach curve? ----

# _1.) Define possibly applicable model functions ----

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
  ypred  <- LwrAsy + TtlGrwth / (1+exp ( -Apprch * (t-Timing)))
  
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


# _2.) Set modelling options ----

saemix.options <- list(algorithms = c(1,1,1), nbiter.saemix = c(200,100), nb.chains=1, save=FALSE, save.graphs=FALSE, seed = 1234, displayProgress = FALSE)

# _3.) Define data (already grouped)----

ModelData.RQ1 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )

# _4.) Define model objects ----

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


# __c) Exponential approach  ----

DecayModel.RQ1 <- saemixModel(model = decay.model,
                              description= "Exponential decay", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,0,0,0,1,0,0,0,1), ncol=3, byrow = TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _5.) Fit models ----

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

# __c) Approach ----

DecayFIT.RQ1 <- saemix(DecayModel.RQ1, ModelData.RQ1, saemix.options)

# Likelihood computed by linearisation
# -2LL= 906.7097 
# AIC = 920.7097 
# BIC = 934.0938 
# 
# Likelihood computed by importance sampling
# -2LL= 906.5984 
# AIC = 920.5984 
# BIC = 933.9826  <- all values lower then both above - so best of the three

# _6.) Plot models ----

# Compute normalized prediction distribution errors (npde) and normalized prediction discrepancies (npd).

# __a) Gompertz ----

plot(GompertzFit.RQ1, plot.type="observations.vs.predictions" )
plot(GompertzFit.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.GompertzFit.RQ1 <- npdeSaemix(GompertzFit.RQ1) # skewness and kurtosis of normalised prediction discrepancies lower then in log model - but no normal distribution ?

# __b) Logistic ----

plot(LogisticFIT.RQ1, plot.type = "observations.vs.predictions")
plot(LogisticFIT.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.LogisticFIT.RQ1 <- npdeSaemix(LogisticFIT.RQ1) # good - normality of residual within range

# __b) Approach model  ----

plot(DecayFIT.RQ1, plot.type = "observations.vs.predictions")
plot(DecayFIT.RQ1, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
npde.DecayFIT.RQ1 <- npdeSaemix(DecayFIT.RQ1) # best - normality of residual within range

# _7.) Answer RQ1 ----

# The approach model is best for modeling body weights, and fits well to the data.

# RQ2: Does Sex have an association with the total weight gain? ----

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
npde.DecayFIT.RQ2 <- npdeSaemix(DecayFIT.RQ2) # residuals still considered normally distributed albeit weakly

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
#  beta_AnimalSex(B)  0.0962  0.09323  96.89 0.302   # insignificant - starting value higher 
#  k                  0.0349  0.00267   7.65 -      
#  beta_AnimalSex(k)  0.0438  0.09197 209.79 0.634   # insignificant- approach steeper  
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

# _6.) Answer RQ2 ----

# Sex has an influence on the body weight, the upper assymptote is higher for males.

# RQ3: Do Sex and litter size have an association with the total weight gain? ----

# _1.) Define data with sex and litter size covariates ----

ModelData.RQ3 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("AnimalSex", "LitterSize")
  )

# _2.) Define model object ----

DecayModel.RQ3 <- saemixModel(model = decay.model,
                              description= "exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1, 1,1 ), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _3.) Fit model ----

DecayFIT.RQ3 <- saemix(DecayModel.RQ3, ModelData.RQ3, saemix.options)

# _4.) Plot model ----

plot(DecayFIT.RQ3, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ3, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ3, plot.type="parameters.vs.covariates", ask = TRUE)
npde.DecayFIT.RQ3 <- npdeSaemix(DecayFIT.RQ3) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _5.) Compare models ----

DecayFIT.RQ2
DecayFIT.RQ3 # <- BIC AIC slightly better then in sex alone

teststatRQ12 <- -2 * (as.numeric(logLik(DecayFIT.RQ2)) - as.numeric(logLik(DecayFIT.RQ3)))
p.val <- pchisq(teststatRQ12, df = 3, lower.tail = FALSE)
p.val #  (significant when adding litter size to sex)

# _6.) Answer RQ3 ----

# _a) In text ----

# Only sex, not litter size influence the body weight significantly. 
# (Litter is borderline insignificant, perhaps influences k parameter)

# _b) In plot ----

coefficients(DecayFIT.RQ3)$fixed

# male mice growth
a = coefficients(DecayFIT.RQ3)$fixed[1] * 1.19165
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black")

# female mice growth
a = coefficients(DecayFIT.RQ3)$fixed[1]
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
  xlab = "days [d]", ylab = "weight [g]", col = "red", add = TRUE)

# male mice growth - increased litter size
a = coefficients(DecayFIT.RQ3)$fixed[1] * 1.19165
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3] * (1-0.09171)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", add = TRUE, lty = "dashed")

# female mice growth  - increased litter size
a = coefficients(DecayFIT.RQ3)$fixed[1]
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3] * (1-0.09171)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "red", add = TRUE, lty = "dashed")

legend(65, 20.5, legend=c("male", "male - larger litter", "female", "female - larger litter"),
       col=c("black", "black", "red", "red"), lty = c(1,2,1,2), cex=0.8)


# RQ4: What are the effects of diet within each sex and, and dependent on litter sizes? ----

# _1.) Define data with sex covariate ----

ModelData.RQ4 <- saemixData(
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("AnimalSex", "LitterSize", "FatherDiet", "MotherDiet")
  )

# _2.) Define model object ----

DecayModel.RQ4 <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _3.) Fit model ----

DecayFIT.RQ4 <- saemix(DecayModel.RQ4, ModelData.RQ4, saemix.options)

# _4.) Inspect model ----

DecayFIT.RQ4
# Fixed effects CV(%) is mostly bad

# Parameter          Estimate   SE    CV(%)  p-value 
# A                  21.55263 0.9979    4.63 -       
# beta_AnimalSex(A)   0.17613 0.0161    9.13 0.00e+00 *
# beta_LitterSize(A)  0.02362 0.0099   41.89 1.70e-02 * 0.0170
# beta_FatherDiet(A) -0.05120 0.0191   37.33 7.39e-03 * 0.00739
# beta_MotherDiet(A) -0.00151 0.0180 1197.06 9.33e-01 
# B                   3.24019 0.8274   25.54 -       
# beta_AnimalSex(B)   0.10238 0.0692   67.62 1.39e-01
# beta_LitterSize(B) -0.22072 0.0509   23.04 1.42e-05 * 0.0000142
# beta_FatherDiet(B)  0.13391 0.0790   59.03 9.03e-02 
# beta_MotherDiet(B) -0.29198 0.0794   27.20 2.36e-04 * 0.000236
# k                   0.09233 0.0200   21.69 -       
# beta_AnimalSex(k)   0.04126 0.0718  173.94 5.65e-01 
# beta_LitterSize(k) -0.21285 0.0448   21.07 2.07e-06 *
# beta_FatherDiet(k)  0.14874 0.0893   60.05 9.59e-02
# beta_MotherDiet(k) -0.20287 0.0852   42.00 1.73e-02 * 
# a.1                 0.77995 0.0402    5.15 -     

# _5.) Plot model ----

plot(DecayFIT.RQ4, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ4, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ4, plot.type="parameters.vs.covariates", ask=TRUE)
npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ4) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _6.) Compare models ----

DecayFIT.RQ3
DecayFIT.RQ4 # <- BIC AIC with diets worse then sex and litter size alone

teststatRQ12 <- -2 * (as.numeric(logLik(DecayFIT.RQ3)) - as.numeric(logLik(DecayFIT.RQ4)))
p.val <- pchisq(teststatRQ12, df = 3, lower.tail = FALSE)
p.val # model with diets is not distinctly different from when not adding diets

# _7.) Answer RQ4 ----

# Adding diet to the model with sex and litter size does not improve, but worsens it.
# Adding diet to the model with sex and litter does not make it significantly different
# from the respective null-model. The model is likely overfitted.

# RQ5: What are the effects of diet within each sex individually dependent on litter sizes? ----

# _1.) Split data by sexes ----

# To avoid overparametrisations.

ModelData.RQ5.male <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "m")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("LitterSize", "FatherDiet", "MotherDiet")
  )

ModelData.RQ5.female <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "f")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("LitterSize", "FatherDiet", "MotherDiet")
  )

# _2.) Define model object ----

DecayModel.RQ5 <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _3.) Fit models ----

DecayFIT.RQ5.male <- saemix(DecayModel.RQ5, ModelData.RQ5.male, saemix.options)

# Fixed effects
# Parameter          Estimate    SE    CV(%)  p-value
# A                  28.566939 2.5112    8.79 -      
# beta_LitterSize(A)  0.000838 0.0184 2199.03 0.9637 
# beta_FatherDiet(A) -0.035472 0.0285   80.24 0.2127 
# beta_MotherDiet(A) -0.053695 0.0285   53.04 0.0594 *
# B                   1.906307 0.7931   41.60 -      
# beta_LitterSize(B) -0.090955 0.0844   92.83 0.2814 
# beta_FatherDiet(B)  0.012596 0.1280 1016.57 0.9216 
# beta_MotherDiet(B)  0.059983 0.1354  225.69 0.6577 
# k                   0.045201 0.0160   35.44 -      
# beta_LitterSize(k) -0.059031 0.0739  125.14 0.4242 
# beta_FatherDiet(k)  0.060994 0.1139  186.66 0.5921 
# beta_MotherDiet(k)  0.173567 0.1191   68.64 0.1451 
# a.1                 0.834962 0.0589    7.05 - 

DecayFIT.RQ5.female <- saemix(DecayModel.RQ5, ModelData.RQ5.female, saemix.options)

# Fixed effects
# Parameter          Estimate   SE    CV(%) p-value
# A                  22.41907 1.2216   5.45 -      
# beta_LitterSize(A)  0.00652 0.0116 178.40 0.57510
# beta_FatherDiet(A) -0.02591 0.0204  78.61 0.20333
# beta_MotherDiet(A) -0.00363 0.0194 535.06 0.85174
# B                   2.87427 0.9309  32.39 -      
# beta_LitterSize(B) -0.13740 0.0688  50.09 0.04590
# beta_FatherDiet(B) -0.03240 0.1159 357.75 0.77984
# beta_MotherDiet(B) -0.29874 0.1095  36.65 0.00637 * but large variability across subjects
# k                   0.07869 0.0207  26.27 -      
# beta_LitterSize(k) -0.11206 0.0558  49.81 0.04469
# beta_FatherDiet(k) -0.04751 0.1003 211.08 0.63567
# beta_MotherDiet(k) -0.21824 0.0945  43.28 0.02086 * but large variability across subjects
# a.1                 0.52514 0.0433   8.25 -      

# _4.) Plot model fits ----

# all looking good
plot(DecayFIT.RQ5.female, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ5.male, plot.type="observations.vs.predictions" )

plot(DecayFIT.RQ5.female, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ5.male, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)

plot(DecayFIT.RQ5.female, plot.type="parameters.vs.covariates", ask=TRUE)
plot(DecayFIT.RQ5.male, plot.type="parameters.vs.covariates", ask=TRUE)

npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ5.female) # skewness and kurtosis of normalised prediction discrepancies lower then in log model
npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ5.male) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _5.) Plot model predictions ----

# **continue here** ----

# testing male model

# male mice growth
a = 28.566939
b = 1.906307
k = 0.045201

curve(
  a * (1 - b * exp(-k * x)),
  from = min(mice_f1_slct$MeasurementDay),
  to = max(mice_f1_slct$MeasurementDay),
  xlab = "days [d]",
  ylab = "weight [g]",
  main = "male mice weight gain (C57BL6/NTac)"
)

# male mice growth - mother on HFD
a <- a
b <- b 
k <- k + (0.173567 * k)

curve(
  a * (1 - b * exp(-k * x)),
  add = TRUE,
  col = "red"
)




# female mice growth

a = 22.41907
b = 2.87427
k = 0.07869

curve(
  a * (1 - b * exp(-k * x)),
  from = min(mice_f1_slct$MeasurementDay),
  to = max(mice_f1_slct$MeasurementDay),
  xlab = "days [d]",
  ylab = "weight [g]",
  main = "female mice weight gain (C57BL6/NTac)"
)

# female mice growth - with mother on HFD

a <- a
b <- b + (- 0.29874 * b)
k <- k + (- 0.21824 * k)
curve(
  a * (1 - b * exp(-k * x)),
  add = TRUE,
  col = "red"
)

# female mice growth - with increasing litter size + 1

a <- a
b <- b + (- 0.13740 * b)
k <- k
curve(
  a * (1 - b * exp(-k * x)),
  add = TRUE,
  col = "blue"
)



# _6.) Compare models ----

# no null model for comparison built yet, needed would be above data sets without diets

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "015_r_use_saemix.R.Rdata"))
renv::snapshot()

