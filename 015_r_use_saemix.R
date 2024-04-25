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

library("GGally")  # pair plot

library("saemix")  # model non-linear mixed-effect dependencies
library("npde")    # model non-linear mixed-effect dependencies

# _3.) Functions ----

get_p_from_seamix_lrt = function(sfit1, sfit2){
  
  # compute LRT test of two models - see 10.5964/meth.7061
  ts <- -2 * (as.numeric(logLik(sfit1)) - as.numeric(logLik(sfit2)))
  pv <-  pchisq(ts, df = 3, lower.tail = FALSE)
  
  return(pv)
  
}

# Get and check data ----

# _1.) Read in data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Add and correct litter sizes  ----

# __a) Add litter sizes ----

mice_f1_slct <- left_join(mice_f1_slct, {mice_f0_slct %>% dplyr::select(AnimalId, MatingWith, LitterSize) %>% distinct}, by = c("MotherId" = "AnimalId", "FatherId" = "MatingWith"))

# __b) Indicate descriptive status of variable ----

mice_f1_slct %<>% rename(LitterSizeDescription =  LitterSize)

# __c) Isolate pup counts for each sex ----

mice_f1_slct %<>% mutate(LitterSizeMale = as.double(sub("\\..*", "", LitterSizeDescription)))
mice_f1_slct %<>% mutate(LitterSizeFemale = as.double(sub(".*\\.", "", LitterSizeDescription)))

# __d) Redefine litter size  ----

# to match with previous model formulae

mice_f1_slct %<>% mutate(LitterSize = LitterSizeMale +  LitterSizeFemale)

# _3.) Check data ----

# __a) Table ----

mice_f1_slct %>% dplyr::select(MeasurementDay, BodyWeightG, AnimalSex, AnimalId, LitterSize) %>% arrange(AnimalSex, AnimalId, MeasurementDay) %>% print(n = Inf)

# __b) Drop superfluous factor levels ----

mice_f1_slct$MotherDiet <- droplevels(mice_f1_slct$MotherDiet)

# __c) XY plots ----

plot_data_check <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, groups = AnimalSex, data = mice_f1_slct, xlab = "day [d]", ylab = "body weight [g]", auto.key = list(title = "sex"))
plot_data_check

ggsave("015_r_use_saemix__data_check.pdf", plot = ggarrange(plot_data_check), path = here("../manuscript/display_items"),
       scale = 1, width = 9, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# __d) Check correlation ----

splom({mice_f1_slct %>% dplyr::select(LitterSize, FatherDiet,  MotherDiet, AnimalSex,  BodyWeightG)}, main = "F1 mice")

ggpairs({mice_f1_slct %>% dplyr::select(LitterSize, FatherDiet,  MotherDiet, AnimalSex,  BodyWeightG)}, 
        title = "F1 mice", 
        axisLabels = "show") 

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
  name.data = mice_f1_slct, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("AnimalSex", "LitterSize", "FatherDiet", "MotherDiet")
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

# __b) Logistic ----

LogisticFIT.RQ1 <- saemix(LogisticModel.RQ1, ModelData.RQ1, saemix.options)

# __c) Approach ----

DecayFIT.RQ1 <- saemix(DecayModel.RQ1, ModelData.RQ1, saemix.options)

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

# _7.) Compare models ----

compare.saemix(GompertzFit.RQ1, LogisticFIT.RQ1, DecayFIT.RQ1)  # DecayFIT.RQ1 is best
get_p_from_seamix_lrt(LogisticFIT.RQ1, DecayFIT.RQ1)            # and significantly so from previous one

# _8.) Answer RQ1 ----

# The approach model is best for modeling body weights, and fits well to the data.

# RQ2: Across all animals does sex have an association with the total weight gain? ----

# _1.) Define model object ----

DecayModel.RQ2 <- saemixModel(model = decay.model,
                              description= "Exponential decay", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _2.) Fit model ----

DecayFIT.RQ2 <- saemix(DecayModel.RQ2, ModelData.RQ1, saemix.options)

# _3.) Plot model ----

plot(DecayFIT.RQ2, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ2, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ2, plot.type="parameters.vs.covariates")
npde.DecayFIT.RQ2 <- npdeSaemix(DecayFIT.RQ2) # residuals still considered normally distributed albeit weakly

# _4.) Compare models ----

compare.saemix(DecayFIT.RQ1, DecayFIT.RQ2)        # DecayFIT.RQ2 is better ...
get_p_from_seamix_lrt(DecayFIT.RQ1, DecayFIT.RQ2) # ... and significantly so

# _5.) Check coefficients of improved model  ----

summary(DecayFIT.RQ2)

# _6.) Answer RQ2 ----

# Across all animals does sex has an influence on the body weight, the upper asymptote is higher for males by 17.5%.

# RQ3: Across all animals do sex and litter size have an association with the total weight gain? ----

# _1.) Define model object ----

DecayModel.RQ3 <- saemixModel(model = decay.model,
                              description= "exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1, 1,1 ), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _2.) Fit model ----

DecayFIT.RQ3 <- saemix(DecayModel.RQ3, ModelData.RQ1, saemix.options)

# _3.) Plot model ----

plot(DecayFIT.RQ3, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ3, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ3, plot.type="parameters.vs.covariates", ask = TRUE)
npde.DecayFIT.RQ3 <- npdeSaemix(DecayFIT.RQ3) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _4.) Compare models ----

compare.saemix(DecayFIT.RQ1, DecayFIT.RQ2, DecayFIT.RQ3)  # Last model is best....
get_p_from_seamix_lrt(DecayFIT.RQ2, DecayFIT.RQ3)         # ...and significantly so

# _5.) Check estimates ----

summary(DecayFIT.RQ3)

# _6.) Answer RQ3 ----

# Across all animals does adding litter to sex improves the model markedly, litter sizes reduces a by 2%. 

# RQ4: Across all animals what are the effects of parental diets sex and litter sizes? ----

# _1.) Define model object ----

DecayModel.RQ4 <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _2.) Fit model ----

DecayFIT.RQ4 <- saemix(DecayModel.RQ4, ModelData.RQ1, saemix.options)

# _3.) Plot model ----

plot(DecayFIT.RQ4, plot.type = "observations.vs.predictions" )
plot(DecayFIT.RQ4, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ4, plot.type = "parameters.vs.covariates", ask = TRUE) # <- use this in paper
npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ4) # skewness and kurtosis of normalized prediction discrepancies lower then in log model

# _4.) Compare models ----

compare.saemix(DecayFIT.RQ1, DecayFIT.RQ2, DecayFIT.RQ3, DecayFIT.RQ4) # last one is best, the most complex one
               
get_p_from_seamix_lrt(DecayFIT.RQ2, DecayFIT.RQ3) # ... it improved significantly - 0.0055129
get_p_from_seamix_lrt(DecayFIT.RQ3, DecayFIT.RQ4) # ... it improved significantly - 0.0003406611

# _5.) Show coefficients ----

summary(DecayFIT.RQ4) #  AIC = 855.1701, AIC = 855.6157 

#   Parameter Estimate     SE  CV(%) p-value
# 1                   A   26.746 1.4247   5.33       -
# 2   beta_AnimalSex(A)    0.181 0.0159   8.80   0.000
# 3  beta_LitterSize(A)   -0.017 0.0078  46.89   0.033 *
# 4  beta_FatherDiet(A)   -0.031 0.0204  65.54   0.127
# 5  beta_MotherDiet(A)   -0.035 0.0161  45.92   0.029 * 
# 6                   B    1.944 0.6240  32.09       -
# 7   beta_AnimalSex(B)    0.013 0.0951 749.43   0.894
# 8  beta_LitterSize(B)   -0.073 0.0431  59.22   0.091
# 9  beta_FatherDiet(B)    0.166 0.1005  60.55   0.099
# 10 beta_MotherDiet(B)    0.011 0.0909 855.42   0.907
# 11                  k    0.058 0.0154  26.51       -
# 12  beta_AnimalSex(k)   -0.043 0.0806 186.14   0.591
# 13 beta_LitterSize(k)   -0.074 0.0375  50.86   0.049 *
# 14 beta_FatherDiet(k)    0.219 0.0949  43.32   0.021 * - but check coeffcient plot (red dashed lines) - wonky estimate
# 15 beta_MotherDiet(k)    0.082 0.0788  96.35   0.299
# 16                a.1    0.727 0.0394   5.42       -
  
# _7.) Answer RQ4 ----

# Most complex model is best, including sex, litter, and both diets model, mother diet depresses a by 3.5%.

# _8) Plot complex and best model ----

# __a) Show values for plotting ----

coefficients(DecayFIT.RQ4)$fixed
summary(DecayFIT.RQ4)

# reference class for covariate AnimalSex :  f 
# reference class for covariate FatherDiet :  CD 
# reference class for covariate MotherDiet :  CD

# __b) Male mice growth ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 + 0.181)
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black")

# __c) Female mice growth ----

a = coefficients(DecayFIT.RQ3)$fixed[1]
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", add = TRUE)


# __d) Male mice growth - increased litter size by one----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 + 0.181) * (1 - 0.017)
b = coefficients(DecayFIT.RQ4)$fixed[2] 
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 - 0.074)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "blue", add = TRUE)

# __e) Female mice growth - increased litter size by one ----

a = coefficients(DecayFIT.RQ3)$fixed[1] * (1 + -0.017)
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3] * (1 - 0.074)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "blue", add = TRUE)

# **continue here with plotting ** ----

# __f) Male mice growth - mother diet ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 + 0.181)
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "darkorange")

# __g) Female mice growth  - mother diet ----

a = coefficients(DecayFIT.RQ3)$fixed[1]
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "darkorange", add = TRUE)


# __h) Male mice growth - father diet ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 + 0.181)
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "red")

# __i) Female mice growth  - mother diet ----

a = coefficients(DecayFIT.RQ3)$fixed[1]
b = coefficients(DecayFIT.RQ3)$fixed[2]
k = coefficients(DecayFIT.RQ3)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "red", add = TRUE)

# adjust and export  / "all effects are additive"
legend(65, 20.5, legend=c("male", "male - larger litter", "female", "female - larger litter"),
       col=c("black", "black", "red", "red"), lty = c(1,2,1,2), cex=0.8)


# **continue here with testing interactions ** ----



# RQ5: Can the complex model be improved by backward selection? ----

# _1.) Define model object ----

DecayModel.RQ5 <- backward.procedure(DecayFIT.RQ4, trace = TRUE)

# _2.) Fit model ----

DecayFIT.RQ5 <- saemix(DecayModel.RQ5, ModelData.RQ1, saemix.options)

# _3.) Plot model ----

plot(DecayFIT.RQ5, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ5, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE) # poor fit
plot(DecayFIT.RQ5, plot.type="parameters.vs.covariates", ask = TRUE)
DecayFIT.RQ5 <- npdeSaemix(DecayFIT.RQ5) # poor fit

# _4.) Compare models ----

compare.saemix(DecayFIT.RQ1, DecayFIT.RQ2, DecayFIT.RQ3, DecayFIT.RQ4, DecayFIT.RQ5)  # Last model is worst....
get_p_from_seamix_lrt(DecayFIT.RQ5, DecayFIT.RQ4)         # ...and significantly so

# _5.) Check estimates ----

summary(DecayFIT.RQ5)

# _6.) Answer RQ5 ----

# Backward selection did not find a better-fitting model
























# *** below - check for factor intercations ****

# RQ5 - relevant for DEG are only males -  Part 1: What are the effects of diet within each sex individually dependent on litter sizes? ----

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
summary(DecayFIT.RQ5.male)

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
summary(DecayFIT.RQ5.female)

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

npde.DecayFIT.RQ5 <- npdeSaemix(DecayFIT.RQ5.female) # skewness and kurtosis of normalised prediction discrepancies lower then in log model
npde.DecayFIT.RQ5 <- npdeSaemix(DecayFIT.RQ5.male) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _5.) Plot model predictions ----

# __a) females ----

# female mice, both parents on low caloric diet 
a = coefficients(DecayFIT.RQ5.female)$fixed[1]
b = coefficients(DecayFIT.RQ5.female)$fixed[2]
k = coefficients(DecayFIT.RQ5.female)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", lty = "solid")

# female mice with increasing litter size, both parents on low caloric diet 
a = coefficients(DecayFIT.RQ5.female)$fixed[1]
b = coefficients(DecayFIT.RQ5.female)$fixed[2] * (1-0.13740) 
k = coefficients(DecayFIT.RQ5.female)$fixed[3] * (1-0.11206)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "red", lty = "solid", add = TRUE)

# female mice, mother on HCD 
a = coefficients(DecayFIT.RQ5.female)$fixed[1]
b = coefficients(DecayFIT.RQ5.female)$fixed[2] * (1-0.29874)
k = coefficients(DecayFIT.RQ5.female)$fixed[3] * (1-0.21824)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", lty = "dashed", add = TRUE)

# female mice with increasing litter size, mother on HCD
a = coefficients(DecayFIT.RQ5.female)$fixed[1]
b = coefficients(DecayFIT.RQ5.female)$fixed[2] * (1-0.13740) * (1-0.29874)
k = coefficients(DecayFIT.RQ5.female)$fixed[3] * (1-0.11206) * (1-0.21824)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "red", lty = "dashed", add = TRUE)

legend(42, 18, legend=c("female", "female, larger litter", "female, mother on HCD", "female, larger litter, mother on HCD"),
       col=c("black", "red", "black", "red"), lty = c(1,1,2,2), cex=0.8)

# __b) males ----

# male mice, both parents on low caloric diet 
a = coefficients(DecayFIT.RQ5.male)$fixed[1]
b = coefficients(DecayFIT.RQ5.male)$fixed[2]
k = coefficients(DecayFIT.RQ5.male)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", lty = "solid")

# male mice, mother on HCD 
a = coefficients(DecayFIT.RQ5.male)$fixed[1] * (1-0.053695)
b = coefficients(DecayFIT.RQ5.male)$fixed[2]
k = coefficients(DecayFIT.RQ5.male)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", lty = "dashed", add = TRUE)

legend(42, 18, legend=c("male", "male, mother on HCD"),
       col=c("black", "black"), lty = c(1,2), cex=0.8)

# _6.) Compare models ----

# see below in RQ8

# RQ6 - relevant for DEG are only males - Part 2: Disentangle Litter Size and Diet effect on Body weight ----

# Same model as in RQ5 but without diets, but only with litter size

ModelData.RQ6.male <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "m")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("LitterSize")
)

ModelData.RQ6.female <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "f")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("LitterSize")
)

# _2.) Define model object ----

DecayModel.RQ6 <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _3.) Fit models ----

DecayFIT.RQ6.male <- saemix(DecayModel.RQ6, ModelData.RQ6.male, saemix.options)
summary(DecayFIT.RQ6.male)

# Litter size not significant

DecayFIT.RQ6.female <- saemix(DecayModel.RQ6, ModelData.RQ6.female, saemix.options)
summary(DecayFIT.RQ6.female)

# Litter size significant for b and k

# _4.) Plot model fits ----

# all looking good
plot(DecayFIT.RQ6.female, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ6.male, plot.type="observations.vs.predictions" )

plot(DecayFIT.RQ6.female, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ6.male, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)

plot(DecayFIT.RQ6.female, plot.type="parameters.vs.covariates", ask=TRUE)
plot(DecayFIT.RQ6.male, plot.type="parameters.vs.covariates", ask=TRUE)

npde.DecayFIT.RQ6 <- npdeSaemix(DecayFIT.RQ6.female) # residuals not normal - skewed
npde.DecayFIT.RQ6 <- npdeSaemix(DecayFIT.RQ6.male)   # residuals less sked then among females but still skewed

# _5.) Plot model predictions ----

# __a) females ----

# female mice, both parents on low caloric diet 
a = coefficients(DecayFIT.RQ6.female)$fixed[1]
b = coefficients(DecayFIT.RQ6.female)$fixed[2]
k = coefficients(DecayFIT.RQ6.female)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "black", lty = "solid")

# female mice with increasing litter size, both parents on low caloric diet 
a = coefficients(DecayFIT.RQ6.female)$fixed[1]
b = coefficients(DecayFIT.RQ6.female)$fixed[2] * (1-0.14491) 
k = coefficients(DecayFIT.RQ6.female)$fixed[3] * (1-0.13252)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      xlab = "days [d]", ylab = "weight [g]", col = "red", lty = "solid", add = TRUE)

# _6.) Compare models ----

DecayFIT.RQ5.female
DecayFIT.RQ6.female # <- simpler model is better 

teststatRQ56.female <- -2 * (as.numeric(logLik(DecayFIT.RQ5.female)) - as.numeric(logLik(DecayFIT.RQ6.female)))
p.val <- pchisq(teststatRQ56.female, df = 3, lower.tail = FALSE)
p.val # model with diets is not distinctly different from when not adding diets

DecayFIT.RQ5.male
DecayFIT.RQ6.male # <- simpler model is better 

teststatRQ56.male <- -2 * (as.numeric(logLik(DecayFIT.RQ5.male)) - as.numeric(logLik(DecayFIT.RQ6.male)))
p.val <- pchisq(teststatRQ56.male, df = 3, lower.tail = FALSE)
p.val # model with diets is not distinctly different from when not adding diets

# RQ7 - relevant for DEG are only males -  Part 3: - Disentangle Litter Size and Diet effect on Body weight ----

# Same model as in RQ5 but only diets, not with litter size

ModelData.RQ7.male <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "m")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("FatherDiet", "MotherDiet")
)

ModelData.RQ7.female <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "f")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("FatherDiet", "MotherDiet")
)

# _2.) Define model object ----

DecayModel.RQ7 <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _3.) Fit models ----

DecayFIT.RQ7.male <- saemix(DecayModel.RQ7, ModelData.RQ7.male, saemix.options)
summary(DecayFIT.RQ7.male)

# Mother's diet significant on all parameters, father's on some 

DecayFIT.RQ7.female <- saemix(DecayModel.RQ7, ModelData.RQ7.female, saemix.options)
summary(DecayFIT.RQ7.female)

# Litter size significant for b and k

# Mother's diet significant b only

# _4.) Plot model fits ----

# all looking good
plot(DecayFIT.RQ7.female, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ7.male, plot.type="observations.vs.predictions" )

plot(DecayFIT.RQ7.female, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ7.male, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)

plot(DecayFIT.RQ7.female, plot.type="parameters.vs.covariates", ask=TRUE)
plot(DecayFIT.RQ7.male, plot.type="parameters.vs.covariates", ask=TRUE)

npde.DecayFIT.RQ7 <- npdeSaemix(DecayFIT.RQ7.female) # residuals not normal - skewed
npde.DecayFIT.RQ7 <- npdeSaemix(DecayFIT.RQ7.male)   # residuals less sked then among females but still skewed

# _5.) Compare models ----

# __a) IC's ----

DecayFIT.RQ5.female # litter size and diets - same as below
DecayFIT.RQ6.female # litter size - same as above
DecayFIT.RQ7.female # diet only - worst

DecayFIT.RQ5.male # litter size and diets - worst 
DecayFIT.RQ6.male # litter size only - best
DecayFIT.RQ7.male # diet only - second best

# __b) LRTs ----

# females - litter size only vs diet only
teststatRQ67.female <- -2 * (as.numeric(logLik(DecayFIT.RQ6.female)) - as.numeric(logLik(DecayFIT.RQ7.female)))
p.val <- pchisq(teststatRQ67.female, df = 3, lower.tail = FALSE)
p.val # model with diets is not distinctly different from model with litter sizes for females

# males - litter size only vs diet only
teststatRQ67.male <- -2 * (as.numeric(logLik(DecayFIT.RQ6.male)) - as.numeric(logLik(DecayFIT.RQ7.male)))
p.val <- pchisq(teststatRQ67.male, df = 3, lower.tail = FALSE)
p.val # model with diets is not distinctly different from model with litter sizes for females

# females - diet vs diet and litter size
get_p_from_seamix_lrt(DecayFIT.RQ7.female, DecayFIT.RQ5.female) # sign differnce when adding litter
                                                                # model better with litter size
# males - diet vs diet and litter size
teststatRQ57.male <- -2 * (as.numeric(logLik(DecayFIT.RQ5.male)) - as.numeric(logLik(DecayFIT.RQ7.male)))
p.val <- pchisq(teststatRQ57.male, df = 3, lower.tail = FALSE)
p.val # no differnces


# RQ8 - relevant for DEG are only males -  Part 4: Get matching null models for RQ5 thru RQ7  ----

# Same model as in RQ5 thru RQ7 but no covariates

ModelData.RQ8.male.null <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "m")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )

ModelData.RQ8.female.null <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "f")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay"
  )

# _2.) Define model object ----

DecayModel.RQ8.null <- saemixModel(model = decay.model,
                              description= "Exponential approach", 
                              psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _3.) Fit models ----

DecayFIT.RQ8.male.null <- saemix(DecayModel.RQ8.null, ModelData.RQ8.male.null, saemix.options)
summary(DecayFIT.RQ8.male.null) # AIC = 562.4538, AIC = 561.7869

DecayFIT.RQ8.female.null <- saemix(DecayModel.RQ8.null, ModelData.RQ8.female.null, saemix.options)
summary(DecayFIT.RQ8.female.null) # AIC = 287.3807, AIC = 287.5092

# _4.) Compare models ----

# __a) IC's ----

DecayFIT.RQ5.female # litter size and diets - same as below - AIC 286.7304, 287.3118 - better then null
DecayFIT.RQ6.female # litter size - same as above - AIC 286.5673, 284.9584 - better then null
DecayFIT.RQ7.female # diet only - worst - AIC - 289.3388, 289.3458 - worse then null

DecayFIT.RQ5.male # litter size and diets - worst - AIC= 568.6449, AIC= 569.5282 - worse then null
DecayFIT.RQ6.male # litter size only - best - AIC= 563.8716, AIC= 563.9557 - worse then null
DecayFIT.RQ7.male # diet only - second best - AIC= 565.2276, AIC= 564.3591 - worse then null

# __b) LRTs ----

get_p_from_seamix_lrt(DecayFIT.RQ8.male.null, DecayFIT.RQ6.male) # no sign difference - with litter size
get_p_from_seamix_lrt(DecayFIT.RQ8.male.null, DecayFIT.RQ7.male) # sign difference - with diet - gets worse
get_p_from_seamix_lrt(DecayFIT.RQ8.male.null, DecayFIT.RQ5.male) # sign difference - with both parameters - gets worse
get_p_from_seamix_lrt(DecayFIT.RQ7.male, DecayFIT.RQ5.male)      # no sign differences - with or without litter with diets

get_p_from_seamix_lrt(DecayFIT.RQ8.female.null, DecayFIT.RQ6.female) # sign difference - when only litter - gets better
get_p_from_seamix_lrt(DecayFIT.RQ8.female.null, DecayFIT.RQ7.female) # sign difference - when only diet - gets worse
get_p_from_seamix_lrt(DecayFIT.RQ8.female.null, DecayFIT.RQ5.female) # sign difference - when both litter and diet - gets better
get_p_from_seamix_lrt(DecayFIT.RQ6.female, DecayFIT.RQ5.female)      # sign differences - litter only is better

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "015_r_use_saemix.Rdata"))
renv::snapshot()

