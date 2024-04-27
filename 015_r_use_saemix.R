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

mice_f1_slct %<>% mutate(LitterSize = LitterSizeMale + LitterSizeFemale)

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
                              description= "exponential approach", 
                              psi0 = matrix( c(28,2,0.05, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
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
                              description= "exponential approach", 
                              psi0 = matrix( c(28,2,0.05, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
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
                              psi0 = matrix( c(28,2,0.05, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1, 1,1), ncol=3, byrow=TRUE),
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
                              description= "exponential approach", 
                              psi0 = matrix( c(28,2,0.05, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                              transform.par = c(1,1,1), 
                              fixed.estim = c(1,1,1),
                              covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                              covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                              omega.init = matrix(c(1,0,0,0,1,0,0,0,1),ncol=3, byrow=TRUE), 
                              error.model="constant")

# _2.) Fit model ----

DecayFIT.RQ4 <- saemix(DecayModel.RQ4, ModelData.RQ1, saemix.options)

# _3.) Plot model ----

plot(DecayFIT.RQ4, plot.type = "observations.vs.predictions" )
plot(DecayFIT.RQ4, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ4, plot.type = "parameters.vs.covariates", ask = TRUE) # <- use this in paper
npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ4) # skewness and kurtosis of normalized prediction discrepancies lower then in log model

# _4.) Compare models ----

compare.saemix(DecayFIT.RQ1, DecayFIT.RQ2, DecayFIT.RQ3, DecayFIT.RQ4) # third one is best, the second to last complex

get_p_from_seamix_lrt(DecayFIT.RQ1, DecayFIT.RQ2) # adding sex to null is sign. better BIC.cov
get_p_from_seamix_lrt(DecayFIT.RQ2, DecayFIT.RQ3) # adding litter to sex is sign. better BIC.cov
get_p_from_seamix_lrt(DecayFIT.RQ4, DecayFIT.RQ3) # adding diet does not improve fit

# _5.) Show coefficients ----

summary(DecayFIT.RQ4)

# ----------------------------------------------------
#   -----------------  Fixed effects  ------------------
#   ----------------------------------------------------
#               Parameter Estimate     SE  CV(%) p-value
#   1                   A   26.520 1.4402   5.43       -
#   2   beta_AnimalSex(A)    0.171 0.0156   9.15 0.00000 *
#   3  beta_LitterSize(A)   -0.012 0.0084  71.80 0.16371
#   4  beta_FatherDiet(A)   -0.055 0.0224  40.70 0.01400 *
#   5  beta_MotherDiet(A)   -0.046 0.0164  35.79 0.00521 *
#   6                   B    1.786 0.5525  30.93       -
#   7   beta_AnimalSex(B)    0.105 0.0846  80.99 0.21693
#   8  beta_LitterSize(B)   -0.085 0.0408  48.24 0.03819 *
#   9  beta_FatherDiet(B)    0.326 0.0925  28.40 0.00043 *
#  10  beta_MotherDiet(B)    0.071 0.0862 121.46 0.41034 
#  11                   k    0.057 0.0152  26.65       -
#  12   beta_AnimalSex(k)    0.044 0.0749 169.22 0.55456
#  13  beta_LitterSize(k)   -0.097 0.0389  40.05 0.01252 *
#  14  beta_FatherDiet(k)    0.386 0.0998  25.82 0.00011 *
#  15  beta_MotherDiet(k)    0.154 0.0772  50.08 0.04586 *
#  16                 a.1    0.753 0.0399   5.29       -

# _7.) Answer RQ4 ----

# The litter only model is the best fit to the data, diet doesn't add anything.

# _8) Plot complex last model model ----

# __a) Show values for plotting ----

coefficients(DecayFIT.RQ4)$fixed
summary(DecayFIT.RQ4)

# reference class for covariate AnimalSex :  f 
# reference class for covariate FatherDiet : CD 
# reference class for covariate MotherDiet : CD

# __b) Create canvas ----

plot(x = 1,                 
     xlab = "days [d]", 
     ylab = "weight [g]",
     xlim = c(0 , max(mice_f1_slct$MeasurementDay) + 5 ), 
     ylim = c(0 , max(mice_f1_slct$BodyWeightG) + 5),
     # main = "Significant model coefficients' offsets from reference values",
     type = "n")

# __c) Show raw data ----

# shows raw data - beyond this all are at reference level - not necessarily among the data
# points(x = mice_f1_slct$MeasurementDay, y = mice_f1_slct$BodyWeightG, col = "lightgrey", pch = 16)

# __d) Male mice growth ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 + 0.171)
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "black", add = TRUE)

# __e) Female mice growth ----

a = coefficients(DecayFIT.RQ4)$fixed[1]
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3]

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "black", add = TRUE)

# __f) Male mice growth - increased litter size by one----

a = coefficients(DecayFIT.RQ4)$fixed[1]               * (1 + 0.171)
b = coefficients(DecayFIT.RQ4)$fixed[2] * (1 - 0.085)
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 - 0.097)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "blue", add = TRUE)

# __g) Female mice growth - increased litter size by one ----

a = coefficients(DecayFIT.RQ4)$fixed[1] 
b = coefficients(DecayFIT.RQ4)$fixed[2] * (1 - 0.085)
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 - 0.097)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "blue", add = TRUE)

# __h) Male mice growth - mother diet ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 - 0.046) * (1 + 0.171)
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 + 0.154)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "darkorange", add = TRUE)

# __i) Female mice growth - mother diet ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 - 0.046)
b = coefficients(DecayFIT.RQ4)$fixed[2]
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 + 0.154)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "darkorange", add = TRUE)

# __j) Male mice growth - father diet ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 - 0.055) *  (1 + 0.171)
b = coefficients(DecayFIT.RQ4)$fixed[2] * (1 + 0.326)
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 + 0.386)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "red", add = TRUE)

# __k) Female mice growth  - father diet ----

a = coefficients(DecayFIT.RQ4)$fixed[1] * (1 - 0.055)
b = coefficients(DecayFIT.RQ4)$fixed[2] * (1 + 0.326)
k = coefficients(DecayFIT.RQ4)$fixed[3] * (1 + 0.386)

curve(a * (1 - b * exp(-k * x)), from = min(mice_f1_slct$MeasurementDay), to = max(mice_f1_slct$MeasurementDay),
      col = "red", add = TRUE)

# adjust and export  / "all effects are additive"
legend(33, 17, legend=c("sons and daughters, parents at control diet, \nno siblings", "one more sibling in litter", "mother on treatment diet", "father on treatment diet"),
       col=c("black", "blue", "darkorange", "red"), lty = c(1,1,1,1), cex=0.8)

# __l) Plot graphic ----

dev.print(svg, paste0(here("../manuscript/display_items"), "/", "015_r_use_saemix__model_rq4.svg"))
dev.print(pdf, paste0(here("../manuscript/display_items"), "/", "015_r_use_saemix__model_rq4.pdf"))

# RQ5: Can the complex model be improved by backward selection? ----

# # _1.) Define model object ----
# 
# DecayModel.RQ5 <- backward.procedure(DecayFIT.RQ4, trace = TRUE)
# 
# # _2.) Fit model ----
# 
# DecayFIT.RQ5 <- saemix(DecayModel.RQ5, ModelData.RQ1, saemix.options)
# 
# # _3.) Plot model ----
# 
# plot(DecayFIT.RQ5, plot.type="observations.vs.predictions" )
# plot(DecayFIT.RQ5, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE) # poor fit
# plot(DecayFIT.RQ5, plot.type="parameters.vs.covariates", ask = TRUE)
# DecayFIT.RQ5 <- npdeSaemix(DecayFIT.RQ5) # poor fit
# 
# # _4.) Compare models ----
# 
# compare.saemix(DecayFIT.RQ1, DecayFIT.RQ2, DecayFIT.RQ3, DecayFIT.RQ4, DecayFIT.RQ5)  # Last model is worst....
# get_p_from_seamix_lrt(DecayFIT.RQ5, DecayFIT.RQ4)         # ...and significantly so
 
 # _5.) Check estimates ----
# 
# summary(DecayFIT.RQ5)
# 
# # _6.) Answer RQ5 ----
# 
# # Backward selection did not find a better-fitting model

# RQ6: Get estimates for males only, for DEG discovery comparison ---- 

# _1.) Split data by sexes, in DEG only males are relevant  ----

ModelData.RQ6.male <- saemixData(
  name.data = {mice_f1_slct %>% filter(AnimalSex == "m")}, header = TRUE, name.group = c("AnimalId"), name.predictors = c("MeasurementDay"), name.response = c("BodyWeightG"), name.X = "MeasurementDay",
  name.covariates = c("LitterSize", "FatherDiet", "MotherDiet")
  )

# _2.) Define model object ----

# __a) Intercept only

DecayModel.RQ6.null <- saemixModel(model = decay.model,
                                   description= "Exponential approach", 
                                   psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                                   transform.par = c(1,1,1), 
                                   fixed.estim = c(1,1,1),
                                   covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                                   # covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                                   omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                                   error.model="constant")

# __b) Litter size

DecayModel.RQ6.litter <- saemixModel(model = decay.model,
                                   description= "Exponential approach", 
                                   psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                                   transform.par = c(1,1,1), 
                                   fixed.estim = c(1,1,1),
                                   covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                                   covariate.model = matrix(c(1,1,1), ncol=3, byrow=TRUE),
                                   omega.init = matrix(c(1,0,0,0, 1,0,0,0,1), ncol=3, byrow=TRUE), 
                                   error.model="constant")

# __c) Litter size and diets

DecayModel.RQ6.litter.diet <- saemixModel(model = decay.model,
                                   description= "Exponential approach", 
                                   psi0 = matrix( c(700,0.9,0.02, 0,0,0), ncol=3, byrow = TRUE, dimnames = list(NULL, c("A","B","k"))),
                                   transform.par = c(1,1,1), 
                                   fixed.estim = c(1,1,1),
                                   covariance.model= matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow = TRUE),
                                   covariate.model = matrix(c(1,1,1, 1,1,1, 1,1,1), ncol=3, byrow=TRUE),
                                   omega.init = matrix(c(1,0,0,0, 1,0,0,0,1),ncol=3, byrow=TRUE), 
                                   error.model="constant")
# _3.) Fit models ----

DecayFIT.RQ6.null <- saemix(DecayModel.RQ6.null, ModelData.RQ6.male, saemix.options)

DecayFIT.RQ6.litter <- saemix(DecayModel.RQ6.litter, ModelData.RQ6.male, saemix.options)

DecayFIT.RQ6.litter.diet <- saemix(DecayModel.RQ6.litter.diet, ModelData.RQ6.male, saemix.options)

# _4.) Evaluate and rank models ----

compare.saemix(DecayFIT.RQ6.null, DecayFIT.RQ6.litter, DecayFIT.RQ6.litter.diet)  # Second model appears most informaive ...

get_p_from_seamix_lrt(DecayFIT.RQ6.null, DecayFIT.RQ6.litter) # ... litter more informative then 0

get_p_from_seamix_lrt(DecayFIT.RQ6.litter, DecayFIT.RQ6.litter.diet) # ... adding diet doesn't improve model

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "015_r_use_saemix.Rdata"))
renv::snapshot()

