#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis testing - using {nlme}"
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

library("here")
library("dplyr")
library("ggplot2")

library("nlme")
library("performance") # Model inspection

# Setup data and model ----

# _1.) Read in data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Add litter size to f1 ----

mice_f1_slct <- left_join(mice_f1_slct, {mice_f0_slct %>% dplyr::select(AnimalId, MatingWith, LitterSize) %>% distinct}, by = c("MotherId" = "AnimalId", "FatherId" = "MatingWith"))

# _3.) Check data ----

mice_f1_slct %>% dplyr::select(MeasurementDay, BodyWeightG, AnimalId, AnimalSex) %>% arrange(AnimalSex, AnimalId, MeasurementDay) %>% print(n = Inf)

ggplot(data = mice_f1_slct, aes(x = "MeasurementDay", y="BodyWeightG",  group = "AnimalId")) +
  geom_line( mapping = aes( x = MeasurementDay, y = BodyWeightG, group = AnimalId, color = AnimalId)) + 
  facet_grid(. ~ AnimalSex) +
  theme_bw()

# _4.) Define possibly applicable model functions ----

# __a) Exponential approach as in {015_r_use_saemix.R}  ----

decay.formula <-  as.formula(y ~ a * (1 - b * exp( -k * x)))
    

# __b) Testing function to get starting values ----

# from estimated previous values - equivalent model (RQ1)
a = 25.3005
b = 1.2644
k = 0.0392
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkgray", xlim =c(30, 100), ylim=c(15, 25))

# RQ1: Fit exponential approach to data to get a null model for comparison ----

# __a) Get null model as in RQ1 {015_r_use_saemix.R} ----

exp.appr.fit.null <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                          data = mice_f1_slct,
                          fixed  = a + b + k ~ 1,
                          random = a ~ 1,
                          groups = ~ AnimalId,
                          start = c(22.41907, 2.87427, 0.07869),
                          na.action = na.exclude,
                          control = nlmeControl(maxIter = 300, msVerbose = FALSE))

summary(exp.appr.fit.null) 

#      AIC      BIC    logLik
# 937.2057 955.7246 -463.6028

# __b) Plot null model ----

# null model coefficients of exp.appr.fit.null
a = 25.302337
b = 1.312241
k = 0.040536

curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "red", xlim =c(30, 100), ylim=c(15, 25),  add = TRUE)



# RQ2: Does Sex have an association with the total weight gain? ----

# __a) Get null model as in RQ1 {015_r_use_saemix.R} ----

# https://stats.stackexchange.com/questions/536364/help-understanding-fixed-effects-interaction-terms-in-nlme


exp.appr.fit.sex <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                         data = mice_f1_slct,
                         fixed  = a + b + k ~ AnimalSex,
                         random = a ~ 1,
                         groups = ~ AnimalId,
                         na.action = na.exclude,
                         start = c(25.30,  1.31,  0.040,
                                    0.17,  0.09,   0.04),
                         control = nlmeControl(maxIter = 300, msVerbose = FALSE))

summary(exp.appr.fit.sex) 

# AIC      BIC    logLik
# 879.7389 909.3691 -431.8694 - better

plot(exp.appr.fit.sex)

# __b) Plot null and sex model ----

# null model coefficients of exp.appr.fit.null
a = 25.302337
b = 1.312241
k = 0.040536

curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkgrey", xlim =c(30, 100), ylim=c(15, 25))

# females in exp.appr.fit.sex

a = 22.6461642335
b = 1.2684247933
k = 0.0408665699

curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "black", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

a = a + 4.5663630884 # <- the only significant change according to model output
b = b + 0.0640978693
k = k - 0.0005498267

curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "red", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# __c) Compare null and sex model ----

anova(exp.appr.fit.null,exp.appr.fit.sex) # <- sex model is better



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

# _4.) Plot model ----

plot(DecayFIT.RQ4, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ4, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ4, plot.type="parameters.vs.covariates", ask=TRUE)
npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ4) # skewness and kurtosis of normalised prediction discrepancies lower then in log model


# _5.) Compare models ----

DecayFIT.RQ3
DecayFIT.RQ4 # <- BIC AIC with diets worse then sex and litter size alone

teststatRQ12 <- -2 * (as.numeric(logLik(DecayFIT.RQ3)) - as.numeric(logLik(DecayFIT.RQ4)))
p.val <- pchisq(teststatRQ12, df = 3, lower.tail = FALSE)
p.val # model with diets is not distinctly differnet from when not adding diets


# RQ5: What are the effects of diet within each sex and, and dependent on litter sizes? ----

# _1.) Define data with sex covariate ----

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

# _3.) Fit model ----

DecayFIT.RQ5.male <- saemix(DecayModel.RQ5, ModelData.RQ5.male, saemix.options)

# Fixed effects
# Parameter          Estimate    SE    CV(%)  p-value
# A                  28.566939 2.5112    8.79 -      
#   beta_LitterSize(A)  0.000838 0.0184 2199.03 0.9637 
#   beta_FatherDiet(A) -0.035472 0.0285   80.24 0.2127 
#   beta_MotherDiet(A) -0.053695 0.0285   53.04 0.0594 
# B                   1.906307 0.7931   41.60 -      
#   beta_LitterSize(B) -0.090955 0.0844   92.83 0.2814 
#   beta_FatherDiet(B)  0.012596 0.1280 1016.57 0.9216 
#   beta_MotherDiet(B)  0.059983 0.1354  225.69 0.6577 
# k                   0.045201 0.0160   35.44 -      
#   beta_LitterSize(k) -0.059031 0.0739  125.14 0.4242 
#   beta_FatherDiet(k)  0.060994 0.1139  186.66 0.5921 
#   beta_MotherDiet(k)  0.173567 0.1191   68.64 0.1451 * 
# a.1                 0.834962 0.0589    7.05 -      


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


DecayFIT.RQ5.female <- saemix(DecayModel.RQ5, ModelData.RQ5.female, saemix.options)

# Fixed effects
# Parameter          Estimate   SE    CV(%) p-value
# A                  22.41907 1.2216   5.45 -      
#   beta_LitterSize(A)  0.00652 0.0116 178.40 0.57510
#   beta_FatherDiet(A) -0.02591 0.0204  78.61 0.20333
#   beta_MotherDiet(A) -0.00363 0.0194 535.06 0.85174
# B                   2.87427 0.9309  32.39 -      
#   beta_LitterSize(B) -0.13740 0.0688  50.09 0.04590 * 
#   beta_FatherDiet(B) -0.03240 0.1159 357.75 0.77984
#   beta_MotherDiet(B) -0.29874 0.1095  36.65 0.00637 * 
# k                   0.07869 0.0207  26.27 -      
#   beta_LitterSize(k) -0.11206 0.0558  49.81 0.04469
#   beta_FatherDiet(k) -0.04751 0.1003 211.08 0.63567
#   beta_MotherDiet(k) -0.21824 0.0945  43.28 0.02086 * 
# a.1                 0.52514 0.0433   8.25 - 



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


# _4.) Plot model ----

plot(DecayFIT.RQ5.female, plot.type="observations.vs.predictions" )
plot(DecayFIT.RQ5.male, plot.type="observations.vs.predictions" )

plot(DecayFIT.RQ5.female, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)
plot(DecayFIT.RQ5.male, plot.type = "both.fit",  ilist = 1:9, smooth = TRUE)

plot(DecayFIT.RQ5.female, plot.type="parameters.vs.covariates", ask=TRUE)
plot(DecayFIT.RQ5.male, plot.type="parameters.vs.covariates", ask=TRUE)

npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ5.female) # skewness and kurtosis of normalised prediction discrepancies lower then in log model
npde.DecayFIT.RQ4 <- npdeSaemix(DecayFIT.RQ5.male) # skewness and kurtosis of normalised prediction discrepancies lower then in log model

# _5.) Compare models ----

# no model for comparison built yet, needed would be above data sets without diets

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "017_r_use_nlme.Rdata"))
renv::snapshot()

