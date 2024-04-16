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

# RQ1: What are the effects of diet regardless of sex ? ----

# _1.) Get a suitable null model which is comparable to the subsequent one ----

exp.appr.fit.diet.nosex.null <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                                     data = mice_f1_slct,
                                     fixed  = a + b + k ~  1,
                                     random = a  ~  AnimalSex | AnimalId,
                                     na.action = na.exclude,
                                     start = c(25.30,  1.31,  0.04),
                                     control = nlmeControl(msMaxIter = 1000, msVerbose = TRUE))

# _2.) Check null model ----

summary(exp.appr.fit.diet.nosex.null)
# AIC      BIC    logLik
# 932.4092 958.3356 -459.2046

# _3.) Plot null model residuals ----

plot(exp.appr.fit.diet.nosex.null) # residuals seem ok

# _4.) Get diet model as in RQ4 in `015_r_use_saemix.R` but without litter size ----

# Diet interactions were not significant

exp.appr.fit.diet.nosex <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                                data = mice_f1_slct,
                                fixed  = a + b + k ~  FatherDiet + MotherDiet,
                                random = a  ~ AnimalSex | AnimalId,
                                na.action = na.exclude,
                                start = c(25.30,  1.31,  0.04,
                                          0.17,  0.09,  0.04,
                                          0.01,  0.08,  0.02),
                                control = nlmeControl(msMaxIter = 1000, msVerbose = TRUE))

# _5.) Check diet model ----

summary(exp.appr.fit.diet.nosex)
# AIC      BIC    logLik
# 924.2106 972.3597 -449.1053

# _6.) Plot diet model residuals ----

plot(exp.appr.fit.diet.nosex) # residuals seem ok

# _7.) Compare models ----

anova(exp.appr.fit.diet.nosex.null, exp.appr.fit.diet.nosex) # adding diet imporves model

# _8.) Plot model predictions ----

# __a) Curve plots ----

# null model curve - all mice ragrdless of sex
a = fixed.effects(exp.appr.fit.diet.nosex.null)[1]
b = fixed.effects(exp.appr.fit.diet.nosex.null)[2]
k = fixed.effects(exp.appr.fit.diet.nosex.null)[3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkgray", xlim =c(30, 100), ylim=c(15, 25))


# chow diet curve  - all mice ragrdless of sex
a = fixed.effects(exp.appr.fit.diet.nosex)[1]
b = fixed.effects(exp.appr.fit.diet.nosex)[4]
k = fixed.effects(exp.appr.fit.diet.nosex)[7]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "black", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# father hfd diet curve  - all mice ragrdless of sex

af = a + fixed.effects(exp.appr.fit.diet.nosex)[2] # * 
bf = b + fixed.effects(exp.appr.fit.diet.nosex)[5]
kf = k + fixed.effects(exp.appr.fit.diet.nosex)[8]
curve(af * (1 - bf * exp( -kf * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkorange", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# mother hfd diet curve  - all mice ragrdless of sex
af = a + fixed.effects(exp.appr.fit.diet.nosex)[3] # * 
bf = b + fixed.effects(exp.appr.fit.diet.nosex)[6]
kf = k + fixed.effects(exp.appr.fit.diet.nosex)[9]
curve(af * (1 - bf * exp( -kf * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "red", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# __b) Prediction plots ----

# RQ2: Get a null model to investigate the sex specific effcet more ----

# _1.) Get null model as in RQ1 of `015_r_use_saemix.R` ----

exp.appr.fit.null <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                          data = mice_f1_slct,
                          fixed  = a + b + k ~ 1,
                          random = a ~ 1 | AnimalId,
                          start = c(22.41907, 2.87427, 0.07869),
                          na.action = na.exclude,
                          control = nlmeControl(maxIter = 300, msVerbose = FALSE))

summary(exp.appr.fit.null) 

#      AIC      BIC    logLik
# 937.2057 955.7246 -463.6028

# RQ3: Does Sex have an association with the total weight gain? ----

# _1.) Get sex model as in RQ1 of `015_r_use_saemix.R` ----

# https://stats.stackexchange.com/questions/536364/help-understanding-fixed-effects-interaction-terms-in-nlme

exp.appr.fit.sex <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                         data = mice_f1_slct,
                         fixed  = a + b + k ~ AnimalSex,
                         random = a ~ 1 | AnimalId,
                         na.action = na.exclude,
                         start = c(25.30,  1.31,   0.04,
                                    0.17,  0.09,   0.04),
                         control = nlmeControl(msMaxIter = 1000, msVerbose = FALSE))

# _2.) Check sex model ----

summary(exp.appr.fit.sex) 

# AIC      BIC    logLik
# 879.7389 909.3691 -431.8694 - better

# _3.) Plot sex model residuals ----

plot(exp.appr.fit.sex)

# _3.) Compare null and sex model ----

anova(exp.appr.fit.null,exp.appr.fit.sex) # <- sex model is better

# _4.) Plot null and sex model ----

# __a) Curve plots ----

# null model coefficients of exp.appr.fit.null
a = fixef(exp.appr.fit.null)[1]
b = fixef(exp.appr.fit.null)[2]
k = fixef(exp.appr.fit.null)[3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkgrey", xlim =c(30, 100), ylim=c(15, 25))

# females in exp.appr.fit.sex
a = fixef(exp.appr.fit.sex)[1]
b = fixef(exp.appr.fit.sex)[3]
k = fixef(exp.appr.fit.sex)[5]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "black", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# male in exp.appr.fit.sex
a = a + fixef(exp.appr.fit.sex)[2]
b = b + fixef(exp.appr.fit.sex)[4]
k = k + fixef(exp.appr.fit.sex)[6]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "red", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# RQ4: Do Sex and litter size associate with the total weight gain? ----

# _1.) Get litter model as in RQ3 of `015_r_use_saemix.R` ----


# exp.appr.fit.litter <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
#                          data = mice_f1_slct,
#                          fixed  = a + b + k ~ AnimalSex + LitterSize,
#                          random = a  ~ 1 | AnimalId,
#                          na.action = na.exclude,
#                          start = c(25.30,  1.31,  0.04,
#                                     0.17,  0.09,   0.04,
#                                     0.01,  0.08,   0.02),
#                          control = nlmeControl(msMaxIter = 1000, msVerbose = TRUE))

# summary(exp.appr.fit.litter) 

# Effect of litter size is to small to be estimated with nlme package and was insignificant in {seamix}

# RQ5: What are the effects of diet within each sex ? ----

# _1.) Get sex/diet model as in RQ4 in `015_r_use_saemix.R` but without litter size ----

exp.appr.fit.diet <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                          data = mice_f1_slct,
                          fixed  = a + b + k ~ AnimalSex + FatherDiet + MotherDiet,
                          random = a  ~ 1 | AnimalId,
                          na.action = na.exclude,
                          start = c(25.30,  1.31,  0.04,
                                     0.17,  0.09,  0.04,
                                     0.01,  0.08,  0.02,
                                     0.01,  0.08,  0.02),
                          control = nlmeControl(msMaxIter = 1000, msVerbose = TRUE))
 
# _2.) Check diet model----

summary(exp.appr.fit.diet) 

# AIC      BIC    logLik
# 873.7234 925.5763 -422.8617 <- better 

# _3.) Plot diet model resiudals ----

plot(exp.appr.fit.diet)

# _4.) Compare sex and diet model ----

anova(exp.appr.fit.sex,exp.appr.fit.diet) # <- diet model is better

# _5.) Plot sex/diet model ----

# __a) Curve plots ----

# chow diet females 
a = fixed.effects(exp.appr.fit.diet)[1]
b = fixed.effects(exp.appr.fit.diet)[5]
k = fixed.effects(exp.appr.fit.diet)[9]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkgrey", xlim =c(30, 100), ylim=c(15, 25))

# chow diet males
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+1]
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+1]
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+1]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "lightgrey", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# **continue here** ----

# females and father diet

afd = af - 0.927409
bfd = bf + 0.038206
kfd = kf - 0.002177

curve(afd * (1 - bfd * exp( -kfd * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "orange", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# females  and mother diet

afd = af - 0.926447
bfd = bf - 0.005355
kfd = kf - 0.003048

curve(afd * (1 - bfd * exp( -kfd * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkorange", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# males in exp.appr.fit.diet

am = af + 4.491561 # <- the only significant change according to model output
bm = bf + 0.065160
km = kf - 0.000464

curve(am * (1 - bm * exp( -km * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "lightgrey", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# males in and father diet

amd = am - 0.927409
bmd = bm + 0.038206
kmd = km - 0.002177

curve(amd * (1 - bmd * exp( -kmd * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "red", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# males in and mother diet

amd = am - 0.926447
bmd = bm - 0.005355
kmd = km - 0.003048

curve(amd * (1 - bmd * exp( -kmd * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", col = "darkred", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# __b) Prediction plots ----

# not yet implemented


# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "017_r_use_nlme.Rdata"))
renv::snapshot()

