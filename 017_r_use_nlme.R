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

library("here")    # handle path names
library("dplyr")   # handle data more easily

library("lattice") # create trellis graphs
library("ggpubr")  # save trellis graphs
library("ggplot2") # save trellis graphs

library("nlme")    # model non-linear mixed-effect dependencies

# Setup data and model ----

# _1.) Read in data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Add litter size to f1 ----

mice_f1_slct <- left_join(mice_f1_slct, {mice_f0_slct %>% dplyr::select(AnimalId, MatingWith, LitterSize) %>% distinct}, by = c("MotherId" = "AnimalId", "FatherId" = "MatingWith"))

# _3.) Check data ----

plot_data_check <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, groups = AnimalSex, data = mice_f1_slct, xlab = "day [d]", ylab = "body weight [g]", auto.key = list(title = "sex"))
plot_data_check

ggsave("017_r_use_nlme__data_check.pdf", plot = ggarrange(plot_data_check), path = here("../manuscript/display_items"),
  scale = 1, width = 9, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# _4.) Define possibly applicable model functions ----

# __a) Exponential approach as in "015_r_use_saemix.R" ----

approach.formula <- as.formula(y ~ a * (1 - b * exp( -k * x)))
    
# __b) Testing function to get starting values ----

# from estimated previous values - equivalent model (RQ1)

a = 25.3005; b = 1.2644; k = 0.0392
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab = "day [d]", ylab = "body weight [g]", col = "darkgray", xlim =c(30, 100), ylim=c(15, 25))

# RQ1: What are the effects of diet regardless of sex ? ----

# _1.) Get a suitable null model ----

# __a) Intended model:

plot_rq1_null_model <- xyplot(BodyWeightG ~ MeasurementDay | AnimalSex, groups = AnimalId , data = mice_f1_slct, xlab = "day [d]", ylab = "body weight [g]",
                              panel = function(x, y) {
                                panel.xyplot(x, y)
                                panel.loess(x, y)
                                })

ggsave("017_r_use_nlme__rq1_null_model.pdf", plot = ggarrange(plot_rq1_null_model), path = here("../manuscript/display_items"),
       scale = 1, width = 5, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# __a) Build model:

exp.appr.fit.diet.nosex.null <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                                     data = mice_f1_slct,
                                     fixed  = a + b + k ~  1,
                                     random = a  ~  AnimalSex | AnimalId,
                                     na.action = na.exclude,
                                     start = c(25.30,  1.31,  0.04),
                                     control = nlmeControl(msMaxIter = 50, msVerbose = FALSE))

# _2.) Check null model ----

summary(exp.appr.fit.diet.nosex.null)

# AIC      BIC    logLik
# 932.4092 958.3356 -459.2046

# Fixed effects:  a + b + k ~ 1 
# Value  Std.Error  DF  t-value p-value
# a.(Intercept) 23.112100 0.25354413 248 91.15613       0
# b              1.311560 0.06272352 248 20.91018       0
# k              0.040515 0.00176842 248 22.91022       0

# _3.) Plot null model fits  ----

plot_rq1_null_model_fit <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, fit = exp.appr.fit.diet.nosex.null,
       strip = TRUE, aspect = "xy", grid = TRUE,
       panel = function(x, y, ..., fit, subscripts) {
         panel.xyplot(x, y, ...)
         ypred <- fitted(fit)[subscripts]
         panel.lines(x, ypred, col = "black")
       },
       xlab = "day [d]", ylab = "body weight [g]")

plot_rq1_null_model_fit

ggsave("017_r_use_nlme__plot_rq1_null_model_fit.pdf", plot = ggarrange(plot_rq1_null_model_fit), path = here("../manuscript/display_items"),
       scale = 1, width = 12, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# _4.) Plot null model residuals ----

plot(exp.appr.fit.diet.nosex.null) # residuals seem ok

# ** continue here **

# _5.) Get matching diet model as in RQ4 in `015_r_use_saemix.R` but without litter size ----

# (Diet interactions were not significant)

exp.appr.fit.diet.nosex <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                                data = mice_f1_slct,
                                fixed  = a + b + k ~  FatherDiet + MotherDiet,
                                random = a  ~  AnimalSex | AnimalId,
                                na.action = na.exclude,
                                start = c(25.30,  1.31,  0.04,
                                          0.17,  0.09,  0.04,
                                          0.01,  0.08,  0.02),
                                control = nlmeControl(msMaxIter = 50, msVerbose = FALSE))

# _5.) Check diet model ----

summary(exp.appr.fit.diet.nosex)

# AIC      BIC    logLik
# 924.2106 972.3597 -449.1053 <- better then null model - see LRT below

# Fixed effects:  a + b + k ~ FatherDiet + MotherDiet 
# Value Std.Error  DF  t-value p-value
# a.(Intercept)   24.262459 0.4877282 242 49.74586  0.0000 *
# a.FatherDietHFD -1.037515 0.4913831 242 -2.11142  0.0358
# a.MotherDietHFD -0.943602 0.4729634 242 -1.99508  0.0472 * only a.MotherDietHFD significant
# b.(Intercept)    1.329206 0.1167432 242 11.38572  0.0000
# b.FatherDietHFD  0.001529 0.1247160 242  0.01226  0.9902
# b.MotherDietHFD -0.030819 0.1234833 242 -0.24958  0.8031
# k.(Intercept)    0.038565 0.0032705 242 11.79192  0.0000
# k.FatherDietHFD  0.001229 0.0035098 242  0.35004  0.7266
# k.MotherDietHFD  0.002341 0.0034729 242  0.67419  0.5008

# _6.) Plot diet model fits  ----

plot_rq1_diet_model_fit <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, fit = exp.appr.fit.diet.nosex,
                                  strip = TRUE, aspect = "xy", grid = TRUE,
                                  panel = function(x, y, ..., fit, subscripts) {
                                    panel.xyplot(x, y, ...)
                                    ypred <- fitted(fit)[subscripts]
                                    panel.lines(x, ypred, col = "black")
                                  },
                                  xlab = "day [d]", ylab = "body weight [g]")

plot_rq1_diet_model_fit

ggsave("017_r_use_nlme__plot_rq1_diet_model_fit.pdf", plot = ggarrange(plot_rq1_diet_model_fit), path = here("../manuscript/display_items"),
       scale = 1, width = 12, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# _7.) Plot diet model residuals ----

plot(exp.appr.fit.diet.nosex) # residuals seem ok

# _8.) Compare models ----

anova(exp.appr.fit.diet.nosex.null, exp.appr.fit.diet.nosex) # adding diet improves model

# _9.) Plot diet model predictions ----

# null model curve - all mice regardless of sex
a = fixed.effects(exp.appr.fit.diet.nosex.null)[1]
b = fixed.effects(exp.appr.fit.diet.nosex.null)[2]
k = fixed.effects(exp.appr.fit.diet.nosex.null)[3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab = "day [d]", ylab = "body weight [g]",
      col = "black", xlim =c(30, 100), ylim=c(15, 25), lty = "dashed")

# chow diet curve  - all mice regardless of sex
a = fixed.effects(exp.appr.fit.diet.nosex)[1]
b = fixed.effects(exp.appr.fit.diet.nosex)[4]
k = fixed.effects(exp.appr.fit.diet.nosex)[7]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab = "day [d]", ylab = "body weight [g]", 
      col = "black", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# father hfd diet curve  - all mice regardless of sex
af = a + fixed.effects(exp.appr.fit.diet.nosex)[2] 
bf = b + fixed.effects(exp.appr.fit.diet.nosex)[5]
kf = k + fixed.effects(exp.appr.fit.diet.nosex)[8]
curve(af * (1 - bf * exp( -kf * x)), from = 35, to = 100, xlab = "day [d]", ylab = "body weight [g]", 
      col = "red", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dashed")

# mother hfd diet curve  - all mice regardless of sex
af = a + fixed.effects(exp.appr.fit.diet.nosex)[3] # * 
bf = b + fixed.effects(exp.appr.fit.diet.nosex)[6]
kf = k + fixed.effects(exp.appr.fit.diet.nosex)[9]
curve(af * (1 - bf * exp( -kf * x)), from = 35, to = 100, 
      xlab = "day [d]", ylab = "body weight [g]", col = "red", xlim =c(30, 100), ylim=c(15, 25),
      add = TRUE)

legend(67, 19, legend=c("null model", "low-caloric", "father high-caloric", "mother high-caloric"),
       col=c("black", "black", "red", "red"), lty = c(2,1,2,1), cex=0.8)

# RQ2: What is the sex specific effect on body weight over time? ----

# _1.) Get null model as in RQ1 of `015_r_use_saemix.R` ----

# eliminating sex from the random effect structure

exp.appr.fit.null <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                          data = mice_f1_slct,
                          fixed  = a + b + k ~ 1,
                          random = a ~ 1 | AnimalId,
                          start = c(22.41907, 2.87427, 0.07869),
                          na.action = na.exclude,
                          control = nlmeControl(maxIter = 50, msVerbose = FALSE))

# _2.) Check null model ----

summary(exp.appr.fit.null) 

#      AIC      BIC    logLik
# 937.2057 955.7246 -463.6028

# _3.) Get sex model as in RQ1 of `015_r_use_saemix.R` ----

# https://stats.stackexchange.com/questions/536364/help-understanding-fixed-effects-interaction-terms-in-nlme

# adding sex as fixed effect

exp.appr.fit.sex <- nlme(BodyWeightG ~ a * (1 - b * exp( -k * MeasurementDay)),
                         data = mice_f1_slct,
                         fixed  = a + b + k ~ AnimalSex,
                         random = a ~ 1 | AnimalId,
                         na.action = na.exclude,
                         start = c(25.30,  1.31,   0.04,
                                    0.17,  0.09,   0.04),
                         control = nlmeControl(msMaxIter = 1000, msVerbose = FALSE))

# _4.) Check sex model ----

summary(exp.appr.fit.sex) 

# AIC      BIC    logLik
# 879.7389 909.3691 -431.8694 < - better then matching null model

# Fixed effects:  a + b + k ~ AnimalSex 
# Value Std.Error  DF  t-value p-value
# a.(Intercept) 22.646164 0.3937909 245 57.50809  0.0000
# a.AnimalSexm   4.566363 0.5185962 245  8.80524  0.0000 * 
# b.(Intercept)  1.268425 0.1087662 245 11.66194  0.0000
# b.AnimalSexm   0.064098 0.1325391 245  0.48361  0.6291
# k.(Intercept)  0.040867 0.0031790 245 12.85519  0.0000
# k.AnimalSexm  -0.000550 0.0038096 245 -0.14433  0.8854

# _5.) Plot sex model residuals ----

plot(exp.appr.fit.sex)

# _6.) Compare null and sex model ----

anova(exp.appr.fit.null, exp.appr.fit.sex) # <- sex model is then matching null model better

# _7.) Plot null and sex model ----

# __a) Trellis plots ----

plot_rq1_sex_model_fit <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, fit = exp.appr.fit.sex,
                                  strip = TRUE, aspect = "xy", grid = TRUE,
                                  panel = function(x, y, ..., fit, subscripts) {
                                    panel.xyplot(x, y, ...)
                                    ypred <- fitted(fit)[subscripts]
                                    panel.lines(x, ypred, col = "black")
                                  },
                                  xlab = "day [d]", ylab = "body weight [g]")

plot_rq1_sex_model_fit

ggsave("017_r_use_nlme__plot_rq2_sex_model_fit.pdf", plot = ggarrange(plot_rq1_sex_model_fit), path = here("../manuscript/display_items"),
       scale = 1, width = 12, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# __b) Curve plots ----

# null model coefficients of exp.appr.fit.null
a = fixef(exp.appr.fit.null)[1]
b = fixef(exp.appr.fit.null)[2]
k = fixef(exp.appr.fit.null)[3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "black", xlim =c(30, 100), ylim=c(15, 25), lty = "dashed")

# females in exp.appr.fit.sex
a = fixef(exp.appr.fit.sex)[1]
b = fixef(exp.appr.fit.sex)[3]
k = fixef(exp.appr.fit.sex)[5]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "red", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# male in exp.appr.fit.sex
a = a + fixef(exp.appr.fit.sex)[2]
b = b + fixef(exp.appr.fit.sex)[4]
k = k + fixef(exp.appr.fit.sex)[6]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "black", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

legend(67, 19, legend=c("null model", "male offsrping", "female offsrping"),
       col=c("black", "black", "red"), lty = c(2,1,1), cex=0.8)


# RQ3: Do Sex and litter size associate with the total weight gain? ----

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

# RQ4: What are the effects of diet, within each sex, on body weight ? ----

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
                          control = nlmeControl(msMaxIter = 50, msVerbose = FALSE))
 
# _2.) Check diet model ----

summary(exp.appr.fit.diet) 

# AIC      BIC    logLik
# 873.7234 925.5763 -422.8617 <- better 

# Fixed effects:  a + b + k ~ AnimalSex + FatherDiet + MotherDiet 
# Value Std.Error  DF  t-value p-value
# a.(Intercept)   23.781456 0.6059554 239 39.24622  0.0000
# a.AnimalSexm     4.491561 0.5067669 239  8.86317  0.0000
# a.FatherDietHFD -0.927409 0.5172474 239 -1.79297  0.0742 . - probably not significant
# a.MotherDietHFD -0.926447 0.5033514 239 -1.84056  0.0669 . - probably not significant
# b.(Intercept)    1.251589 0.1459295 239  8.57667  0.0000
# b.AnimalSexm     0.065160 0.1298100 239  0.50197  0.6162
# b.FatherDietHFD  0.038206 0.1220970 239  0.31292  0.7546
# b.MotherDietHFD -0.005355 0.1210320 239 -0.04424  0.9647
# k.(Intercept)    0.037969 0.0042242 239  8.98847  0.0000
# k.AnimalSexm    -0.000464 0.0037177 239 -0.12492  0.9007
# k.FatherDietHFD  0.002177 0.0034452 239  0.63182  0.5281
# k.MotherDietHFD  0.003048 0.0033996 239  0.89669  0.3708

# _3.) Plot diet model resiudals ----

plot(exp.appr.fit.diet)

# _4.) Compare sex and diet model ----

anova(exp.appr.fit.sex, exp.appr.fit.diet) # <- diet model is better

# _5.) Plot sex/diet model ----

# __a) Trellis plot of model fit

plot_rq4_sex__diet_model_fit <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, fit = exp.appr.fit.diet,
                                 strip = TRUE, aspect = "xy", grid = TRUE,
                                 panel = function(x, y, ..., fit, subscripts) {
                                   panel.xyplot(x, y, ...)
                                   ypred <- fitted(fit)[subscripts]
                                   panel.lines(x, ypred, col = "black")
                                 },
                                 xlab = "day [d]", ylab = "body weight [g]")

plot_rq4_sex__diet_model_fit

ggsave("017_r_use_nlme__plot_rq4_sex__diet_model_fit.pdf", plot = ggarrange(plot_rq4_sex__diet_model_fit), path = here("../manuscript/display_items"),
       scale = 1, width = 12, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# __b) Curve plots ----

# chow diet females 
a = fixed.effects(exp.appr.fit.diet)[1]
b = fixed.effects(exp.appr.fit.diet)[5]
k = fixed.effects(exp.appr.fit.diet)[9]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "red", xlim =c(30, 100), ylim=c(15, 25))

# fathers hfd diet females 
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+2] 
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+2] 
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+2]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "darkorange", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dashed")

# mothers hfd diet females 
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+3] 
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+3] 
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "darkorange", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dotted")

# mothers and fathers hfd diet females 
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+3] + fixed.effects(exp.appr.fit.diet)[1+2]
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+3] + fixed.effects(exp.appr.fit.diet)[5+2] 
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+3] + fixed.effects(exp.appr.fit.diet)[9+2]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "darkorange", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dotdash")

# chow diet males
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+1]
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+1]
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+1]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG",
      col = "black", xlim =c(30, 100), ylim=c(15, 25), add = TRUE)

# fathers hfd diet males
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+1] + fixed.effects(exp.appr.fit.diet)[1+2] 
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+1] + fixed.effects(exp.appr.fit.diet)[5+2] 
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+1] + fixed.effects(exp.appr.fit.diet)[9+2]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "darkgray", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dashed")

# mothers hfd diet males
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+1] + fixed.effects(exp.appr.fit.diet)[1+3] 
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+1] + fixed.effects(exp.appr.fit.diet)[5+3] 
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+1] + fixed.effects(exp.appr.fit.diet)[9+3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "darkgray", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dotted")

# mothers and fathers hfd diet males
a = fixed.effects(exp.appr.fit.diet)[1] + fixed.effects(exp.appr.fit.diet)[1+1] + fixed.effects(exp.appr.fit.diet)[1+2] + fixed.effects(exp.appr.fit.diet)[1+3]
b = fixed.effects(exp.appr.fit.diet)[5] + fixed.effects(exp.appr.fit.diet)[5+1] + fixed.effects(exp.appr.fit.diet)[5+2] + fixed.effects(exp.appr.fit.diet)[5+3] 
k = fixed.effects(exp.appr.fit.diet)[9] + fixed.effects(exp.appr.fit.diet)[9+1] + fixed.effects(exp.appr.fit.diet)[9+2] + fixed.effects(exp.appr.fit.diet)[9+3]
curve(a * (1 - b * exp( -k * x)), from = 35, to = 100, xlab="MeasurementDay", ylab="BodyWeightG", 
      col = "darkgray", xlim =c(30, 100), ylim=c(15, 25), add = TRUE, lty = "dotdash")

legend(62, 20.3, legend=c("sons: parents low-caloric", "father high-caloric", "mother high-caloric", "mother & father high-caloric",
                        "daughters: parents low caloric",  "father high-caloric", "mother high-caloric", "mother & father high-caloric"),
       col=c("black", "darkgrey", "darkgrey","darkgrey",
             "red", "darkorange", "darkorange", "darkorange"), 
       lty = c("solid", "dashed", "dotted", "dotdash",
               "solid", "dashed", "dotted", "dotdash"), cex=0.7)

ggsave("017_r_use_nlme__foo.pdf", path = here("../manuscript/display_items"),
       scale = 1, width = 12, height = 5, units = c("in"), dpi = 300, limitsize = TRUE)

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "017_r_use_nlme.Rdata"))
renv::snapshot()

