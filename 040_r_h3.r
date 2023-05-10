#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis 3"
#' author: "Paul Czechowski ``paul.czechowski@helmholtz-munich.de``"
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

# Prepare environment  ---- 

# _1.) Collect garbage ----

rm(list=ls())
gc()

# _2.) Packages ----

library("here")     # environment management - use package "here" in conjunction with a RStudio project
library("renv")     # environment can be snap shot

library("readxl")   # table I/O
library("janitor")  # data cleaning 

library("dplyr")    # pipes and more
library("magrittr") # even more pipes

library("plotrix")
library("lattice")
library("ggplot2")
library("ggrepel")  

library("modeest")

library("lme4")     # Linear Mixed Effects Models
library("lmerTest") # Tests in Linear Mixed Effects Models
library("mgcv")     # General Additive Model

library("effects")     # Model inspection
library("performance") # Model inspection
library("cAIC4")       # Model selection 

library("gtsummary")
library("modelsummary")

# 3.) Functions ----

get_model_prdictions_with_mgcv = function(model_fit, model_data, ...){
  
  
  # get predictions from merTools
  model_preditcions <- predict.gam(object = model_fit, newdata = model_data,
                                   se.fit = TRUE, terms = NULL,
                                   exclude=NULL,
                                   na.action = na.omit, unconditional = TRUE,
  )
  # diagnostic
  message("Dim. of model pred.= ", dim(model_preditcions)[1], " ", dim(model_preditcions)[2])
  message("Dim. of input data = ", dim(model_data)[1], " ", dim(model_data)[2] )
  
  # merge measured and predicted values for plotting
  model_data_with_predictions <- cbind(model_data, model_preditcions)
  
  return(model_data_with_predictions)
}



# Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_H2variables.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_H2variables.rds"))


# Logistic Regression: Select and shape data ----

# Hypotheses are
# - "Obesity of (fe)male offspring is dependent on both parents obesity"
# - "Obesity of (fe)male offspring is dependent on mothers' obesity"
# - "Obesity of (fe)male offspring is dependent on fathers' obesity"

# _1.) Isolate data for modelling ----

mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, AnimalSex, ObesityLgcl, ObeseParents) 

# _2.) Check balance of modelling data and get a graphical or table summary  ----

mice_f1_model_data %>% select(ObesityLgcl, ObeseParents) %>% count(ObesityLgcl, ObeseParents, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, ObeseParents) %>% table()

# Logistic Regression: Model offspring' obesity as function of parents obesity  ----

# _1.) Intercept-only model ----

mod_0 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
summary(mod_0)

# _2.) Animal sex only ----

mod_1 <- lme4::glmer(ObesityLgcl ~ AnimalSex +  (1 | AnimalId), data = mice_f1_model_data, family = binomial)
summary(mod_1)  # AnimalSexm *** (due to weight curve Obesity definistion)

# _3.) Parents obesity only  ----

mod_2 <- lme4::glmer(ObesityLgcl ~ ObeseParents + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
summary(mod_2) # ObeseParentsMotherFatherObese ***

# _4.) Full model ----

mod_3 <- lme4::glmer(ObesityLgcl ~ AnimalSex + ObeseParents + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
summary(mod_3)

exp(fixef(mod_3))
round(exp(fixef(mod_3)), digits = 3)

# _5.) Test logistic rgressions' fixed effects using ANOVA function and cAIC ----

# cAIC(mod_0) -- 0.34
# cAIC(mod_1) -- 0.17
# cAIC(mod_2) -- 0.30
# cAIC(mod_3) -- 0.13 - lowest  

anova(mod_0, mod_1) # 5.286e-12 ***
anova(mod_0, mod_2) # 0.0236 *
anova(mod_0, mod_3) # 4.047e-10 ***

anova(mod_1, mod_2) #

anova(mod_1, mod_3) # 
anova(mod_2, mod_3) # 2.192e-10 ***

# GAM and body weight: Select and shape data ----

# Hypotheses are
# - "Body weight of (fe)male offspring is dependent on both parents obesity"
# - "Body weight of (fe)male offspring is dependent on mothers' obesity"
# - "Body weisght of (fe)male offspring is dependent on fathers' obesity"

# _1.) Isolate data for modelling ----

mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, AnimalSex, MeasurementDay, BodyWeightG, ObeseParents) 

# _2.) Check balance of modelling data and get a graphical or table summary  ----

mice_f1_model_data %>% select(MeasurementDay, ObeseParents) %>% count(MeasurementDay, ObeseParents, sort = TRUE)
mice_f1_model_data %>% select(MeasurementDay, ObeseParents) %>% table()

plotrix::sizetree(mice_f1_model_data %>% dplyr::select(AnimalSex, ObeseParents)) 
plotrix::sizetree(mice_f1_model_data %>% dplyr::select(AnimalSex, ObeseParents, MeasurementDay)) 
plotrix::sizetree(mice_f1_model_data %>% dplyr::select(AnimalSex, ObeseParents, BodyWeightG)) 


# GAM and body weight: Model offspring' body weight as function of parents obesity  ----

# https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# https://rdrr.io/cran/mgcv/man/gam.selection.html

# _1.) Define models ----

# __a) No random effects ----
mod_1 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents, data = mice_f1_model_data, method = "ML", family = "gaussian")

# __b) Subject only as random effect ----
mod_2 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents + s(AnimalId, bs = 're'), data = mice_f1_model_data, method = "ML", family = "gaussian")

# __c) Subject at each time  as random effect ----
mod_3 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents + s(AnimalId, MeasurementDay, bs = 're'), data = mice_f1_model_data, method = "ML", family = "gaussian")

# __d) Adding time correlation to Subjscts random effects instead of complicated ranodm effect ----
mod_4 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents + s(AnimalId, bs = 're'), correlation = coAR1(form = ~ MeasurementDay | AnimalId), data = mice_f1_model_data, method = "ML",  family = "gaussian")

# _2.) Model summaries ----

summary(mod_1)
summary(mod_2)
summary(mod_3)
summary(mod_4)

# _3.) Model appraisals ----

gratia::appraise(mod_1)                            # no more heteroscedasticity 
gratia::appraise(mod_2)
gratia::appraise(mod_3) # fits well- some inverse correlation in responses vs fitted values - perhaps negligible
gratia::appraise(mod_4) # also quite good


# _4.) Inspect smoother and confidence intervals ----

# show confidence intervals graphically (for now only)
# - see https://stats.stackexchange.com/questions/33327/confidence-interval-for-gam-model
# - see dahed lines in plot below

plot.gam(mod_3, residuals = TRUE, rug = TRUE, pages = 1, all.terms = TRUE)


# _5.) Inspect model predictions ----

mod_3_predictions  <- get_model_prdictions_with_mgcv(mod_3, mice_f1_model_data)

ggplot(data = mod_3_predictions, aes(x = MeasurementDay, y = BodyWeightG, colour = AnimalId)) +
  geom_point(aes(y = BodyWeightG, group = AnimalId), alpha = 0.5) +
  geom_line(aes(y = fit, group = AnimalId), alpha = 0.5, linewidth = 0.2) +
  facet_wrap(ObeseParents ~ AnimalSex, ncol = 2) + 
  theme_bw() + 
  labs(title = "Offsprings body weight by sex and parents obesity statuts", 
       subtitle = paste("R model formula: ", as.character(paste(deparse(formula(mod_3), width.cutoff = 500), collapse=""))),
       x="age [d]", y = "body weight [g]")
  

# Save finished data ----
saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_H3variables.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_H3variables.rds"))

# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "020_r_h1.RData"))
renv::snapshot()
