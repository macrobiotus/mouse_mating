#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis 3: Offspring-specific sex effect of inheriting obesity.Superseded by modelling with {saemix}."
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

#' # Prepare environment

# Prepare environment  ---- 

# _1.) Collect garbage ----

rm(list=ls())
gc()
stop("This code is outdated and likely not needed anymore. Before running name diets in code to match data - see chnaglog file on 27.05.2024")

# _2.) Packages ----

library("here")     # environment management - use package "here" in conjunction with a RStudio project
library("renv")     # environment can be snap shot

library("readxl")   # table I/O
library("janitor")  # data cleaning 
library("hablar")
library("tidyr")

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

# _3.) Functions ----

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

#' # Get data

# Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_H2variables.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_H2variables.rds"))

#' # Logistic Regression: Select and shape data

# Logistic Regression: Select and shape data ----

# Hypotheses are
# - "Obesity of (fe)male offspring is dependent on both parents obesity"
# - "Obesity of (fe)male offspring is dependent on mothers' obesity"
# - "Obesity of (fe)male offspring is dependent on fathers' obesity"

# _1.) Isolate data for modelling ----

mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, AnimalSex, ObesityLgcl, ObeseParents) 

# _2.) Check balance of modelling data and get a graphical or table summary  ----

mice_f1_model_data %>% select(ObesityLgcl, ObeseParents) %>% dplyr::count(ObesityLgcl, ObeseParents, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, ObeseParents) %>% table()

#' # Logistic Regression: Model offspring' obesity as function of parents obesity

# Logistic Regression: Model offspring' obesity as function of parents obesity  ----

# _1.) Intercept-only model ----

if (!length(unique(mice_f1_model_data[["ObesityLgcl"]])) == 1){
  
  mod_0 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
  equatiomatic::extract_eq(mod_0)
  summary(mod_0)

}

# _2.) Animal sex only ----

if (!length(unique(mice_f1_model_data[["ObesityLgcl"]])) == 1){
 
   mod_1 <- lme4::glmer(ObesityLgcl ~ AnimalSex +  (1 | AnimalId), data = mice_f1_model_data, family = binomial)
  equatiomatic::extract_eq(mod_1)
  summary(mod_1)  # AnimalSexm *** (due to weight curve Obesity definistion)
  
}
  
# _3.) Parents obesity only  ----

if (!length(unique(mice_f1_model_data[["ObesityLgcl"]])) == 1){
  
  mod_2 <- lme4::glmer(ObesityLgcl ~ ObeseParents + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
  equatiomatic::extract_eq(mod_2)
  summary(mod_2) # ObeseParentsMotherFatherObese ***

}

# _4.) Full model ----

if (!length(unique(mice_f1_model_data[["ObesityLgcl"]])) == 1){

mod_3 <- lme4::glmer(ObesityLgcl ~ AnimalSex + ObeseParents + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_3)
summary(mod_3)

exp(fixef(mod_3))
round(exp(fixef(mod_3)), digits = 3)

}

# _5.) Test logistic regressions' fixed effects using ANOVA function and cAIC ----

if (!length(unique(mice_f1_model_data[["ObesityLgcl"]])) == 1){

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

} 

#' # GAM and body weight: Select and shape data

# GAM and body weight: Select and shape data ----

# Hypotheses are
# - "Body weight of (fe)male offspring is dependent on both parents obesity"
# - "Body weight of (fe)male offspring is dependent on mothers' obesity"
# - "Body weisght of (fe)male offspring is dependent on fathers' obesity"

# _1.) Isolate data for modelling ----

mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, AnimalSex, MeasurementDay, BodyWeightG, ObesityLgcl, ObeseParents, MotherDiet, FatherDiet) 
mice_f1_model_data$AnimalSex <- relevel(mice_f1_model_data$AnimalSex, "f")
mice_f1_model_data$ObeseParents <- relevel(mice_f1_model_data$ObeseParents, "MotherFatherNotObese")

# _2.) Check balance of modelling data and get a graphical or table summary  ----

mice_f1_model_data %>% select(MeasurementDay, ObeseParents) %>% dplyr::count(MeasurementDay, ObeseParents, sort = TRUE)
mice_f1_model_data %>% select(MeasurementDay, ObeseParents) %>% table()

plotrix::sizetree(mice_f1_model_data %>% dplyr::select(AnimalSex, ObeseParents)) 
plotrix::sizetree(mice_f1_model_data %>% dplyr::select(AnimalSex, ObeseParents, MeasurementDay)) 
plotrix::sizetree(mice_f1_model_data %>% dplyr::select(AnimalSex, ObeseParents, BodyWeightG)) 

# _3.) Data summary ----

mice_f1_model_data
summary(mice_f1_model_data)

#' # GAM and body weight: Model offspring' body weight as function of parents obesity

# GAM and body weight: Model offspring' body weight as function of parents obesity  ----

# https://fromthebottomoftheheap.net/2021/02/02/random-effects-in-gams/
# https://rdrr.io/cran/mgcv/man/gam.selection.html

# _1.) Define models ----

# __a) No random effects ----
mod_1 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents, data = mice_f1_model_data, method = "ML", family = "gaussian")
equatiomatic::extract_eq(mod_1)

# __b) No random effects, but correlated errors ----
mod_2 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents, data = mice_f1_model_data, method = "ML", correlation = coAR1(form = ~ MeasurementDay | AnimalId), family = "gaussian")
equatiomatic::extract_eq(mod_2)

# __c) Subject only as random effect ----
mod_3 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents + s(AnimalId, bs = 're'), data = mice_f1_model_data, method = "ML", family = "gaussian")
equatiomatic::extract_eq(mod_3)

# __d) Subject at each time  as random effect ----
mod_4 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents + s(AnimalId, MeasurementDay, bs = 're'), data = mice_f1_model_data, method = "ML", family = "gaussian")
equatiomatic::extract_eq(mod_4)

# __e) Adding time correlation to Subjscts random effects instead of complicated ranodm effect ----
mod_5 <- gam(BodyWeightG ~  s(MeasurementDay, k=5, bs="tp") + AnimalSex + ObeseParents + s(AnimalId, bs = 're'), correlation = coAR1(form = ~ MeasurementDay | AnimalId), data = mice_f1_model_data, method = "ML",  family = "gaussian")
equatiomatic::extract_eq(mod_5)

# __f) Reverting to d) but for each sex ----
mod_6 <- gam(BodyWeightG ~  s(MeasurementDay, by = AnimalSex, k=5, bs="fs", m=2)  + ObeseParents + s(AnimalId, MeasurementDay, bs = 're'), data = mice_f1_model_data, method = "ML", family = "gaussian")
equatiomatic::extract_eq(mod_6)

# _2.) Model rankings ----

AIC(mod_1)
AIC(mod_2)
AIC(mod_3)
AIC(mod_4) # lowest
AIC(mod_5)
AIC(mod_6) # even lower

# _3.) Model summaries ----

summary(mod_1)
summary(mod_2)
summary(mod_3)
summary(mod_4) # Deviance explained = 96.6%
summary(mod_5)
summary(mod_6) # Deviance explained = 96.9%

# _4.) Model appraisals ----

gratia::appraise(mod_1)                            # no more heteroscedasticity 
gratia::appraise(mod_2)
gratia::appraise(mod_3) 
gratia::appraise(mod_4) # fits well- some inverse correlation in responses vs fitted values - perhaps negligible
gratia::appraise(mod_5) 
gratia::appraise(mod_6) 

# _5.) Inspect smoother and confidence intervals ----

# show confidence intervals graphically (for now only)
# - see https://stats.stackexchange.com/questions/33327/confidence-interval-for-gam-model
# - see dahed lines in plot below

plot.gam(mod_6, residuals = TRUE, rug = TRUE, pages = 1, all.terms = TRUE)
confint(mod_4, parm = NULL, level = 0.95)
anova(mod_3, mod_4)

# _6.) Inspect model predictions ----

# unconditional = TRUE
mod_4_predictions  <- get_model_prdictions_with_mgcv(mod_4, mice_f1_model_data)

ggplot(data = mod_4_predictions, aes(x = MeasurementDay, y = BodyWeightG, colour = AnimalId)) +
  geom_point(aes(y = BodyWeightG, group = AnimalId), alpha = 0.5) +
  geom_line(aes(y = fit, group = AnimalId), alpha = 0.5, linewidth = 0.2) +
  facet_wrap(ObeseParents ~ AnimalSex, ncol = 2) + 
  theme_bw() + 
  labs(title = "Offsprings body weight by sex and parents obesity statuts", 
       subtitle = paste("R model formula: ", as.character(paste(deparse(formula(mod_4), width.cutoff = 500), collapse=""))),
       x="age [d]", y = "body weight [g]")

mod_6_predictions  <- get_model_prdictions_with_mgcv(mod_6, mice_f1_model_data)

gam6_plot <- ggplot(data = mod_6_predictions, aes(x = MeasurementDay, y = BodyWeightG, colour = AnimalId)) +
  geom_point(aes(y = BodyWeightG, group = AnimalId), alpha = 0.5) +
  geom_line(aes(y = fit, group = AnimalId), alpha = 0.5, linewidth = 0.2) +
  facet_wrap(ObeseParents ~ AnimalSex, ncol = 2) + 
  theme_bw() + 
  labs(title = "Offsprings body weight by sex and parents obesity statuts", 
       subtitle = paste("R model formula: ", as.character(paste(deparse(formula(mod_6), width.cutoff = 500), collapse=""))),
       x="age [d]", y = "body weight [g]")

gam6_plot

ggsave(device = cairo_pdf, plot = gam6_plot, width = 210, height = 210, units = c("mm"), dpi = 300,scale = 1.2, path = here("../manuscript/display_items"), filename = "040_r_h3__gam_6_F1_weight_trajectories.pdf")

# _7.) Check which data can be used for RNAseq ----

# __a) Read in RNA seq metadata, check, and format ----

mice_f1_rna_seq <- readxl::read_excel("/Users/paul/Documents/HM_MouseMating/communication/190916 Probenliste Clariom S.xlsx")
mice_f1_rna_seq %>% print(n = Inf) # inspect
mice_f1_rna_seq %<>% clean_names(case = "upper_camel")
mice_f1_rna_seq %<>% select(-c("X6"))
mice_f1_rna_seq %<>% rename(DietGroup = X7)
mice_f1_rna_seq %<>% tidyr::fill("Sex","ParentalDietMoFa", "DietGroup")
mice_f1_rna_seq %<>% convert(fct(Animal, Tissue, Sex, ParentalDietMoFa, DietGroup))
mice_f1_rna_seq_no_tissues <- mice_f1_rna_seq %>% select(-c(Sample, Tissue)) %>% distinct()

# __b) Format data from modelling frame  ----

mice_f1_model_data_uniques <- mice_f1_model_data %>% select(-c(BodyWeightG, MeasurementDay)) %>% distinct()

# __c) Join both data sets to see what can be done with RNAseq data ---

mice_f1_model_data_rna_seqed <- left_join(mice_f1_model_data_uniques, mice_f1_rna_seq_no_tissues, by = c("AnimalId" = "Animal")) 
mice_f1_modeled_data_with_rna_seq_data <- mice_f1_model_data_rna_seqed %>% filter(Sex != "NA" & ParentalDietMoFa != "NA" & DietGroup != "NA")

# __d) Check data prior to export

# This is SI Table 2 (30-Aug-2023):
mice_f1_modeled_data_with_rna_seq_data %>% select(ObesityLgcl, ObeseParents) %>% table()
mice_f1_modeled_data_with_rna_seq_data

# __f) Add parental information

mice_f1_modeled_data_with_rna_seq_data <- left_join(mice_f1_modeled_data_with_rna_seq_data,  distinct(mice_f1_slct[c("AnimalId", "MotherId", "FatherId")]))

# _8.) Export data as Excel file for RNAseq analysis ----

# __a) In below object explore (1) differential expression, (2) GO and KEGG terms, 
#     between individuals for (a) female or (b) male  animal sex or (c) both sexes in unison
#     always between the factor levles "MotherFatherObese" and "MotherFatherNotObese" of factor "ObeseParents".
#     for each obese or not obese offspring. See notes 12.05.2023 in Communications folder and README.md 

saveRDS(mice_f1_model_data, file = here("rds_storage", "040_r_h3__mice_f1_model_data.rds"))
saveRDS(mice_f1_modeled_data_with_rna_seq_data, file = here("rds_storage", "040_r_h3__mice_f1_modeled_data_with_rna_seq_data.rds"))

# __b) Generate a useful SI Table 1
mf1_summary_xlsx <- mice_f1_modeled_data_with_rna_seq_data
mf1_summary_xlsx %<>% select(AnimalId, AnimalSex, ObesityLgcl, ObeseParents, MotherDiet, FatherDiet, MotherId, FatherId)
openxlsx::write.xlsx(mf1_summary_xlsx, paste0(here("tables"), "/", "040_r_h3__rna_seq_sample.xlsx"), asTable = TRUE, overwrite = TRUE)

# __b) Generate a useful summary of SI Table 1 for the main text

mf1_summary_xlsx %>% select(MotherDiet, AnimalSex) %>% table()


#' # Save finished data

# Save finished data ----
saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_H3variables.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_H3variables.rds"))

#' # Snapshot environment

# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "040_r_h3.RData"))
renv::snapshot()
