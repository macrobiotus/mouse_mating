#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis 2: parents sex-specific obesity inheritance – effect of mothers or fathers’ obesity"
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

# Prepare environment ---- 

#' ## Collect garbage

# _1.) Collect garbage ----

rm(list=ls())
gc()

#' # Packages

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

#' # Get data

# Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_H1variables.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_H1variables.rds"))

#' # Select and shape data for Subquestion 1: "Obesity in fathers"

# Select and shape data for Subquestion 1: "Obesity in fathers" ----

#' ## Get a factor that defines "parents’ obesity for each father"

# _1.) Get a factor that defines "parents’ obesity for each father" ----

# - check values used for logical factor definition 

levels(mice_f1_slct$ObeseParents)
mice_f1_slct %>% select(ObeseParents) %>% dplyr::count(ObeseParents, sort = TRUE)

# - define logical factors for parents

mice_f1_slct %<>% mutate(ObeseFatherLgcl = case_when(
                        (ObeseParents == "FatherObese") ~ TRUE,
                        (ObeseParents == "MotherFatherNotObese") ~ FALSE,
                        (ObeseParents == "MotherObese") ~ FALSE, 
                        (ObeseParents == "MotherFatherObese") ~ TRUE))

#' ## Get a factor that defines "obesity of any-sex-offspring"

# _2.) Get a factor that defines "obesity of any-sex-offspring" ----

# - the one already there is cumbersome

levels(mice_f1_slct$Obesity)

# - already defined one in the previous script                     

#' ## Isolate data for modelling

# _3.) Isolate data for modelling ----
  
mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, ObesityLgcl, ObeseFatherLgcl) 

#' ## Check balance of modelling data and get a graphical or table summary

# _4.) Check balance of modelling data and get a graphical or table summary  ----

mice_f1_model_data %>% select(ObesityLgcl, ObeseFatherLgcl) %>% dplyr::count(ObesityLgcl, ObeseFatherLgcl, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, ObeseFatherLgcl) %>% table()

#' # Logistic Regression: Model offsprings' obesity as function of Father's obesity

# Logistic Regression: Model offsprings' obesity as function of father's obesity  ----

#' ## Intercept-only model

# _1.) Intercept-only model ----

mod_0 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_0)
summary(mod_0)

#' ## Actual model

# _2.) Actual model ----

mod_1 <- lme4::glmer(ObesityLgcl ~ ObeseFatherLgcl + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_1)
summary(mod_1)

#' ## Test effect of parents obesity status

# _3.) Test effect of parents obesity status ----

anova(mod_0, mod_1)

#' # Select and shape data for Subquestion 2: "Obesity in mothers"

# Select and shape data for Subquestion 2: "Obesity in mothers" ----

#' ## Get a factor that defines "parents’ obesity for each mother"

# _1.) Get a factor that defines "parents’ obesity for each mother" ----

# - check values used for logical factor definition 

levels(mice_f1_slct$ObeseParents)
mice_f1_slct %>% select(ObeseParents) %>% dplyr::count(ObeseParents, sort = TRUE)

# - define logical factors for parents

mice_f1_slct %<>% mutate(ObeseMotherLgcl = case_when(
  (ObeseParents == "FatherObese") ~ FALSE,
  (ObeseParents == "MotherFatherNotObese") ~ FALSE,
  (ObeseParents == "MotherObese") ~ TRUE, 
  (ObeseParents == "MotherFatherObese") ~ TRUE))

#' ## Get a factor that defines "obesity of any-sex-offspring"

# _2.) Get a factor that defines "obesity of any-sex-offspring" ----

# - the one already there is cumbersome

levels(mice_f1_slct$Obesity)

# - already defined one in the previous script                     

#' ## Isolate data for modelling

# _3.) Isolate data for modelling ----

mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, ObesityLgcl, ObeseMotherLgcl) 

#' ## Check balance of modelling data and get a graphical or table summary

# _4.) Check balance of modelling data and get a graphical or table summary ----

mice_f1_model_data %>% select(ObesityLgcl, ObeseMotherLgcl) %>% dplyr::count(ObesityLgcl, ObeseMotherLgcl, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, ObeseMotherLgcl) %>% table()

#' # Logistic Regression: Model offsprings' obesity as function of mothers obesity

# Logistic Regression: Model offsprings' obesity as function of mothers obesity  ----

#' ## Intercept-only model

# _1.) Intercept-only model ----

mod_3 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_3)
summary(mod_0)

#' ## Actual model

# _2.) Actual model ----

mod_4 <- lme4::glmer(ObesityLgcl ~ ObeseMotherLgcl + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_4)
summary(mod_4)

#' ## Test effect of parents obesity status

# _3.) Test effect of parents obesity status ----

anova(mod_3, mod_4)

#' # Save finished data

# Save finished data ----
saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_H2variables.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_H2variables.rds"))

#' # Snapshot environment

# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "030_r_h2.RData"))
renv::snapshot()

