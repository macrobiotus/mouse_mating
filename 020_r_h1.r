#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis 1"
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

# Data read-in, cleaning, formatting ----

# _1.) Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# _2.) Select and shape data for question ----

# Question is "obesity of any-sex-offspring is dependent on any-sex parents’ obesity"

# __a) Get a factor that defines "any-sex parents’ obesity" ----

# - check values used for logical factor definition 

levels(mice_f1_slct$ObeseParents)
mice_f1_slct %>% select(ObeseParents) %>% count(ObeseParents, sort = TRUE)

# - define logical factors for parents

mice_f1_slct %<>% mutate(ObeseParentsLgcl = case_when(
                        (ObeseParents == "FatherObese") ~ TRUE,
                        (ObeseParents == "MotherFatherNotObese") ~ FALSE,
                        (ObeseParents == "MotherObese") ~ TRUE, 
                        (ObeseParents == "MotherFatherObese") ~ TRUE))


# __b) Get a factor that defines "obesity of any-sex-offspring" ----

# - the one already there is cumbersome

levels(mice_f1_slct$Obesity)

# - define logical factors for offspring                    

mice_f1_slct %<>% mutate(ObesityLgcl = case_when(
  (Obesity == "Obese") ~ TRUE,
  (Obesity == "NotObese") ~ FALSE))

# __c) Isolate data for modelling ----
  
mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, ObesityLgcl, ObeseParentsLgcl) 

# Check balance of modelling data and get a graphical or table summary  ----

mice_f1_model_data %>% select(ObesityLgcl, ObeseParentsLgcl) %>% count(ObesityLgcl, ObeseParentsLgcl, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, ObeseParentsLgcl) %>% table()

# Model offsprings' obesity as function of parents obesity  ----

# 1.) Most simple case: logistic regression ----

# __a) Intercept-only model ----

mod_0 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
summary(mod_0)

# __b) Actual model ----
mod_1 <- lme4::glmer(ObesityLgcl ~ ObeseParentsLgcl + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
summary(mod_1)

# __c) Test effect of parents obesity status ----

anova(mod_0, mod_1)


# Save finished data ----
saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_H1variables.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_H1variables.rds"))

# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "020_r_h1.RData"))
renv::snapshot()




