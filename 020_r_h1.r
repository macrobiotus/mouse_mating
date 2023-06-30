#' ---
#' title: "Mice Mating Study"
#' subtitle: "Modelling for Hypothesis 1: Inheritance of Obesity – Any or both parents’ obesity"
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

#' # Prepare environment

# Prepare environment ---- 

#' ## Collect garbage

# _1.) Collect garbage ----

rm(list=ls())
gc()

#' ## Packages

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

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

#' # Select and shape data for Subquestion 1: "Any-sex or both parents’ obesity"

# Select and shape data for Subquestion 1: "Any-sex or both parents’ obesity" ----

# Hypothesis is "Obesity of any-sex-offspring is dependent on any-sex parents’ obesity"

#' ## Get a factor that defines "any-sex parents’ obesity"

# _1.) Get a factor that defines "any-sex parents’ obesity" ----

# - check values used for logical factor definition 

levels(mice_f1_slct$ObeseParents)
mice_f1_slct %>% select(AnimalId, ObeseParents) %>% distinct() %>% dplyr::count(ObeseParents, sort = TRUE)

# - define logical factors for parents

mice_f1_slct %<>% mutate(ObeseParentsLgcl = case_when(
                        (ObeseParents == "FatherObese") ~ TRUE,
                        (ObeseParents == "MotherFatherNotObese") ~ FALSE,
                        (ObeseParents == "MotherObese") ~ TRUE, 
                        (ObeseParents == "MotherFatherObese") ~ TRUE))

# - check result of previous command

mice_f1_slct %>% select(AnimalId, ObeseParentsLgcl) %>% distinct() %>% dplyr::count(ObeseParentsLgcl, sort = TRUE)


#' ## Get a factor that defines "obesity of any-sex-offspring"

# _2.) Get a factor that defines "obesity of any-sex-offspring" ----

# - the one already there is cumbersome

levels(mice_f1_slct$Obesity)

# - define logical factors for offspring                    

mice_f1_slct %<>% mutate(ObesityLgcl = case_when(
  (Obesity == "Obese") ~ TRUE,
  (Obesity == "NotObese") ~ FALSE))

# - check result of previous command

mice_f1_slct %>% select(AnimalId, ObesityLgcl) %>% distinct() %>% dplyr::count(ObesityLgcl, sort = TRUE)


#' ## Isolate data for modelling

# _3.) Isolate data for modelling ----
  
mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, ObesityLgcl, ObeseParentsLgcl) 

#' ## Check balance of modelling data and get a graphical or table summary

# _4.) Check balance of modelling data and get a graphical or table summary ----

mice_f1_model_data %>% select(ObesityLgcl, ObeseParentsLgcl) %>% dplyr::count(ObesityLgcl, ObeseParentsLgcl, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, ObeseParentsLgcl) %>% table()

#' # Logistic regression: Model offsprings' obesity as function of parents obesity

# Logistic regression: Model offsprings' obesity as function of parents obesity ----

#' ## Intercept-only model

# _1.) Intercept-only model ----

mod_0 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_0)
summary(mod_0)




#' ## Actual model

# _2.) Actual model ----

mod_1 <- lme4::glmer(ObesityLgcl ~ ObeseParentsLgcl + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_1)
summary(mod_1)

#' ## Test effect of parents obesity status

# _3.) Test effect of parents obesity status ----

anova(mod_0, mod_1)

#' # Select and shape data for Subquestion 2: "Both parents’ obesity only"

# Select and shape data for Subquestion 2: "Both parents’ obesity only" ----

# Hypothesis is "obesity of any-sex-offspring is dependent on both parents’ obesity"

#' ## Get a factor that defines "both parents’ obesity"

# _1.) Get a factor that defines "both parents’ obesity" ----

# - check values used for logical factor definition 

levels(mice_f1_slct$ObeseParents)
mice_f1_slct %>% select(ObeseParents) %>% dplyr::count(ObeseParents, sort = TRUE)

# - define logical factors for parents

mice_f1_slct %<>% mutate(BothObeseParentsLgcl = case_when(
  (ObeseParents == "FatherObese") ~ FALSE,
  (ObeseParents == "MotherFatherNotObese") ~ FALSE,
  (ObeseParents == "MotherObese") ~ FALSE, 
  (ObeseParents == "MotherFatherObese") ~ TRUE))

#' ## Isolate data for modelling

# _2.) Isolate data for modelling ----

mice_f1_model_data <- mice_f1_slct %>% select(AnimalId, ObesityLgcl, BothObeseParentsLgcl) 

#' ## Check balance of modelling data and get a graphical or table summary

# _3.) Check balance of modelling data and get a graphical or table summary ----

mice_f1_model_data %>% select(ObesityLgcl, BothObeseParentsLgcl) %>% dplyr::count(ObesityLgcl, BothObeseParentsLgcl, sort = TRUE)
mice_f1_model_data %>% select(ObesityLgcl, BothObeseParentsLgcl) %>% table()

#' # Logistic Regression: Model offsprings' obesity as function of parents obesity

# Logistic Regression: Model offsprings' obesity as function of parents obesity ----

#' ## Intercept-only model

# _1.) Intercept-only model ----

mod_2 <- lme4::glmer(ObesityLgcl ~ 1 + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_2)
summary(mod_2)
# modelsummary(mod_2)

#' ## Actual model

# _2.) Actual model ----

mod_3 <- lme4::glmer(ObesityLgcl ~ BothObeseParentsLgcl + (1 | AnimalId), data = mice_f1_model_data, family = binomial)
equatiomatic::extract_eq(mod_3)
summary(mod_3)

round(exp(fixef(mod_3)), digits = 2)

#' ## Test effect of parents obesity status

# _3.) Test effect of parents obesity status ----

anova(mod_2, mod_3)

#' # Save finished data

# Save finished data ----

saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_H1variables.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_H1variables.rds"))

#' # Snapshot environment

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "020_r_h1.RData"))
renv::snapshot()

