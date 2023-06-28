#' ---
#' title: "Mice Mating Study"
#' subtitle: "Raw data formatting"
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

# Prepare environment  ---- 

rm(list=ls())
gc()

# _1.) Load Packages ----

library("here")     # environment management - use package "here" in conjunction with a RStudio project
library("renv")     # environment can be snap shot

library("readxl")   # table I/O

library("janitor")  # data cleaning 
library("dplyr")    # pipes and more
library("magrittr") # even more pipes

library("hablar") # Type conversions


# Data read-in, cleaning, formatting ----

# _1.) Get data ----

xlsx_path <- here("raw_data", "221214_mice_pairing_reshaped.xlsx")
mice_f0 <- read_excel(xlsx_path, sheet = "Parental (F0)",  na = c("","-",  "NA")) %>% clean_names()
mice_f1 <- read_excel(xlsx_path, sheet = "Offspring (F1)", na = c("","-", "NO", "NA")) %>% clean_names()


# _2.) Type correction  ----

# __a) Check types ----

glimpse(mice_f0)
glimpse(mice_f1)

# __b) Correct types ----

mice_f0 %<>% mutate(across(where(is.character), factor))
mice_f0 %<>% mutate(week = factor(week))
mice_f0 %<>% mutate(animal_id = factor(animal_id))
mice_f0 %<>% mutate(mating_with = factor(mating_with))

mice_f1 %<>% mutate(across(where(is.character), factor))
mice_f1 %<>% mutate(week = factor(week))
mice_f1 %<>% mutate(mother_id = factor(mother_id))
mice_f1 %<>% mutate(father_id = factor(father_id))
mice_f1 %<>% mutate(animal_id = factor(animal_id))

# __c) Check types ----

glimpse(mice_f0)
glimpse(mice_f1)

# _3.) Check data for completeness ----

# Why are there no males in the f0?
fathers_animal_id <- c("8344", "8005", "8005", "8345", "8345", "8335", "8335", "8346", "8337")
mice_f0$animal_id %in% fathers_animal_id

# -> only coded as Partner

# _3.) Data re-coding ----

# __a) Recode variables to CamelCase ---- 

mice_f0 %<>% clean_names(case = "upper_camel")
mice_f1 %<>% clean_names(case = "upper_camel")

# __b) Add generational information to data ---- 

mice_f0 %<>% mutate("Generation" = as.factor("f0"))
mice_f1 %<>% mutate("Generation" = as.factor("f1"))

# __c) Re-code measurement timing - factor week to numerical days ----

mice_f0 %<>% mutate("MeasurementDay" = readr::parse_number(as.character(Week)) * 7)
mice_f1 %<>% mutate("MeasurementDay" = readr::parse_number(as.character(Week)) * 7)

glimpse(mice_f0)
glimpse(mice_f1)


# __d) Check which animals have been used for clinical assessments and save original state ----

# f0 generation - check 
glimpse(mice_f0)
mice_f0_animal_ids <- mice_f0 %>% dplyr::select(AnimalId, AnimalSex, Group, Diet, PartnerDiet) %>% distinct()
mice_f0_animal_ids %>% print(n = Inf)

# f0 generation - save original state
openxlsx::write.xlsx(mice_f0_animal_ids, paste0(here("tables"), "/", "000_r_format_data__f0_mice_with_all_diets.xlsx"), asTable = TRUE, overwrite = TRUE)
saveRDS(mice_f0, file = here("rds_storage", "mice_f0.rds"))

# f1 generation - check 
glimpse(mice_f1)
mice_f1_animal_ids <- mice_f1 %>% dplyr::select(AnimalId, AnimalSex, MotherDiet, FatherDiet) %>% distinct()
mice_f1_animal_ids %>% print(n = Inf)
openxlsx::write.xlsx(mice_f1_animal_ids, paste0(here("tables"), "/", "000_r_format_data__f1_mice_with_all_diets.xlsx"), asTable = TRUE, overwrite = TRUE)

# f1 generation - save original state
openxlsx::write.xlsx(mice_f1_animal_ids, paste0(here("tables"), "/", "000_r_format_data__f0_mice_with_all_diets.xlsx"), asTable = TRUE, overwrite = TRUE)
saveRDS(mice_f1, file = here("rds_storage", "mice_f1.rds"))

# __e) Remove Group column among f0 and mixed group among f1  ----

# to not get confused downstream and beacuse ther is no RNAseq data for mixed f1 
mice_f0 %<>% select(-c("Group")) 
mice_f1 %<>% filter(MotherDiet != "Mix")

# f0 generation - save altered state
openxlsx::write.xlsx(mice_f0, paste0(here("tables"), "/", "000_r_format_data__f0_mice_cleaned_detail.xlsx"), asTable = TRUE, overwrite = TRUE)
saveRDS(mice_f1, file = here("rds_storage", "mice_f0_cleaned.rds"))

openxlsx::write.xlsx(mice_f0 %>% dplyr::select(AnimalId, AnimalSex, Diet, PartnerDiet) %>% distinct(), paste0(here("tables"), "/", "000_r_format_data__f0_mice_cleaned.xlsx"), asTable = TRUE, overwrite = TRUE)
saveRDS(mice_f1, file = here("rds_storage", "mice_f0_cleaned.rds"))


# f1 generation - save altered state
openxlsx::write.xlsx(mice_f1, paste0(here("tables"), "/", "000_r_format_data__f1_mice_cleaned_detail.xlsx"), asTable = TRUE, overwrite = TRUE)
saveRDS(mice_f1, file = here("rds_storage", "mice_f1_cleaned.rds"))

openxlsx::write.xlsx(mice_f1 %>% dplyr::select(AnimalId, AnimalSex, MotherDiet, FatherDiet) %>% distinct(), paste0(here("tables"), "/", "000_r_format_data__f1_mice_cleaned.xlsx"), asTable = TRUE, overwrite = TRUE)
saveRDS(mice_f1, file = here("rds_storage", "mice_f1_cleaned.rds"))

# __d) Check which animals have been used for RNA sequencing ----

mice_f1_rna_seq <- readxl::read_excel("/Users/paul/Documents/HM_MouseMating/communication/190916 Probenliste Clariom S.xlsx")
mice_f1_rna_seq %>% print(n = Inf) # inspect
mice_f1_rna_seq %<>% clean_names(case = "upper_camel")
mice_f1_rna_seq %<>% dplyr::select(-c("X6"))
mice_f1_rna_seq %<>% rename(DietGroup = X7)
mice_f1_rna_seq %<>% tidyr::fill("Sex","ParentalDietMoFa", "DietGroup")
mice_f1_rna_seq %<>% convert(fct(Animal, Tissue, Sex, ParentalDietMoFa, DietGroup))
mice_f1_rna_seq_no_tissues <- mice_f1_rna_seq %>% select(-c(Sample, Tissue)) %>% distinct()


# __e) Output a table that show which animals had been used for RNA sequencing  ----

mice_f1_modeled_data_with_rna_seq_data <- left_join(mice_f1_animal_ids, mice_f1_rna_seq_no_tissues, by = c("AnimalId" = "Animal"))
mice_f1_modeled_data_with_rna_seq_data %<>% dplyr::select(-c("Sex", "ParentalDietMoFa"))
mice_f1_modeled_data_with_rna_seq_data %<>% mutate(RNAseq = case_when(is.na(DietGroup)  ~ FALSE, !is.na(DietGroup) ~ TRUE))
mice_f1_modeled_data_with_rna_seq_data %<>% dplyr::select(-c("DietGroup"))

openxlsx::write.xlsx(mice_f1_modeled_data_with_rna_seq_data, paste0(here("tables"), "/", "000_r_format_data__rna_seq_sample.xlsx"), asTable = TRUE, overwrite = TRUE)

# __f) Select all variables that could be relevant for modelling ----

mice_f0_slct <- mice_f0 # 
mice_f1_slct <- mice_f1 # %>%  dplyr::select(animal_id, animal_sex, body_weight_g, mother_diet, father_diet, week, mother_id, father_id) 

# check selected variables
glimpse(mice_f0_slct)
glimpse(mice_f1_slct)

# Save finished data ----

saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct.rds"))

# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "000_r_format_data.RData"))
renv::snapshot()
