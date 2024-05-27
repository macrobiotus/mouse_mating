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

#' # Prepare environment

# Prepare environment  ---- 

rm(list=ls())
gc()

#' ## Load Packages

# _1.) Load Packages ----

library("here")     # environment management - use package "here" in conjunction with a RStudio project
library("renv")     # environment can be snap shot

library("readxl")   # table I/O

library("janitor")  # data cleaning 
library("dplyr")    # pipes and more
library("magrittr") # even more pipes

library("hablar") # Type conversions

#' ## Functions

# _2.) Functions ----

'%!in%' <- function(x,y)!('%in%'(x,y))


#' # Data read-in, cleaning, formatting

# Data read-in, cleaning, formatting ----

#' ## Get data

# _1.) Get data ----

xlsx_path <- here("raw_data", "221214_mice_pairing_reshaped.xlsx")
mice_f0 <- read_excel(xlsx_path, sheet = "Parental (F0)",  na = c("","-",  "NA")) %>% clean_names()
mice_f1 <- read_excel(xlsx_path, sheet = "Offspring (F1)", na = c("","-", "NO", "NA")) %>% clean_names()

#' ## Type correction

# _2.) Type correction  ----

#' ### Check types 

# __a) Check types ----

glimpse(mice_f0)
glimpse(mice_f1)

#' ### Correct types

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

#' ### Check types

# __c) Check types ----

glimpse(mice_f0)
glimpse(mice_f1)

#' ## Check data for completeness

# _3.) Check data for completeness ----

# Why are there no males in the f0?
fathers_animal_id <- c("8344", "8005", "8005", "8345", "8345", "8335", "8335", "8346", "8337")
mice_f0$animal_id %in% fathers_animal_id

# -> only coded as Partner

#' #  Data re-coding

# Data re-coding ----

#' ## Recode variables to CamelCase

# _1.) Recode variables to CamelCase ---- 

mice_f0 %<>% clean_names(case = "upper_camel")
mice_f1 %<>% clean_names(case = "upper_camel")

#' ## Add generational information to data

# _2.) Add generational information to data ---- 

mice_f0 %<>% mutate("Generation" = as.factor("f0"))
mice_f1 %<>% mutate("Generation" = as.factor("f1"))

#' ## Re-code measurement timing - factor week to numerical days

# _3.) Re-code measurement timing - factor week to numerical days ----

mice_f0 %<>% mutate("MeasurementDay" = readr::parse_number(as.character(Week)) * 7)
mice_f1 %<>% mutate("MeasurementDay" = readr::parse_number(as.character(Week)) * 7)

glimpse(mice_f0)
glimpse(mice_f1)

# _4.) Rencoding dietary variables for second submission (first revision) ----

# Inserted 27.05.2024 - r

warning("Re-enecoded dietary variables here on 26.05.2024")

mice_f0 %<>% mutate(Diet = recode(Diet, "HFD" = "WD", "CD" = "CD"))
mice_f0 %<>% mutate(PartnerDiet = recode(PartnerDiet, "HFD" = "WD", "CD" = "CD"))
mice_f0 %<>% mutate(Group = recode(Group, "HFD" = "WD", "CD" = "CD"))

mice_f1 %<>% mutate(MotherDiet = recode(MotherDiet, "HFD" = "WD", "CD" = "CD"))
mice_f1 %<>% mutate(FatherDiet = recode(FatherDiet, "HFD" = "WD", "CD" = "CD"))

#' ## Check which animals have been used for clinical assessments and save original state

# _5.) Check which animals have been used for clinical assessments and save original state ----

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

#' ## Deal with "Mixed" parental diets by removing them

# _6.) Remove "Mixed" parental diets and save altered state  ----

#' ### Show f0 animals that received mixed diets

# __a) Show f0 animals that received mixed diets ----

mice_f0 %>% filter(AnimalId %in% c(8994, 8995, 8996, 8997)) %>% arrange(AnimalId) %>% dplyr::select(AnimalId, AnimalSex,  Diet, Group) %>% print(n = Inf)

#' ### Show f0 animals that did not receive mixed diets

# __b) Show f0 animals that did not receive mixed diets ----

mice_f0 %>% filter(AnimalId %!in% c(8994, 8995, 8996, 8997)) %>% arrange(AnimalId) %>% dplyr::select(AnimalId, AnimalSex, Diet, Group) %>% print(n = Inf)

#' ### Keep only f0 animals that did not receive mixed diets

# __c) Keep only f0 animals that did not receive mixed diets ----

mice_f0 %<>% filter(AnimalId %!in% c(8994, 8995, 8996, 8997))

#' ### Remove f1 individuals from parents with mixed diets, if any

# __d) Remove f1 individuals from mothers with "mixed" diets, if any  ----

# -- testing approach --

# there are 50 f1 mice prior to filtering
mice_f1 %>% pull(AnimalId) %>% unique() %>% length 

# there is data for 8 female mice in the f0 
mice_f0 %>% dplyr::select(AnimalId, AnimalSex) %>% distinct
mice_f0 %>% dplyr::select(AnimalId, AnimalSex) %>% pull(AnimalId) %>% unique() %>% length 

# there are 50 f1 mice after filtering
mice_f1 %>% filter(MotherId %in% (mice_f0 %>% pull(AnimalId) %>% unique())) %>% pull(AnimalId) %>% unique() %>% length 

# -- doing the filtering --
mice_f1 %<>% filter(MotherId %in% (mice_f0 %>% pull(AnimalId) %>% unique()))

#' ### Remove Group column among f0 and mixed group among f1

# __e) Remove Group column among f0 and mixed group among f1  ----

# to not get confused downstream and beacuse ther is no RNAseq data for mixed f1 
mice_f0 %<>% dplyr::select(-c("Group")) 
mice_f1 %<>% filter(MotherDiet != "Mix")

#' ### Save altered state

# __f) Save altered state ----

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

#' ## Check which animals have been used for RNA sequencing

# _7.) Check which animals have been used for RNA sequencing ----

mice_f1_rna_seq <- readxl::read_excel("/Users/paul/Documents/HM_MouseMating/communication/190916 Probenliste Clariom S.xlsx")
mice_f1_rna_seq %>% print(n = Inf) # inspect
mice_f1_rna_seq %<>% clean_names(case = "upper_camel")
mice_f1_rna_seq %<>% dplyr::select(-c("X6"))
mice_f1_rna_seq %<>% rename(DietGroup = X7)
mice_f1_rna_seq %<>% tidyr::fill("Sex","ParentalDietMoFa", "DietGroup")
mice_f1_rna_seq %<>% convert(fct(Animal, Tissue, Sex, ParentalDietMoFa, DietGroup))
mice_f1_rna_seq_no_tissues <- mice_f1_rna_seq %>% dplyr::select(-c(Sample, Tissue)) %>% distinct()

#' ## Output a table that show which animals had been used for RNA sequencing

# _8.) Tabularize which animals have been used for RNA sequencing  ----

mice_f1_modeled_data_with_rna_seq_data <- left_join(mice_f1_animal_ids, mice_f1_rna_seq_no_tissues, by = c("AnimalId" = "Animal"))
mice_f1_modeled_data_with_rna_seq_data %<>% dplyr::select(-c("Sex", "ParentalDietMoFa"))
mice_f1_modeled_data_with_rna_seq_data %<>% mutate(RNAseq = case_when(is.na(DietGroup)  ~ FALSE, !is.na(DietGroup) ~ TRUE))
mice_f1_modeled_data_with_rna_seq_data %<>% dplyr::select(-c("DietGroup"))

openxlsx::write.xlsx(mice_f1_modeled_data_with_rna_seq_data, paste0(here("tables"), "/", "000_r_format_data__rna_seq_sample.xlsx"), asTable = TRUE, overwrite = TRUE)

mice_f1_modeled_data_with_rna_seq_data %>% print(n = Inf)

#' ## Select all variables that could be relevant for modelling

# _9.) Select all variables that could be relevant for modelling ----

mice_f0_slct <- mice_f0 # 
mice_f1_slct <- mice_f1 # %>%  dplyr::select(animal_id, animal_sex, body_weight_g, mother_diet, father_diet, week, mother_id, father_id) 

# check selected variables
glimpse(mice_f0_slct)
glimpse(mice_f1_slct)

#' # Save finished data

# Save finished data ----

saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct.rds"))

#' # Snapshot environment

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "000_r_format_data.RData"))
renv::snapshot()
