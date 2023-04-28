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





#'
#' # Data read-in, cleaning, formatting
#'

# Data read-in, cleaning, formatting ----

#' ## Get data

# _1.) Get data ----

xlsx_path <- here("raw_data", "221214_mice_pairing_reshaped.xlsx")
mice_f0 <- read_excel(xlsx_path, sheet = "Parental (F0)",  na = c("","-",  "NA")) %>% clean_names()
mice_f1 <- read_excel(xlsx_path, sheet = "Offspring (F1)", na = c("","-", "NO", "NA")) %>% clean_names()

#' ## Type correction

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

# _3.) Save intermediate data ----

saveRDS(mice_f0, file = here("rds_storage", "mice_f0.rds"))
saveRDS(mice_f1, file = here("rds_storage", "mice_f1.rds"))


# _4.) Data merging ----

# [use section if f0 data needs to be joined with f1 data]

#' ## Data re-coding

# _5.) Data re-coding ----

# __a) Subset variables for further use ----

# check which variables are available for selection
glimpse(mice_f1)

# select all variables that could be relavnt for modelling
mice_f1_slct <- mice_f1 %>% 
  dplyr::select(animal_id, animal_sex, body_weight_g, mother_diet, father_diet, week, mother_id, father_id) 

# check selected variables
glimpse(mice_f1_slct)
mice_f1_slct %>% print(n = Inf)

# __b) Re-code factor week as numerical days ----
mice_f1_slct %<>% mutate(day = readr::parse_number(as.character(week)) * 7) 

