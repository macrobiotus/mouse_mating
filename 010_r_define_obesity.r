#' ---
#' title: "Mice Mating Study"
#' subtitle: "Defining obesity in mice"
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

library("lme4")     # Linear Mixed Effects Models
library("lmerTest") # Tests in Linear Mixed Effects Models
library("mgcv")     # General Additive Model

library("effects")     # Model inspection
library("performance") # Model inspection
library("cAIC4")       # Model selection 

# Data read-in, cleaning, formatting ----

# _1.) Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1.rds"))


# OLD CODE BELOW ---


# -> Find F0 mice with weight loss and add to F1 data as factor.

# _1.) Add analysis-specific variables ----

# __a) Distinguish mothers with and without weight loss prior to pregnancy ----

# find mothers that have lost weight: 
glimpse(mice_f0) 

# - finding potential mothers - done.
mice_f0_females <- mice_f0 %>% filter(animal_sex == "f")

# - defining start of pregnancy for each mouse - in weeks - relative to birth date - done.
mice_f0_females <- mice_f0_females %>% mutate(mating_week = difftime(mating_date[2], animal_birth[1], units = "weeks")) %>%
  mutate(mating_week = round(as.numeric(mating_week), digits = 0)) %>% 
  relocate(mating_week, .before = mating_date)
# - check if all mothers are the same age at mating - yes.
mice_f0_females %>% pull(mating_week) %>% unique() 


# - plotting individual weights and overall weights of mice up until mating - some may have lost weight - but not many.  
lattice::xyplot(body_weight_g ~ week, groups = animal_id, data = mice_f0_females, panel =   
                  function(x, y, ..., subscripts, groups) {
                    for (lev in levels(groups)) {
                      ok <- groups == lev
                      panel.xyplot(x[ok], y[ok], type = "smooth", col = lev)
                    }
                    panel.xyplot(x, y, type = "smooth", col = "red", lwd = 3)
                  })

# - how to identify mice who have lost more weight then others? Compare smooths?
#   see here ? https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/

#  attempting to scale weight 
center_this <- function(x){(x - mean(x, na.rm=TRUE))}
center_and_scale_that <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}


# - again for better overview - plotting scaled individual weights and overall weights of mice up until mating - on scaled data

mice_f0_females %<>% mutate(body_weight_g_scaled = center_this(body_weight_g))
mice_f0_females %<>% mutate(body_weight_g_scaled = center_and_scale_that(body_weight_g))


mice_f0_females %<>% group_by(animal_id) %>% mutate(body_weight_g_scaled = center_this(body_weight_g))
mice_f0_females %<>% group_by(animal_id) %>% mutate(body_weight_g_scaled = center_this(body_weight_g))

lattice::xyplot(body_weight_g_scaled ~ week, groups = animal_id, data = mice_f0_females, panel =   
                  function(x, y, ..., subscripts, groups) {
                    for (lev in levels(groups)) {
                      ok <- groups == lev
                      panel.xyplot(x[ok], y[ok], type = "smooth", col = lev)
                    }
                    panel.xyplot(x, y, type = "smooth", col = "red", lwd = 3)
                  })

# -> sort mothers mice grow curves by second derivatives of polynomials - half the data that is more curvy is the weight loss group

# https://stackoverflow.com/questions/48692150/calculate-derivatives-curvature-of-a-polynomial
# https://stackoverflow.com/questions/40438195/function-for-polynomials-of-arbitrary-order-symbolic-method-preferred/40442584#40442584

# - get list og polynomial models
#   using raw polynomials so thet I can use the function for getting second derivatives of data

body_weight_models = mice_f0_females %>% group_by(animal_id) %>% do(model = lm(body_weight_g_scaled ~ as.numeric(week)+I(as.numeric(week)^2), data = .))

#   check models
lapply(body_weight_models$model, summary) # not all data are curves, but that is ok
# lapply(body_weight_models$model, plot)    # far from perfect, but will do for data selection

#   isolate coefficients from model list
pc_lst <- lapply(body_weight_models[[2]], coef) 
week_lst <- split(as.numeric(mice_f0_females$week), mice_f0_females$animal_id)

#   get curevature of polynomial at each time pont
result <- Map(polynom_curvature, week_lst, pc_lst)

#   sum curvature of polynomial at each time point - smallest values have most curved grow curve
result_unlist <- lapply(result, sum) %>% unlist 

#    check if these orders make sense - yes they do. These numerical values ...
sort(result_unlist) 

# ... approximately correspond to order of plot labels 
lattice::xyplot(body_weight_g_scaled ~ week, groups = animal_id, data = mice_f0_females, panel =   
                  function(x, y, ..., subscripts, groups) {
                    for (lev in levels(groups)) {
                      ok <- groups == lev
                      panel.xyplot(x[ok], y[ok], type = "smooth", col = lev)
                      panel.text(x[ok], y[ok], labels = lev, pos=3)
                    }
                    panel.xyplot(x, y, type = "smooth", col = "red", lwd = 3)
                  })

# these numerical values - can be split 50:50 to define mice who have lost weight - i.e. - have the most curved body weight curve
sort(result_unlist) 

# cut into above and below median
summary(result_unlist)
f0_mothers_hi_gain <- names(result_unlist)[which(result_unlist >= median(result_unlist)) ]   
f0_mothers_lo_gain <- names(result_unlist)[which(result_unlist <= median(result_unlist)) ]

# mark mice in f0 table - based on hi and lo weight gain
mice_f0_females %<>% mutate(gain = case_when(
  any(animal_id %in% f0_mothers_lo_gain) == TRUE ~ "lo",
  any(animal_id %in% f0_mothers_hi_gain) == TRUE ~ "hi")) %>% 
  relocate(gain, .after = animal_id) %>% mutate(gain = as.factor(gain))

mice_f0_females$gain %>% summary()

# __b) Distinguish between offsprings whose mothers have or havn't lost weight ----

# get data for this objective 
mice_f1_slct_o1 <- mice_f1_slct

# mark mice in f1 table - based on hi and lo weight gain
mice_f1_slct_o1 <- left_join(mice_f1_slct_o1, {mice_f0_females %>% dplyr::select(animal_id, gain) %>% distinct}, by = c("mother_id" =  "animal_id"))
summary(mice_f1_slct_o1$gain)
mice_f1_slct_o1 %<>% rename(mother_gain = gain) 

# Is analysis of mother_gain possible - is it defined at all? 
mice_f1_slct_o1 %>% pull(mother_gain) %>% summary() # only low mice there - no data to analyse? 



# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "010_r_define_obesity.RData"))
renv::snapshot()


