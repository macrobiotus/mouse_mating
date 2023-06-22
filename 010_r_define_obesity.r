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

library("modeest")

library("lme4")     # Linear Mixed Effects Models
library("lmerTest") # Tests in Linear Mixed Effects Models
library("mgcv")     # General Additive Model

library("effects")     # Model inspection
library("performance") # Model inspection
library("cAIC4")       # Model selection 

library("gtsummary")

# _3.) Functions ----

# Calculate derivatives of polynomials ----
# - from https://stackoverflow.com/questions/40438195/function-for-polynomials-of-arbitrary-order-symbolic-method-preferred/40442584#40442584

g <- function (x, pc, nderiv = 0L) {
  ## check missing aruments
  if (missing(x) || missing(pc)) stop ("arguments missing with no default!")
  ## polynomial order p
  p <- length(pc) - 1L
  ## number of derivatives
  n <- nderiv
  ## earlier return?
  if (n > p) return(rep.int(0, length(x)))
  ## polynomial basis from degree 0 to degree `(p - n)`
  X <- outer(x, 0:(p - n), FUN = "^")
  ## initial coefficients
  ## the additional `+ 1L` is because R vector starts from index 1 not 0
  beta <- pc[n:p + 1L]
  ## factorial multiplier
  beta <- beta * factorial(n:p) / factorial(0:(p - n))
  ## matrix vector multiplication
  drop(X %*% beta)
}

polynom_curvature <- function (x, pc) {
  d1 <- g(x, pc, 1L)           # 1st derivative
  d2 <- g(x, pc, 2L)           # 2nd derivative
  d2 / (1 + d1 * d1) ^ (3 / 2) # curvature: 
}


# Data read-in, cleaning, formatting ----

# _1.) Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct.rds"))

# _2.) Inspect data ----

# Plot f0 and f1  weights by measurement day 

xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, type = "b", sub="f0 weight at measuerment age")
xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, type = "b", sub="f1 weight at measuerment age")


# Define Obesity ----

# _1.) Model 2nd degree polynomials ----

# not using poly() to coomodate function calculateing derivatives
F0_BodyWeight_Models <-  mice_f0_slct %>% group_by(AnimalId) %>% do(model = lm(BodyWeightG ~ as.numeric(MeasurementDay)+I(as.numeric(MeasurementDay)^2), data = .))
F1_BodyWeight_Models <-  mice_f1_slct %>% group_by(AnimalId) %>% do(model = lm(BodyWeightG ~ as.numeric(MeasurementDay)+I(as.numeric(MeasurementDay)^2), data = .)) 

# _2.) Inspect models ----

# - not all models are perfect, or the polynomial warranted, but as expected from curves
lapply(F0_BodyWeight_Models$model, summary) 
lapply(F1_BodyWeight_Models$model, summary) 


# - graphically, can only be done in lattice using poly() but result shoul be similar
xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, strip = FALSE,
      aspect = "xy", pch = 16, grid = TRUE,
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ...)
        fm <- lm(y ~ poly(x, 2))
        panel.lines(x, fitted(fm), col.line = "black")
      },
      xlab = "Standardized age", ylab = "Height (cm)")

# - graphically, can only be done in lattice using poly() but result shoul be similar
xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, strip = FALSE,
       aspect = "xy", pch = 16, grid = TRUE,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         fm <- lm(y ~ poly(x, 2))
         panel.lines(x, fitted(fm), col.line = "black")
       },
       xlab = "Standardized age", ylab = "Height (cm)")

# _3.) Isolate values for curvature calculations ----

# For function `polynom_curvature`

# - Isolate coefficients from model lists 

F0_BodyWeight_Models_Coefficients <- lapply(F0_BodyWeight_Models[[2]], coefficients) 
F1_BodyWeight_Models_Coefficients <- lapply(F1_BodyWeight_Models[[2]], coefficients) 

# - Isolate measurement days 

F0_BodyWeight_Models_MDays <- split(as.numeric(mice_f0_slct$MeasurementDay), mice_f0_slct$AnimalId)
F1_BodyWeight_Models_MDays <- split(as.numeric(mice_f1_slct$MeasurementDay), mice_f1_slct$AnimalId)

# _4.) Get curvatures of polynomials at each time point ----

# With function `polynom_curvature`
F0_PCurvature <- Map(polynom_curvature, F0_BodyWeight_Models_MDays, F0_BodyWeight_Models_Coefficients)
F1_PCurvature <- Map(polynom_curvature, F1_BodyWeight_Models_MDays, F1_BodyWeight_Models_Coefficients)

# _5.) Sum curvature of animals' polynomials across all time point ----

# Smallest values have most curved grow curve
F0_PCurvature_Results <- lapply(F0_PCurvature, sum) %>% unlist 
F1_PCurvature_Results <- lapply(F1_PCurvature, sum) %>% unlist 

# _6.) Check if curvature calculatoins make sense

# Show animal's id's and summed curvatures
sort(F0_PCurvature_Results)
sort(F1_PCurvature_Results)

# Plot f0 and f1 weights by measurement day

# Including sums of 2nd curvatures
xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, type = "b", sub="f0 weight at measuerment age, inlcuding curvature summary",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(80,17, cex = 0.75, labels = signif(F0_PCurvature_Results[panel.number()]), digits = 4) })

xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, type = "b", sub="f1 weight at measuerment age, inlcuding curvature summary",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(80,17, cex = 0.75, labels = signif(F1_PCurvature_Results[panel.number()]), digits = 4) })

# Separate obese from non-obese mice ----

# _1.) Finding modes for splitting data in obese and non-obese ----

plot(density(F0_PCurvature_Results), main = "F0 weight gain curvatures - Estimated mode")
abline(v = naive(F0_PCurvature_Results,bw = 0.001), col = 2)

plot(density(F1_PCurvature_Results), main = "F1 weight gain curvatures - Estimated mode")
abline(v = naive(F1_PCurvature_Results, bw = 0.01), col = 2)

# Storing modes

F0_Mode <- naive(F0_PCurvature_Results,bw = 0.001)
F1_Mode <- naive(F1_PCurvature_Results, bw = 0.01)

# _2.) Defining F0 animals with high or low weight gain ----

# Not yet obesity statuts - see methods

F0_HiGain <- names(F0_PCurvature_Results)[which(F0_PCurvature_Results >= F0_Mode) ]   
F0_LoGain <- names(F0_PCurvature_Results)[which(F0_PCurvature_Results < F0_Mode) ]   

# _3.) Defining F1 animals with high or low weight gain 

# Not yet obesity status - see methods

F1_HiGain <- names(F1_PCurvature_Results)[which(F1_PCurvature_Results >= F1_Mode) ]   
F1_LoGain <- names(F1_PCurvature_Results)[which(F1_PCurvature_Results < F1_Mode) ]   

# _4.) Mark mice in f0 table ----

# based on hi and low weight gain - can't get to work mutate and case_when

mice_f0_slct["WeightGain"] <- NA
mice_f0_slct[which(as.character(mice_f0_slct[["AnimalId"]]) %in% F0_HiGain), ]["WeightGain"] <- "hi"
mice_f0_slct[which(as.character(mice_f0_slct[["AnimalId"]]) %in% F0_LoGain), ]["WeightGain"] <- "lo"
mice_f0_slct[["WeightGain"]] <- as.factor(mice_f0_slct[["WeightGain"]])
mice_f0_slct$WeightGain %>% summary()

# _5.) Mark mice in f1 table ----

# based on hi and low weight gain - can't get to work mutate and case_when

mice_f1_slct["WeightGain"] <- NA
mice_f1_slct[which(as.character(mice_f1_slct[["AnimalId"]]) %in% F1_HiGain), ]["WeightGain"] <- "hi"
mice_f1_slct[which(as.character(mice_f1_slct[["AnimalId"]]) %in% F1_LoGain), ]["WeightGain"] <- "lo"
mice_f1_slct[["WeightGain"]] <- as.factor(mice_f1_slct[["WeightGain"]])
mice_f1_slct$WeightGain %>% summary()

# _6.) Define obese f0 mice ----

# Mark f0 with any HFD (also mixed diet) and hi weight gain are classified obese - all others are not obese! 

mice_f0_slct %<>% mutate(Obesity = case_when(
  (Diet == "HFD" & WeightGain == "hi") ~ "Obese",
  (Diet == "Mix" & WeightGain == "hi") ~ "Obese",
  (Diet == "CD" & WeightGain == "lo") ~ "NotObese", 
  TRUE ~ "NotObese"))
mice_f0_slct[["Obesity"]] <- as.factor(mice_f0_slct[["Obesity"]])
mice_f0_slct$Obesity %>% summary()

# _7.) Define obese f1 mice ----

# Mark f1 with  hi weight gain as obese - all others are not obese! 

mice_f1_slct %<>% mutate(Obesity = case_when(
  (WeightGain == "hi") ~ "Obese",
  (WeightGain == "lo") ~ "NotObese", 
  TRUE ~ "NotObese"))
mice_f1_slct[["Obesity"]] <- as.factor(mice_f1_slct[["Obesity"]])
mice_f1_slct$Obesity %>% summary()

# _8.) Add f0 obesity status to f1 data ----

# ***Weight data is available for mothers - using HFD AND weight gain as  proxy variables for obesity ***

F0_ObeseMothers <- mice_f0_slct %>% filter(Obesity == "Obese" & AnimalSex == "f") %>% pull("AnimalId")  %>% unique

# **** No weight data appears to be available for fathers - using HFD as a proxy variable for obesity ****

F0_ObeseFathers <- mice_f0_slct %>% filter(PartnerDiet == "HFD" ) %>% pull("MatingWith")  %>% unique

# - adding mothers' obesity status to F1

mice_f1_slct["ObeseMother"] <- FALSE
mice_f1_slct[which(as.character(mice_f1_slct[["MotherId"]]) %in% F0_ObeseMothers), ]["ObeseMother"] <- TRUE
mice_f1_slct[["ObeseMother"]] <- as.logical(mice_f1_slct[["ObeseMother"]])
mice_f1_slct$ObeseMother %>% summary()

# - adding farthers' obesity status to F1

mice_f1_slct["ObeseFather"] <- FALSE
mice_f1_slct[which(as.character(mice_f1_slct[["FatherId"]]) %in% F0_ObeseFathers), ]["ObeseFather"] <- TRUE
mice_f1_slct[["ObeseFather"]] <- as.logical(mice_f1_slct[["ObeseFather"]])
mice_f1_slct$ObeseFather %>% summary()

# - adding parents' obesity status to F1

mice_f1_slct %<>% mutate(ObeseParents = case_when(
  (ObeseMother == TRUE  & ObeseFather == FALSE) ~ "MotherObese",
  (ObeseMother == FALSE & ObeseFather == TRUE)  ~ "FatherObese",
  (ObeseMother == TRUE  & ObeseFather == TRUE)  ~ "MotherFatherObese",
  (ObeseMother == FALSE & ObeseFather == FALSE) ~ "MotherFatherNotObese",
  TRUE ~ NA))
mice_f1_slct[["ObeseParents"]] <- as.factor(mice_f1_slct[["ObeseParents"]])
mice_f1_slct$ObeseParents %>% summary()

# Summarizing final data ----

# _1.) Plotting f0 weight at measurement age, including sex and obesity status ----

xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, type = "b", sub="f0 weight at measurement age, inlcuding sex and obesity status",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(80,18, cex = 0.75, labels = mice_f0_slct$AnimalSex[panel.number()])
         panel.text(80,17, cex = 0.75, labels = mice_f0_slct$Obesity[panel.number()]) })

# _2.) Plotting f1 weight at measurement age, including sex and obesity status ----

xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, type = "b", sub="f1 weights at measurement ages, with sex and obesity status, and parents obesity status",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(80,13, cex = 0.75, labels = mice_f1_slct$AnimalSex[panel.number()])
         panel.text(80,16, cex = 0.75, labels = mice_f1_slct$Obesity[panel.number()])
         panel.text(80,19, cex = 0.75, labels = mice_f1_slct$ObeseParents[panel.number()])
         })

# _3.) Table summaries ----

mice_f0_slct %>% 
  select(Group, BodyWeightG, MeasurementDay, AnimalSex, WeightGain, Obesity) %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/010_r_define_obesity__mice_f0_slct__summary.docx") 


mice_f1_slct %>% 
  select(BodyWeightG, MeasurementDay, AnimalSex, WeightGain, Obesity, MotherDiet, FatherDiet, ObeseMother, ObeseFather) %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/010_r_define_obesity__mice_f1_slct__summary.docx") 

# How many HFD mothers were not obese?
mice_f0_slct %>% select(Group, Obesity) %>% count(Group, Obesity, sort = TRUE)

# Save finished data ----
saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

# Snapshot environment ----
sessionInfo()
save.image(file = here("scripts", "010_r_define_obesity.RData"))
renv::snapshot()

