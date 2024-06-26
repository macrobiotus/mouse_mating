#' ---
#' title: "Mice Mating Study"
#' subtitle: "Defining obesity in mice"
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
library("tibble") 
library("hablar") 
library("vctrs")

library("plotrix")
library("lattice")
library("ggplot2")
library("ggrepel")
library("ggpubr")

library("modeest")

library("lme4")     # Linear Mixed Effects Models
library("lmerTest") # Tests in Linear Mixed Effects Models
library("mgcv")     # General Additive Model

library("effects")     # Model inspection
library("performance") # Model inspection
library("cAIC4")       # Model selection 

library("gtsummary")

#' ## Functions

# _3.) Functions ----

'%!in%' <- function(x,y)!('%in%'(x,y))

# Calculate derivatives of polynomials 
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

#' # Data read-in, cleaning, formatting

# Data read-in, cleaning, formatting ----

# _1.) Get data ----

mice_f0_slct <- readRDS(file = here("rds_storage", "mice_f0_slct.rds"))
mice_f1_slct <- readRDS(file = here("rds_storage", "mice_f1_slct.rds"))

# _2.) Inspect data ----

# __a) Show F0 Information ----

mice_f0_slct %>% dplyr::select(AnimalId, AnimalSex, Diet, MatingWith, PartnerDiet) %>% distinct 

# __b) Show F1 Information ----

mice_f1_slct %>% 
  dplyr::select(AnimalId, AnimalSex, MotherDiet, FatherDiet) %>% 
  distinct %>% 
  arrange(AnimalSex, MotherDiet,FatherDiet) %>% print(n=Inf)


# _3.) Remove week 16 from F1 ----

# - Many NAs may affect growth curve analysis by yielding 0
# - F0 has data only available to week 14


# __a) Check for missing Measurement Days ----

mice_f0_slct %>% dplyr::select(AnimalId,AnimalSex, Week, BodyWeightG) %>% arrange(AnimalId) %>% print(n = Inf)
mice_f1_slct %>% dplyr::select(AnimalId,AnimalSex, Week, BodyWeightG) %>% arrange(AnimalId) %>% print(n = Inf)

# __b) Remove week 16 from F1

mice_f1_slct %<>% filter(Week != 16)

# __c) Check for missing Measurement Days ----

mice_f0_slct %>% dplyr::select(AnimalId,AnimalSex, Week, BodyWeightG) %>% arrange(AnimalId) %>% print(n = Inf)
mice_f1_slct %>% dplyr::select(AnimalId,AnimalSex, Week, BodyWeightG) %>% arrange(AnimalId) %>% print(n = Inf)

# __d) Plot F0 and F1  weights by measurement day ----

mice_f0_slct_xyplot <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, type = "b", sub="F0 weight at measurement age")
mice_f1_slct_xyplot <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, type = "b", sub="F1 weight at measurement age")

mice_f0_slct_xyplot
mice_f1_slct_xyplot

mice_slct_xyplots <- ggarrange(mice_f0_slct_xyplot, mice_f1_slct_xyplot, labels = c("a", "b"), ncol = 1, nrow = 2, heights = c(1.1, 1.9))

ggsave(plot = mice_slct_xyplots, path = here("../manuscript/display_items"), filename = "010_r_define_obesity__mice_weights.pdf")

#' # Define Obesity

# Define Obesity ----

# _1.) Model 2nd degree polynomials ----

# not using poly() to accommodate function calculating derivatives
F0_BodyWeight_Models <-  mice_f0_slct %>% group_by(AnimalId) %>% do(model = lm(BodyWeightG ~ as.numeric(MeasurementDay)+I(as.numeric(MeasurementDay)^2), data = .))
F1_BodyWeight_Models <-  mice_f1_slct %>% group_by(AnimalId) %>% do(model = lm(BodyWeightG ~ as.numeric(MeasurementDay)+I(as.numeric(MeasurementDay)^2), data = .)) 

# _2.) Inspect models ----

# __a) Show model summaries as text ----

# - not all models are perfect, or the polynomial warranted, but as expected from curves

# will clutter - uncomment only if needed

# lapply(F0_BodyWeight_Models$model, summary) 
# lapply(F1_BodyWeight_Models$model, summary) 

# __b) Show model summaries graphically ----

# - graphically, can only be done in lattice using poly() but result should be similar
mice_f0_slct_xyplot_poly <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct,
      aspect = "xy", pch = 16, grid = TRUE,
      panel = function(x, y, ...) {
        panel.xyplot(x, y, ...)
        fm <- lm(y ~ poly(x, 2))
        panel.lines(x, fitted(fm), col.line = "black")
      },
      xlab = "MeasurementDay", ylab = "BodyWeightG", sub="modelled F0 weight at measurement age")

# - graphically, can only be done in lattice using poly() but result shoul be similar
mice_f1_slct_xyplot_poly <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct,
       aspect = "xy", pch = 16, grid = TRUE,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         fm <- lm(y ~ poly(x, 2))
         panel.lines(x, fitted(fm), col.line = "black")
       },
       xlab = "MeasurementDay", ylab = "BodyWeightG", sub="modelled F1 weight at measurement age")

mice_f0_slct_xyplot_poly
mice_f1_slct_xyplot_poly

mice_slct_xyplots_poly <- ggarrange(mice_f0_slct_xyplot_poly, mice_f1_slct_xyplot_poly, labels = c("a", "b"), ncol = 1, nrow = 2, heights = c(1.1, 1.9))

ggsave(plot = mice_slct_xyplots_poly, path = here("../manuscript/display_items"), filename = "010_r_define_obesity__mice_weights_lm.pdf")

# _3.) Isolate values for curvature calculations ----

# For function `polynom_curvature`

# __a)  Isolate coefficients from model lists ----

F0_BodyWeight_Models_Coefficients <- lapply(F0_BodyWeight_Models[[2]], coefficients) 
F1_BodyWeight_Models_Coefficients <- lapply(F1_BodyWeight_Models[[2]], coefficients) 

# __b) Isolate measurement days ----

# inspect Measurement days
mice_f0_slct %>% dplyr::select(MeasurementDay, AnimalId, AnimalSex) %>% table()
mice_f1_slct %>% dplyr::select(MeasurementDay, AnimalId, AnimalSex) %>% table()

# isolate measurement days
F0_BodyWeight_Models_MDays <- split(as.numeric(mice_f0_slct$MeasurementDay), mice_f0_slct$AnimalId)
F1_BodyWeight_Models_MDays <- split(as.numeric(mice_f1_slct$MeasurementDay), mice_f1_slct$AnimalId)

# remove missing data
F0_BodyWeight_Models_MDays <- list_drop_empty(F0_BodyWeight_Models_MDays)
F1_BodyWeight_Models_MDays <- list_drop_empty(F1_BodyWeight_Models_MDays)

# check filtering
length(F0_BodyWeight_Models_MDays)
length(F1_BodyWeight_Models_MDays)



# _4.) Get curvatures (2nd derivative) of polynomials at each time point ----

# With function `polynom_curvature`
F0_PCurvature <- Map(polynom_curvature, F0_BodyWeight_Models_MDays, F0_BodyWeight_Models_Coefficients)
F1_PCurvature <- Map(polynom_curvature, F1_BodyWeight_Models_MDays, F1_BodyWeight_Models_Coefficients)

# _5.) Sum curvature of animals' polynomials across all time point ----

# Smallest values have most curved grow curve
F0_PCurvature_Results <- lapply(F0_PCurvature, sum) %>% unlist 
F1_PCurvature_Results <- lapply(F1_PCurvature, sum) %>% unlist 

# _6.) Check if curvature calculatoins make sense ----

# __a) Show animal's id's and summed curvatures ----
sort(F0_PCurvature_Results)
sort(F1_PCurvature_Results)

# __b) Plot F0 and F1 weights by measurement day ----

# Including sums of 2nd curvatures
mice_f0_slct_xyplot_curves <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, type = "b", sub="F0 weight at measuerment age, inlcuding curvature summary",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(80,17, cex = 0.75, labels = signif(F0_PCurvature_Results[panel.number()]), digits = 4) })

mice_f1_slct_xyplot_curves <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, type = "b", sub="F1 weight at measuerment age, inlcuding curvature summary",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(80,17, cex = 0.75, labels = signif(F1_PCurvature_Results[panel.number()]), digits = 4) })

# The more steep the curve the higher the shown number
mice_f0_slct_xyplot_curves
mice_f1_slct_xyplot_curves

mice_slct_xyplots_curves <- ggarrange(mice_f0_slct_xyplot_curves, mice_f1_slct_xyplot_curves, labels = c("a", "b"), ncol = 1, nrow = 2, heights = c(1.1, 1.9))

ggsave(plot = mice_slct_xyplots_curves, path = here("../manuscript/display_items"), filename = "010_r_define_obesity__mice_weights_curves.pdf")

#' # Separate obese from non-obese mice

# Separate obese from non-obese mice ----

# _1.) Finding a measure for splitting data in obese and non-obese ----

# __a) Trying to define a mode ----

# modes don't work well to separte bimodal distribution of second derivatives - will use medeian instead
plot(density(F0_PCurvature_Results), main = "F0 weight gain curvatures - Estimated mode")
abline(v = naive(F0_PCurvature_Results, bw = 0.001), col = 2)

plot(density(F1_PCurvature_Results), main = "F1 weight gain curvatures - Estimated mode")
abline(v = naive(F1_PCurvature_Results, bw = 0.001), col = 2)

# __b) Show distribution of 2nd derivates and the median which separates high gain from low gain mice ----

# Use this for reporting if necessary:
median(F0_PCurvature_Results)

mice_f0_derivatives_density <- ggplot(as_tibble_col(F0_PCurvature_Results, column_name = "SecondDerivative"), aes(x=SecondDerivative)) + 
  geom_density() +
  geom_vline(xintercept = median(F0_PCurvature_Results), color="blue") +
  xlab("second derivatives of growth curves") +
  ylab("denisty") +
  ggtitle("F0 distribution of growth curve derivatives, with median value") +
  theme_bw()

# Use this for reporting if necessary:
median(F1_PCurvature_Results)

mice_f1_derivatives_density <- ggplot(as_tibble_col(F1_PCurvature_Results, column_name = "SecondDerivative"), aes(x=SecondDerivative)) + 
  geom_density() +
  geom_vline(xintercept = median(F1_PCurvature_Results), color="blue") +
  xlab("second derivatives of growth curves") +
  ylab("denisty") +
  ggtitle("F1 distribution of growth curve derivatives, with median") +
  theme_bw()

mice_f0_derivatives_density
mice_f1_derivatives_density

mice_derivatives_densities <- ggarrange(mice_f0_derivatives_density, mice_f1_derivatives_density, labels = c("a", "b"), ncol = 1, nrow = 2, heights = c(1, 1))

ggsave(scale = 0.75, plot = mice_derivatives_densities, path = here("../manuscript/display_items"), filename = "010_r_define_obesity__mice_derivatives_densities.pdf")

# Storing modes and medians

F0_Mode <- naive(F0_PCurvature_Results,bw = 0.001)
F1_Mode <- naive(F1_PCurvature_Results, bw = 0.01)

F0_Median <- median(F0_PCurvature_Results)
F1_Median <- median(F1_PCurvature_Results)

# _2.) Defining F0 animals with high or low weight gain ----

# __a) Not used anymore: Define hi and lo gain via median of growth curv derivatives ---- 

# Does not give results that make any sense - makes only sens for F0
#  presumably because ist splits positive and negative values
F0_HiGain <- names(F0_PCurvature_Results)[which(F0_PCurvature_Results >= F0_Median) ]   
F0_LoGain <- names(F0_PCurvature_Results)[which(F0_PCurvature_Results < F0_Median) ]   

# __b) Used 30-Aug-2023e: Define hi and lo gain via positive or negative curvature summary ---- 
F0_HiGain <- names(F0_PCurvature_Results)[which(F0_PCurvature_Results >= 0) ]   
F0_LoGain <- names(F0_PCurvature_Results)[which(F0_PCurvature_Results < 0) ]   

# __c) Inspect results ----

# HiGian F0 are
F0_HiGain

# LoGain F0 are
F0_LoGain

# _3.) Defining F1 animals with high or low weight gain ----

# __a) Not used anymore: Define hi and lo gain via median of growth curv derivatives ---- 

# Does not give results that make any sense - makes only sens for F0
#  presumably because ist splits positive and negative values
F1_HiGain <- names(F1_PCurvature_Results)[which(F1_PCurvature_Results >= F1_Median) ]   
F1_LoGain <- names(F1_PCurvature_Results)[which(F1_PCurvature_Results < F1_Median) ]   

# __b) Used 30-Aug-2023e: Define hi and lo gain via positive or negative curvature summary ---- 
F1_HiGain <- names(F1_PCurvature_Results)[which(F1_PCurvature_Results >= 0) ]   
F1_LoGain <- names(F1_PCurvature_Results)[which(F1_PCurvature_Results < 0) ]   

# __c) Inspect results ----

# HiGian F1 are
F1_HiGain

# LoGain F1 are
F1_LoGain

# _4.) Mark obese mice in F0 table ----

# __a) Mark animals ----

mice_f0_slct["WeightGain"] <- NA
mice_f0_slct[which(as.character(mice_f0_slct[["AnimalId"]]) %in% F0_HiGain), ]["WeightGain"] <- "hi"
mice_f0_slct[which(as.character(mice_f0_slct[["AnimalId"]]) %in% F0_LoGain), ]["WeightGain"] <- "lo"
mice_f0_slct[["WeightGain"]] <- as.factor(mice_f0_slct[["WeightGain"]])

# __b) Check successful marking ----

mice_f0_slct %>% dplyr::select(AnimalId) %>% distinct # next command below should have 12 animals
mice_f0_slct %>% dplyr::select(AnimalId, WeightGain) %>% distinct # 12 animals are marked as expected

# __c ) Export table with successful marking ----

mice_f0_slct %>% dplyr::select(AnimalId, WeightGain) %>% distinct %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(here("../manuscript/display_items/"), "010_r_define_obesity__mice_f0_slct__weight_gain.docx")) 

# _5.) Mark mice in F1 table ----

# __a) Mark animals ----

mice_f1_slct["WeightGain"] <- NA
mice_f1_slct[which(as.character(mice_f1_slct[["AnimalId"]]) %in% F1_HiGain), ]["WeightGain"] <- "hi"
mice_f1_slct[which(as.character(mice_f1_slct[["AnimalId"]]) %in% F1_LoGain), ]["WeightGain"] <- "lo"
mice_f1_slct[["WeightGain"]] <- as.factor(mice_f1_slct[["WeightGain"]])

# __b) Check successful marking ----

mice_f1_slct %>% dplyr::select(AnimalId) %>% distinct # next command below should have 50 animals
mice_f1_slct %>% dplyr::select(AnimalId, WeightGain) %>% distinct # 50 animals are marked as expected

# __c ) Export table with successful marking ----

mice_f1_slct %>% dplyr::select(AnimalId, WeightGain) %>% distinct %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(here("../manuscript/display_items/"), "010_r_define_obesity__mice_f1_slct__weight_gain.docx")) 

# _6.) Define obese F0 mice ----

# __a) Mark F0 with any HFD and hi weight gain as classified obese and all others are not obese! ----

mice_f0_slct %<>% mutate(Obesity = case_when(
  (Diet == "HFD" & WeightGain == "hi") ~ "Obese",
  (Diet == "CD" & WeightGain == "lo") ~ "NotObese", 
  TRUE ~ "NotObese"))
mice_f0_slct[["Obesity"]] <- as.factor(mice_f0_slct[["Obesity"]])

# __b) Check successful marking ----

mice_f0_slct %>% dplyr::select(AnimalId) %>% distinct # next command below should have 8 animals

# One F0 mouse with HFD did not get obese per our definition - 8992

mice_f0_slct %>% dplyr::select(AnimalId, Diet, WeightGain, Obesity) %>% distinct %>% arrange(AnimalId) 

# __c) Export table with successful obesity marking ----

mice_f0_slct %>% dplyr::select(AnimalId, Diet, WeightGain, Obesity) %>% distinct %>% arrange(AnimalId) %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(here("../manuscript/display_items/"), "010_r_define_obesity__mice_f0_slct__obesity.docx")) 

# _7.) Define obese F1 mice ----

# __a) Mark F1 with hi weight gain as obese and all others are not obese ----

mice_f1_slct %<>% mutate(Obesity = case_when(
  (WeightGain == "hi") ~ "Obese",
  (WeightGain == "lo") ~ "NotObese", 
  TRUE ~ "NotObese"))
mice_f1_slct[["Obesity"]] <- as.factor(mice_f1_slct[["Obesity"]])
mice_f1_slct$Obesity %>% summary()

# __b) Check successful marking ----

mice_f1_slct %>% dplyr::select(AnimalId) %>% distinct # next command below should have 50 animals

# obesity marking successful as object has 50 lines 

mice_f1_slct %>% dplyr::select(AnimalId, MotherDiet, FatherDiet,  WeightGain, Obesity) %>% distinct %>% arrange(AnimalId) 

# __c) Export table with successful obesity marking ----

# table with summaries
mice_f1_slct %>% dplyr::select(AnimalId, MotherDiet, FatherDiet,  WeightGain, Obesity) %>% distinct %>% arrange(AnimalId) %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(here("../manuscript/display_items/"), "010_r_define_obesity__mice_f1_slct__obesity_table.docx")) 

# plain table for sanity
mice_f1_slct %>% dplyr::select(AnimalId, MotherDiet, FatherDiet,  WeightGain, Obesity) %>% distinct %>% arrange(AnimalId) %>% 
  openxlsx::write.xlsx(., file = paste0(here("../manuscript/display_items/"), "010_r_define_obesity__mice_f1_slct__obesity.xlsx"), asTable = TRUE) 

# _8.) Add F0 obesity status to F1 data ----

# __a) Check F0 Obesity status ----

# ***Weight data is available for mothers - using HFD AND weight gain as  proxy variables for obesity ***

F0_ObeseMothers <- mice_f0_slct %>% filter(Obesity == "Obese" & AnimalSex == "f") %>% pull("AnimalId")  %>% unique

# **** No weight data appears to be available for fathers - using HFD as a proxy variable for obesity ****
mice_f0_slct %>% filter(Obesity == "Obese" & AnimalSex == "m") %>% pull("AnimalId")  %>% unique

F0_ObeseFathers <- mice_f0_slct %>% filter(PartnerDiet == "HFD" ) %>% pull("MatingWith")  %>% unique

# __b) Adding F0 mothers' obesity status to F1 ----

mice_f1_slct["ObeseMother"] <- FALSE
mice_f1_slct[which(as.character(mice_f1_slct[["MotherId"]]) %in% F0_ObeseMothers), ]["ObeseMother"] <- TRUE
mice_f1_slct[["ObeseMother"]] <- as.logical(mice_f1_slct[["ObeseMother"]])
mice_f1_slct$ObeseMother %>% summary()

# __c) Adding F0 fathers' obesity status to F1 ----

mice_f1_slct["ObeseFather"] <- FALSE
mice_f1_slct[which(as.character(mice_f1_slct[["FatherId"]]) %in% F0_ObeseFathers), ]["ObeseFather"] <- TRUE
mice_f1_slct[["ObeseFather"]] <- as.logical(mice_f1_slct[["ObeseFather"]])
mice_f1_slct$ObeseFather %>% summary()

# __d) Adding parents' obesity status to F1 ----

mice_f1_slct %<>% mutate(ObeseParents = case_when(
  (ObeseMother == TRUE  & ObeseFather == FALSE) ~ "MotherObese",
  (ObeseMother == FALSE & ObeseFather == TRUE)  ~ "FatherObese",
  (ObeseMother == TRUE  & ObeseFather == TRUE)  ~ "MotherFatherObese",
  (ObeseMother == FALSE & ObeseFather == FALSE) ~ "MotherFatherNotObese",
  TRUE ~ NA))
mice_f1_slct[["ObeseParents"]] <- as.factor(mice_f1_slct[["ObeseParents"]])

# __e) Inspecting and exporting parents obesity status ----

mice_f1_slct %>% dplyr::select(AnimalId, MotherDiet, FatherDiet,  WeightGain, Obesity, ObeseParents) %>% distinct %>% arrange(AnimalId) %>% print(n = Inf)

mice_f1_slct %>% dplyr::select(AnimalId, MotherDiet, FatherDiet,  WeightGain, Obesity, ObeseParents) %>% distinct %>% arrange(AnimalId) %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = paste0(here("../manuscript/display_items/"), "010_r_define_obesity__mice_f1_slct__parents_obesity.docx")) 

#' # Summarizing final data

# Summarizing final data ----

# _1.) Plotting F0 weight at measurement age, including sex and obesity status ----

mice_f0_slct_xyplot_final <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f0_slct, type = "b", sub="F0 weight at measurement age, inlcuding sex and obesity status",
       panel = function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(45,27, cex = 0.75, labels = mice_f0_slct$AnimalSex[panel.number()])
         panel.text(51,27, cex = 0.75, labels = mice_f0_slct$Diet[panel.number()])
         panel.text(66,27, cex = 0.75, labels = mice_f0_slct$Obesity[panel.number()])
         panel.text(80,17, cex = 0.75, labels = signif(F0_PCurvature_Results[panel.number()]), digits = 4)
         }
       )

# _2.) Plotting F1 weight at measurement age, including sex and obesity status ----

# __a.) Compound figure, for manuscript ----

mice_f1_slct_xyplot_final <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = mice_f1_slct, type = "b", sub="F1 weights at measurement ages, with sex and obesity status, and parents obesity status",
       panel=function(x, y,...){
         panel.xyplot(x,y,...)
         panel.text(35,30, cex = 0.75, labels = mice_f1_slct$AnimalSex[panel.number()])
         panel.text(70,30, cex = 0.75, labels = mice_f1_slct$Obesity[panel.number()])
         panel.text(70,12, cex = 0.75, labels = mice_f1_slct$ObeseParents[panel.number()])
         panel.text(80,17, cex = 0.75, labels = signif(F1_PCurvature_Results[panel.number()]), digits = 4)
         })

# __b.) Male offspring, for talks ----

mice_f1_slct_xyplot_final_m <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = (mice_f1_slct %>% dplyr::filter(AnimalSex == "m")), type = "b", sub="F1 **male** weights at measurement ages, with sex and obesity status, and parents obesity status",
                                      panel=function(x, y,...){
                                        panel.xyplot(x,y,...)
                                        panel.text(35,30, cex = 0.75, labels = (mice_f1_slct %>% dplyr::filter(AnimalSex == "m"))$AnimalSex[panel.number()])
                                        panel.text(70,30, cex = 0.75, labels = (mice_f1_slct %>% dplyr::filter(AnimalSex == "m"))$Obesity[panel.number()])
                                        panel.text(70,12, cex = 0.75, labels = (mice_f1_slct %>% dplyr::filter(AnimalSex == "m"))$ObeseParents[panel.number()])
                                        panel.text(80,17, cex = 0.75, labels = signif(F1_PCurvature_Results[panel.number()]), digits = 4)
                                      })


# __b.) Female offspring, for talks ----

mice_f1_slct_xyplot_final_f <- xyplot(BodyWeightG ~ MeasurementDay | AnimalId, data = (mice_f1_slct %>% dplyr::filter(AnimalSex == "f")), type = "b", sub="F1 **female** weights at measurement ages, with sex and obesity status, and parents obesity status",
                                      panel=function(x, y,...){
                                        panel.xyplot(x,y,...)
                                        panel.text(35,22, cex = 0.75, labels = (mice_f1_slct %>% dplyr::filter(AnimalSex == "f"))$AnimalSex[panel.number()])
                                        panel.text(70,22, cex = 0.75, labels = (mice_f1_slct %>% dplyr::filter(AnimalSex == "f"))$Obesity[panel.number()])
                                        panel.text(70,12, cex = 0.75, labels = (mice_f1_slct %>% dplyr::filter(AnimalSex == "f"))$ObeseParents[panel.number()])
                                        panel.text(80,17, cex = 0.75, labels = signif(F1_PCurvature_Results[panel.number()]), digits = 4)
                                      })






# _3.) Revision started 30.08.2023 - Adding new plots with weight deltas of F0 ----

# __a) Copy object for new plot ----

mice_f0_slct_mb  <- mice_f0_slct

# __b) Get weight delta ----

# getting time points for delta - min and max wont work as some mice don't ahve data at max weeks
mice_f0_slct_mb %>% hablar::convert(num(Week)) %>% group_by(AnimalId) %>% dplyr::slice(c(which.min(Week), which.max(Week))) %>% arrange 

# show available weeks:
mice_f0_slct_mb %>% convert(num(Week)) %>% pull(Week) %>% unique

# getting time points for delta - min week and week 14
mice_f0_slct_mb %<>% hablar::convert(num(Week)) %>% group_by(AnimalId) %>% dplyr::slice(c(which.min(Week), which(Week == 14)))

# calculate weight gain by subtracting weight at week 14 from weight at week 4 - add this to caption
mice_f0_slct_mb %<>% group_by(AnimalId) %>% mutate(BodyWeightGainDeltaG = BodyWeightG - first(BodyWeightG)) %>% relocate(BodyWeightGainDeltaG, .after = BodyWeightG)

# keep only rows with the relevant weight gain delta - the second column of each group
mice_f0_slct_mb %<>% group_by(AnimalId) %>% dplyr::slice(min(n(), 2))

# sanity check - compare to `010_r_define_obesity__mice_f1_slct__obesity.xlsx` for F1 exported above, should look similar
mice_f0_slct_mb %>% 
  dplyr::select(AnimalId, AnimalSex, BodyWeightGainDeltaG,  Diet,  WeightGain, Obesity) %>% 
  distinct %>% 
  arrange(AnimalSex, BodyWeightGainDeltaG) %>%
  print(n = Inf)

# __c) Prepare plotting ----

# # encode sex into the factor variable to label ggplot 2 facets without much work
# mice_f0_slct_mb %<>% 
#   convert(chr(MotherDiet, FatherDiet)) %>%
#   mutate(MotherDiet = paste0("Father ", MotherDiet)) %>%
#   mutate(FatherDiet = paste0("Mother ", FatherDiet)) %>%
#   convert(chr(MotherDiet, FatherDiet))

# __d) Plot weight delta ----

f0_mice_weights_sex_deltas <- ggplot(data = mice_f0_slct_mb, aes(x = AnimalSex, y = BodyWeightGainDeltaG)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), aes(shape = WeightGain), color = "black", size = 3) +
  geom_boxplot(width = 0.2, alpha = 0.2) +
  facet_wrap(Diet ~ . , ncol = 2) + 
  coord_cartesian(ylim = c(0, 16)) +
  theme_bw(base_size = 12) +
  labs(x = "F0 animal sex", y = "F0 weight gain Δ [g]")
NULL

f0_mice_weights_sex_deltas

# _4.) Revision started 22.08.2023 - Adding new plots with weight deltas of F1 ----

# __a) Copy object for new plot ----

mice_f1_slct_mb  <- mice_f1_slct

# __b) Get weight delta ----

# getting time points for delta - min and max wont work as some mice don't ahve data at max weeks
mice_f1_slct_mb %>% convert(num(Week)) %>% group_by(AnimalId) %>% slice(c(which.min(Week), which.max(Week))) %>% arrange 

# show available weeks:
mice_f1_slct_mb %>% convert(num(Week)) %>% pull(Week) %>% unique

# getting time points for delta - min week and week 14
mice_f1_slct_mb %<>% convert(num(Week)) %>% group_by(AnimalId) %>% slice(c(which.min(Week), which(Week == 14)))

# calculate weight gain by subtracting weight at week 14 from weight at week 4 - add this to caption
mice_f1_slct_mb %<>% group_by(AnimalId) %>% mutate(BodyWeightGainDeltaG = BodyWeightG - first(BodyWeightG)) %>% relocate(BodyWeightGainDeltaG, .after = BodyWeightG)

# keep only rows with the relvant weight gein delta - the second column of each group
mice_f1_slct_mb %<>% group_by(AnimalId) %>% slice(min(n(), 2))

# sanity check - compare to `010_r_define_obesity__mice_f1_slct__obesity.xlsx` expoted above, shoul look similar
mice_f1_slct_mb %>% 
  dplyr::select(AnimalId, AnimalSex, BodyWeightGainDeltaG,  MotherDiet, FatherDiet,  WeightGain, Obesity) %>% 
  distinct %>% 
  arrange(AnimalSex, BodyWeightGainDeltaG) %>%
  print(n = Inf)

# __c) Prepare plotting ----

# encode sex into the factor variable to label ggplot 2 facets without much work
mice_f1_slct_mb %<>% 
  convert(chr(MotherDiet, FatherDiet)) %>%
  mutate(MotherDiet = paste0("Father ", MotherDiet)) %>%
  mutate(FatherDiet = paste0("Mother ", FatherDiet)) %>%
  convert(chr(MotherDiet, FatherDiet))

# __d) Plot weight delta ----

f1_mice_weights_sex_deltas <- ggplot(data = mice_f1_slct_mb, aes(x = AnimalSex, y = BodyWeightGainDeltaG)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), aes(shape = WeightGain), color = "black", size=3) +
  geom_boxplot(width = 0.2, alpha = 0.2) +
  facet_wrap(MotherDiet ~ FatherDiet, ncol = 4) + 
  coord_cartesian(ylim = c(0, 15)) +
  theme_bw(base_size = 12) +
  labs(x = "F1 animal sex", y = "F1 weight gain Δ [g]")
  NULL

f1_mice_weights_sex_deltas

# _5.) Showing and exporting plot ----

# __a) Save weight curve plots ----

mice_f0_slct_xyplot_final
mice_f1_slct_xyplot_final

mice_slct_xyplots_obesity <- ggarrange(mice_f0_slct_xyplot_final, mice_f1_slct_xyplot_final, labels = c("a", "b"), ncol = 1, nrow = 2, heights = c(1, 2))

ggsave(plot = mice_slct_xyplots_obesity, scale = 1.2, path = here("../manuscript/display_items"), filename = "010_r_define_obesity__mice_weights_sex_obesity.pdf")

# __b) Save F0 and F1 weight delta plots ----

f0_f1_mice_weights_sex_deltas <- ggarrange(f0_mice_weights_sex_deltas, f1_mice_weights_sex_deltas, labels = c("a", "b"), ncol = 2, nrow = 1, widths =  c(2, 4))

ggsave(device = cairo_pdf, plot = f0_f1_mice_weights_sex_deltas, width = 210, height = 100, units = c("mm"), dpi = 300, scale = 1.2, path = here("../manuscript/display_items"), filename = "010_r_define_obesity__f0_f1_mice_weights_sex_deltas.pdf")

# _6.) Table summaries ----

# __a) Default tables - of original draft ----

mice_f0_slct %>% 
  dplyr::select(AnimalId, AnimalSex, WeightGain, Obesity) %>% distinct() %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/010_r_define_obesity__mice_f0_slct__summary.docx") 

mice_f1_slct %>% 
  dplyr::select(AnimalId, AnimalSex, WeightGain, Obesity, MotherDiet, FatherDiet, ObeseMother, ObeseFather) %>% distinct() %>% 
  tbl_summary(., missing = "no") %>%
  as_gt() %>%
  gt::gtsave(filename = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/010_r_define_obesity__mice_f1_slct__summary.docx") 

# How many HFD mothers were not obese?
mice_f0_slct %>% dplyr::select(AnimalId, Obesity) %>% distinct 
mice_f0_slct %>% dplyr::select(AnimalId, Obesity) %>% distinct %>% tbl_summary(., missing = "no")

# __b) Tables showing all information of section a) above, but also the relationship between mice ----

# not implemented here see script 40 and README.md 30-Aug-2023

#' #  Save finished data

# Save finished data ----

saveRDS(mice_f0_slct, file = here("rds_storage", "mice_f0_slct_with_obesity.rds"))
saveRDS(mice_f1_slct, file = here("rds_storage", "mice_f1_slct_with_obesity.rds"))

#' # Snapshot environment

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "010_r_define_obesity.RData"))
renv::snapshot()

