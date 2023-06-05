#' ---
#' title: "Mice Mating Study"
#' subtitle: "Arreay analysis provided by AH, adjusted by PC"
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

library(here)
library(renv)
library(magrittr)
library(ggplot2)
library(mgcv)
library(parallel, lib.loc = "/Users/paul/Library/Caches/org.R-project.R/R/renv/sandbox/R-4.2/aarch64-apple-darwin20/84ba8b13")
library(RVAideMemoire)
library("GCSscore")
library("DESeq2")
library("pheatmap")
library("oligoClasses") # to annotate array target
library("pd.clariom.s.mouse")
# for GO analysis
library("org.Hs.eg.db")
library("limma")
library(stringr)
library(tidyr)
library(Biobase)
library(BiocGenerics)
library(gplots)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(KernSmooth)
library(data.table)
# library(swamp)
library(factoextra)
library(FactoMineR)
library(UpSetR)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(usethis) 

# _3.) Increase meomeory if needed ----

# usethis::edit_r_environ()

# _4.) Functions ----

# rewrite metadata in expression data sets
adjust_array_data = function(expression_set, model_variables) {
  require(dplyr)
  require(Biobase)
  
  # adjust column and row names names in expression data
  colnames(expression_set) <- colnames(expression_set) %>%
    str_replace_all("bAT",    "BRAT") %>%
    str_replace_all("ingWAT", "SCAT") %>%
    str_replace_all("liver",  "LIAT") %>%
    str_replace_all("eWAT",   "EVAT")
  
  # adjust "Tissue" column values  in expression data
  pData(expression_set) %<>% mutate(
    Tissue = case_when(
      Tissue == "bAT"    ~ "BRAT",
      Tissue == "ingWAT" ~ "SCAT",
      Tissue == "liver"  ~ "LIAT",
      Tissue == "eWAT"   ~ "EVAT"
    )
  )
  
  # merge obesity variables from modelling  to metadata from array experiments
  pData(expression_set) <- left_join((pData(expression_set) %>% dplyr::select(
    -c(Sex, Parental_diet, Diet_mother, Diet_father, Group)
  )),
  (model_variables),
  by = c("Animal" = "AnimalId"))
  
  return(expression_set)
  
}

# get PCA plots
get_pca_plot = function(expr_data_pca, expr_data_raw, variable, legend_title, plot_title) {
  require("factoextra")
  require("ggplot2")
  
  pca_plot <- fviz_pca_ind(
    expr_data_pca,
    label = "none",
    habillage = factor(pData(expr_data_raw)[[variable]]),
    pointsize = 2,
    palette = c("firebrick3", "purple", "steelblue3", "gold3"),
    legend.title = legend_title,
    invisible = "none",
    addEllipses = FALSE,
    title = plot_title
  ) +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance")
    ) +
    theme_bw() +
    scale_shape_manual(values = c(19, 19, 19, 19))
  
  return(pca_plot)
  
}

# Getting model predictions with mgcv
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

# Convenience function to test Offsprings obesity status in all tissue samples
get_dge_for_offspring_obesity = function(ExpSet){
  
  # message
  message("Convenience function to test Offsprings obesity status in all tissue samples - analysis of factor \"ObesityLgcl\" is hard-coded.")
  
  # Checking sample size
  message("Consider whether the sample size is appropraite for DGE, you shoul have at least six sample per group:")
  pData(ExpSet) %>% dplyr::select(ObesityLgcl) %>% table() %>% print()
  
  # Creating design matrix and contrasts
  design_offspr_obese <- model.matrix(~ ExpSet[["ObesityLgcl"]] - 1)
  colnames(design_offspr_obese) <-c("NotObese", "Obese") 
  contrast_offspr_obese <- makeContrasts("NotObese-Obese", levels = design_offspr_obese)
  
  # Creating and evaluating linear model of gene expression
  fit1 <- lmFit(ExpSet, design_offspr_obese)
  fit1c <- contrasts.fit(fit1, contrast_offspr_obese)
  fit1c <- eBayes(fit1c)
  
  # message
  message("Returning topTable() output:")
  
  return(topTable(fit1c, p.value = 0.05, number = 50))
  
}

# Convenience function to test parents' obesity status in all tissue samples
get_dge_for_parent_obesity = function(ExpSet){
  
  message("Convenience function to test parents' obesity status in all tissue samples - analysis of factor \"ObeseParents\" is hard-coded.")
  
  message("Consider whether the sample size is appropriate for DGE, you shoul have at least six samples per group:")
  
  pData(ExpSet) %>% dplyr::select(ObeseParents) %>% table() %>% print()
  
  message(paste0("Building initial model on data set \"", deparse(substitute(ExpSet))), "\".")   
  
  design_parent_obese <- model.matrix(~ ExpSet[["ObeseParents"]] - 1)
  colnames(design_parent_obese) <-c("MotherFatherNotObese", "FatherObese", "MotherFatherObese", "MotherObese") 
  fit_parent_obese <- lmFit(ExpSet, design_parent_obese)
  
  message("Defining and applying contrasts: One of \"MotherFatherNotObese\", \"FatherObese\", \"MotherFatherObese\", or \"MotherObese\" against all remening three levels.") 
  
  contrast_mfo  <- makeContrasts("MotherFatherObese" = MotherFatherObese - (MotherObese + FatherObese + MotherFatherNotObese)/3 , levels = design_parent_obese)
  contrast_fob   <- makeContrasts("FatherObese" = FatherObese - (MotherObese + MotherFatherObese + MotherFatherNotObese)/3 , levels = design_parent_obese)
  contrast_mob   <- makeContrasts("MotherObese" = MotherObese - (FatherObese + MotherFatherObese + MotherFatherNotObese)/3 , levels = design_parent_obese)
  contrast_nob <- makeContrasts("MotherFatherNotObese" = MotherFatherNotObese - (FatherObese + MotherFatherObese + MotherFatherObese)/3 , levels = design_parent_obese)
  
  fit_mfo <- contrasts.fit(fit_parent_obese, contrast_mfo)
  fit_fob <- contrasts.fit(fit_parent_obese, contrast_fob)
  fit_mob <- contrasts.fit(fit_parent_obese, contrast_mob)
  fit_nob <- contrasts.fit(fit_parent_obese, contrast_nob)
  
  fit_mfo_eb <- eBayes(fit_mfo)
  fit_fob_eb <- eBayes(fit_fob)
  fit_mob_eb <- eBayes(fit_mob)
  fit_nob_eb <- eBayes(fit_nob)
  
  message("returning results list - check list names for top table identification.") 
  
  return(
    list(
      "MotherFatherObese" = topTable(fit_mfo_eb, p.value = 0.05, number = 50), 
      "FatherObese" = topTable(fit_fob_eb, p.value = 0.05, number = 50),
      "MotherObese" = topTable(fit_mob_eb, p.value = 0.05, number = 50), 
      "MotherFatherNotObese" = topTable(fit_nob_eb, p.value = 0.05, number = 50)
    )
  )
}


# _5.) Color code for plotting ----

hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(n=11, name="RdBu"))(50))

# Load and shape data ----

# _1.) Loading normalized data as ExpressionSet type ----

# Data are Clariom S mouse arrays
# I removed strong outliers: A285_liver clustert zu bAT, A339_liver weit weg vom Rest

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/allTissues_normData.RData") # only if you are interested to look into the normalized data

# copy to stick to manuscript naming conventions
FLAT <- normData; rm(normData)

saveRDS(FLAT, file = here("rds_storage", "050_r_array_analysis__normalized_data.rds"))

# _2.) Loading normalized data as ExpressionSet type - normalised individually for each tissue ----

# Ich habe die normalisierung für jedes Gewebe getrennt für die DGE Analysen gemacht, 
# da dies genauer ist (Gewebe zu weit auseinander in PCA)

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/normData4DGE.RData") #für jedes Gewebe die normalisierten Daten

# copy to stick to manuscript naming conventions - corrcted as per AH 25.05.2023
BRAT <- bAT_normData; rm(bAT_normData)        # brown adipose tissue
SCAT <- ingWAT_normData; rm(ingWAT_normData)  # subcutabneous adipose tissue 
LIAT <- Liv_normData; rm(Liv_normData)        # liver adipose tissue
EVAT <- eWAT_normData; rm(eWAT_normData)      # visceral adipose tissue

saveRDS(BRAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_bat.rds"))
saveRDS(SCAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_ewat.rds"))
saveRDS(LIAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_liv.rds" ))
saveRDS(EVAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_ingwat.rds"))

# _3.) Loading metadata from modelling (obesity variables) ----

mice_f1_modeled_data_with_rna_seq_data <- readRDS(file = here("rds_storage", "040_r_h3__mice_f1_modeled_data_with_rna_seq_data.rds"))

# _4.) Adjust variable names and inspect data ----

FLAT <- adjust_array_data(FLAT, mice_f1_modeled_data_with_rna_seq_data)
BRAT <- adjust_array_data(BRAT, mice_f1_modeled_data_with_rna_seq_data)
SCAT <- adjust_array_data(SCAT, mice_f1_modeled_data_with_rna_seq_data)
LIAT <- adjust_array_data(LIAT, mice_f1_modeled_data_with_rna_seq_data)
EVAT <- adjust_array_data(EVAT, mice_f1_modeled_data_with_rna_seq_data)

# check, if you like, using
# pData(BRAT)
# exprs(BRAT)

# _5.) Loading AHs dietary DGE analysis results ----

# __a) Load data ----

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/DGELists.RData") #lädt alle DGE Tabellen

# __b) Define list with dietary variables ----

DGEL_diet <- list(bAT_CD_HFD_VS_CD_CD, bAT_HFD_CD_VS_CD_CD, bAT_HFD_CD_VS_CD_HFD, bAT_HFD_HFD_VS_CD_CD, bAT_HFD_HFD_VS_CD_HFD,    
 bAT_HFD_HFD_VS_HFD_CD, eWAT_CD_HFD_VS_CD_CD, eWAT_HFD_CD_VS_CD_CD, eWAT_HFD_CD_VS_CD_HFD, eWAT_HFD_HFD_VS_CD_CD, 
 eWAT_HFD_HFD_VS_CD_HFD, eWAT_HFD_HFD_VS_HFD_CD, ingWAT_CD_HFD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_HFD,
 ingWAT_HFD_HFD_VS_CD_CD, ingWAT_HFD_HFD_VS_CD_HFD, ingWAT_HFD_HFD_VS_HFD_CD, Liv_CD_HFD_VS_CD_CD, Liv_HFD_CD_VS_CD_CD,
 Liv_HFD_CD_VS_CD_HFD, Liv_HFD_HFD_VS_CD_CD, Liv_HFD_HFD_VS_CD_HFD, Liv_HFD_HFD_VS_HFD_CD) 

names(DGEL_diet) <-  c("BRAT_CD_HFD_VS_CD_CD", "BRAT_HFD_CD_VS_CD_CD", "BRAT_HFD_CD_VS_CD_HFD", "BRAT_HFD_HFD_VS_CD_CD", "BRAT_HFD_HFD_VS_CD_HFD",     
 "BRAT_HFD_HFD_VS_HFD_CD", "EVAT_CD_HFD_VS_CD_CD", "EVAT_HFD_CD_VS_CD_CD", "EVAT_HFD_CD_VS_CD_HFD", "EVAT_HFD_HFD_VS_CD_CD",  
 "EVAT_HFD_HFD_VS_CD_HFD", "EVAT_HFD_HFD_VS_HFD_CD", "SCAT_CD_HFD_VS_CD_CD", "SCAT_HFD_CD_VS_CD_CD", "SCAT_HFD_CD_VS_CD_HFD", 
 "SCAT_HFD_HFD_VS_CD_CD", "SCAT_HFD_HFD_VS_CD_HFD", "SCAT_HFD_HFD_VS_HFD_CD", "LIAT_CD_HFD_VS_CD_CD", "LIAT_HFD_CD_VS_CD_CD", 
 "LIAT_HFD_CD_VS_CD_HFD", "LIAT_HFD_HFD_VS_CD_CD", "LIAT_HFD_HFD_VS_CD_HFD", "LIAT_HFD_HFD_VS_HFD_CD")

# __c) **CONTINUE HERE IF RE-IMPLEMNETATION OF DGE FAILS ** Define list with obesity variables 

DGEL_obes <- DGEL_diet

# As per 01.06.2023 in README commit `5fd8790e5024bce8e05f885d08f219b1c736ef58`:
# modify or exchange dietary variables to contain information regarding offsprings and/or parental obesity as per last scripts overview and manuscript tasks

# names(DGEL_obes) <-  c( "foo", "bar")

# __d) Clean environment and save data ---- 

rm(bAT_CD_HFD_VS_CD_CD, bAT_HFD_CD_VS_CD_CD, bAT_HFD_CD_VS_CD_HFD, bAT_HFD_HFD_VS_CD_CD, bAT_HFD_HFD_VS_CD_HFD,    
   bAT_HFD_HFD_VS_HFD_CD, eWAT_CD_HFD_VS_CD_CD, eWAT_HFD_CD_VS_CD_CD, eWAT_HFD_CD_VS_CD_HFD, eWAT_HFD_HFD_VS_CD_CD, 
   eWAT_HFD_HFD_VS_CD_HFD, eWAT_HFD_HFD_VS_HFD_CD, ingWAT_CD_HFD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_HFD,
   ingWAT_HFD_HFD_VS_CD_CD, ingWAT_HFD_HFD_VS_CD_HFD, ingWAT_HFD_HFD_VS_HFD_CD, Liv_CD_HFD_VS_CD_CD, Liv_HFD_CD_VS_CD_CD,
   Liv_HFD_CD_VS_CD_HFD, Liv_HFD_HFD_VS_CD_CD, Liv_HFD_HFD_VS_CD_HFD, Liv_HFD_HFD_VS_HFD_CD)

saveRDS(DGEL_diet, file = here("rds_storage", "050_r_array_analysis__dge_lists_by_diet.rds"))

# Overall Principal Component Analysis ----

# _1.) Get PC for all tissues ----

PCA_FLAT <- prcomp( t (exprs(FLAT)), scale. = FALSE)
percentVar <- round(100*PCA_FLAT$sdev^2 / sum(PCA_FLAT$sdev^2), 1)

# _2.) Plot PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_flat_a <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a")

# Overall expression differences and obesity status among f1 offspring
plot_pca_flat_b <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b")

# Overall expression differences and obesity parental obesity status among f0
plot_pca_flat_c <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c")


# _3.) Combine and save overall plots ----
plot_pca_flat <- ggarrange(plot_pca_flat_a, plot_pca_flat_b, plot_pca_flat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_flat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_flat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_flat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_flat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# _4.) Plot PCs in 3D (only for tissue types so far)  ----

dataGG <- data.frame(PC1 = PCA_FLAT$x[,1], PC2 = PCA_FLAT$x[,2],PC3 = PCA_FLAT$x[,3])

plot_pca_sptial <- plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  add_markers(color = ~pData(FLAT)$Tissue,colors = c("firebrick3","purple","steelblue3","gold3"),
              symbol = ~pData(FLAT)$Tissue, symbols = c(19,19,19,19))

saveWidget(plot_pca_sptial, file = paste0(here("plots"),"/", "050_r_array_analysis__plot_pca_flat.html"), selfcontained = T, libdir = "lib")

# Tissue-Specific PCA ----

# _1.) Get PCs for BRAT ----

PCA_BRAT <- prcomp( t (exprs(BRAT)), scale. = FALSE)
percentVar <- round(100*PCA_BRAT$sdev^2 / sum(PCA_BRAT$sdev^2), 1)

# _2.) Plot PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_brat_a <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a")

# Overall expression differences and obesity status among f1 offspring
plot_pca_brat_b <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b")

# Overall expression differences and obesity parental obesity status among f0
plot_pca_brat_c <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c")

# _3.) Combine and save overall plots ----
plot_pca_brat <- ggarrange(plot_pca_brat_a, plot_pca_brat_b, plot_pca_brat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_brat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_brat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_brat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_brat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# _4.) Get PCs for SCAT ----

PCA_SCAT <- prcomp( t (exprs(SCAT)), scale. = FALSE)
percentVar <- round(100*PCA_SCAT$sdev^2 / sum(PCA_SCAT$sdev^2), 1)

# _5.) Plot PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_scat_a <- get_pca_plot(expr_data_pca = PCA_SCAT, expr_data_raw = SCAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a")

# Overall expression differences and obesity status among f1 offspring
plot_pca_scat_b <- get_pca_plot(expr_data_pca = PCA_SCAT, expr_data_raw = SCAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b")

# Overall expression differences and obesity parental obesity status among f0
plot_pca_scat_c <- get_pca_plot(expr_data_pca = PCA_SCAT, expr_data_raw = SCAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c")

# _6.) Combine and save overall plots ----

plot_pca_scat <- ggarrange(plot_pca_scat_a, plot_pca_scat_b, plot_pca_scat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_scat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_scat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_scat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_scat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# _7.) Get PCs for EVAT ----

PCA_EVAT <- prcomp( t (exprs(EVAT)), scale. = FALSE)
percentVar <- round(100*PCA_EVAT$sdev^2 / sum(PCA_EVAT$sdev^2), 1)

# _8.) Plot PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_evat_a <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a")

# Overall expression differences and obesity status among f1 offspring
plot_pca_evat_b <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b")

# Overall expression differences and obesity parental obesity status among f0
plot_pca_evat_c <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c")

# _9.) Combine and save overall plots ----

plot_pca_evat <- ggarrange(plot_pca_evat_a, plot_pca_evat_b, plot_pca_evat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_evat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_evat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_evat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_evat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# _10.) Get PCs for LIAT ----

PCA_LIAT <- prcomp( t (exprs(LIAT)), scale. = FALSE)
percentVar <- round(100*PCA_LIAT$sdev^2 / sum(PCA_LIAT$sdev^2), 1)

# _11.) Plot PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_liat_a <- get_pca_plot(expr_data_pca = PCA_LIAT, expr_data_raw = LIAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a")

# Overall expression differences and obesity status among f1 offspring
plot_pca_liat_b <- get_pca_plot(expr_data_pca = PCA_LIAT, expr_data_raw = LIAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b")

# Overall expression differences and obesity parental obesity status among f0
plot_pca_liat_c <- get_pca_plot(expr_data_pca = PCA_LIAT, expr_data_raw = LIAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c")

# _12.) Combine and save overall plots ----

plot_pca_liat <- ggarrange(plot_pca_liat_a, plot_pca_liat_b, plot_pca_liat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_liat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_liat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_liat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_liat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# >>>> Continue here after 02.06.2023 ----

# >>> Unfinished code section below ----

# Get numerical summaries of above PCAs ----

# use the PCA graphics and results above to get numerical summaties of hwat can be sees
# follow https://lauren-blake.github.io/Regulatory_Evol/analysis/gene_exp_corr.html
# in section "Correlations"

# use these results to inform DGE
# see DGE section to see aht models will be tested - try to match those modles here
# report those models used here and in DGE and write down the results
# possibly adjust DGE models

# >>> Unfinished code section above ----

# Re-implement analysis of array intensities  ----

# _1.) Shape and check array intensity data ----

# see https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html

# __a) Check input data formats ----

# see available expression data 
FLAT; BRAT; SCAT; LIAT; EVAT

# check full data set density
pData(FLAT) # metadata - use `ObesityLgcl` and possibly `ObeseParents`
pData(FLAT) %>% dplyr::select(ObesityLgcl, ObeseParents) %>% table()
exprs(FLAT)

# check one of four tissue data sets - brown adipose tissue 
pData(BRAT) # metadata - use `ObesityLgcl` and possibly `ObeseParents`
pData(BRAT) %>% dplyr::select(ObesityLgcl, ObeseParents) %>% table()
exprs(BRAT)

# __b) >>> Not done yet: Add array annoations to ExpreesionSet data ----

# see https://bioconductor.org/packages/release/data/annotation/html/pd.clariom.s.mouse.html


# __c) Covert expression set to data table for inspection ----

# see https://support.bioconductor.org/p/77432/
# see GSCOre manual
# https://www.bioconductor.org/packages/release/bioc/vignettes/GCSscore/inst/doc/GCSscore.pdf

# see expression set 
FLAT                   # ExpressionSet of all data
m_ints  <- exprs(FLAT) # isolate matrix of intensities
d_phen  <- pData(FLAT) # isolat data.frame of phenotypic information.

# get a data table of the exprssion set data above
#  - see https://stackoverflow.com/questions/52431288/r-data-table-how-to-go-from-tibble-to-data-table-to-tibble-back
FLAT_DT <- tibble::rownames_to_column( cbind(d_phen, t(m_ints)), var = "Sample") %>% as_tibble()  
setDT(FLAT_DT)
setkey(FLAT_DT, Sample)

# pivot data table to long - can't be done for tibble due to memory constraints
# - see https://cran.r-project.org/web/packages/data.table/vignettes/datatable-reshape.html

FLAT_DT.m1 <- melt(FLAT_DT,  id.vars = c("Sample", "Animal", "Tissue", "AnimalSex", "ObesityLgcl", "ObeseParents", "MotherDiet", "FatherDiet", "Sex", "ParentalDietMoFa", "DietGroup"), 
  variable.name = "ArrayTarget", value.name = "Intensity")


# __d) Inspect expression data raw intensities distribution  ----

# to check that equal amounts of data are available for comparison - they are not

ggplot(FLAT_DT.m1) +
  # geom_density(aes(Intensity), stat = "bin", bins = 200) +
  geom_density(aes(Intensity, colour = ObesityLgcl), stat = "bin", bins = 200) +
  ylab("Number of measurments across all other variables") +
  xlab("Intensity") +
  ggtitle("Intensities' availibilty and distribution for offsprings obesity") +
  facet_wrap(.~Tissue) + 
  theme_bw()

ggplot(FLAT_DT.m1) +
  # geom_density(aes(Intensity), stat = "bin", bins = 200) +
  geom_density(aes(Intensity, colour = ObeseParents), stat = "bin", bins = 200) +
  ylab("Number of measurments across all other variables") +
  xlab("Intensity") +
  ggtitle("Intensities' availibilty and distribution for parents obesity") +
  facet_wrap(.~Tissue) + 
  theme_bw()

# __e) Inspect expression data raw intensities' density  ----

# To check if distributions are different - hopefully they are a bit - yes perhaps in EVAT when Mother and Father are not Obese

ggplot(FLAT_DT.m1) +
  # geom_density(aes(Intensity), stat = "bin", bins = 200) +
  geom_density(aes(Intensity, colour = ObesityLgcl), stat = "density") +
  ylab("Density of measurments across all other variables") +
  xlab("Intensity") +
  ggtitle("Intensities' density for offsprings obesity") +
  facet_wrap(.~Tissue) + 
  theme_bw()

ggplot(FLAT_DT.m1) +
  # geom_density(aes(Intensity), stat = "bin", bins = 200) +
  geom_density(aes(Intensity, colour = ObeseParents), stat = "density") +
  ylab("Denisty of measurments across all other variables") +
  xlab("Intensity") +
  ggtitle("Intensities' density for parents obesity") +
  facet_wrap(.~Tissue) + 
  theme_bw()

# >>> Construction site below ----

# 2.) DGE analysis using Limma ----

# see DESeq2 tutorial at https://colauttilab.github.io/RNA-Seq_Tutorial.html
# see GCSCore tutaorial at https://www.bioconductor.org/packages/release/bioc/vignettes/GCSscore/inst/doc/GCSscore.pdf
# see Limma slides at https://s3.amazonaws.com/assets.datacamp.com/production/course_6456/slides/chapter1.pdf
# **use this!** - Limma slides at  https://kasperdanielhansen.github.io/genbioconductor/html/limma.html


# __a) Define limma models matching manuscript hypotheses ----

# ____  Test for DGE among obese and non-obese offspring ----

# No DGE  detected across all tissues or in any tissue based on offsprings' obesity status

get_dge_for_offspring_obesity(FLAT)
get_dge_for_offspring_obesity(BRAT)
get_dge_for_offspring_obesity(LIAT)
get_dge_for_offspring_obesity(SCAT)
get_dge_for_offspring_obesity(EVAT)

# __b) Define limma models partially matching manuscript hypotheses (fall-back) ----

# ____  Test for DGE among offspring based on parental obesity ----

# Defining and applying contrasts: One of "MotherFatherNotObese", "FatherObese", "MotherFatherObese", or "MotherObese" 
#  against all remaining three levels.

FLAT_TopTableList <- get_dge_for_parent_obesity(FLAT)
BRAT_TopTableList <- get_dge_for_parent_obesity(BRAT)
LIAT_TopTableList <- get_dge_for_parent_obesity(LIAT)
EVAT_TopTableList <- get_dge_for_parent_obesity(EVAT)

# >>> Construction site above ---- 

# Experimental: DGE-analysis using GAMs - Investigate overall tissue specific expression differences based on obesity variables ----

# from inspection 
# - check EVAT of MotherFatherNotObese across all genes
# - check LIAT and MotherFatherObese across all genes
# - also check Obesity Lgl

# _1.) Set up parallel computing ----

# parallel::makeForkCluster(nnodes = getOption("mc.cores", 6L))

# _2.) Select modelling data ----

# select model data size for mdoel testing and code develpment - decrese number is code is getting too slow
ArrayTargets <- sample(FLAT_DT.m1[["ArrayTarget"]], 500, replace=FALSE)
model_data <- subset(FLAT_DT.m1, ArrayTarget %in% ArrayTargets) # subes to 100 of 20000 genes to speed up compuaterion
unique(model_data[["Tissue"]])

# _3.) Build model ----

# mod_0 <- gam(ObesityLgcl ~ s(Intensity, k=6, bs="tp", m=2) + Tissue + s(ArrayTarget, k=6, bs="tp"), data = model_data, method = "ML", family = binomial(link = "logit"))

# _4.) Save model ----

# _5.) Check model ----

# summary(mod_0)
# gam.check(mod_0)
# gratia::appraise(mod_0)

# _6.) Inspect model predictions  ----

# mod_0_predictions  <- get_model_prdictions_with_mgcv(mod_0, model_data)
# 
# ggplot(data = mod_0_predictions, aes(x = ObesityLgcl, y = Intensity, colour = Tissue)) +
#   geom_point(aes(y = Intensity,), alpha = 0.5)
#   geom_line(aes(y = fit, group = AnimalId), alpha = 0.5, linewidth = 0.2) +
#   facet_wrap(ObeseParents ~ AnimalSex, ncol = 2) + 
#   theme_bw() + 
#   labs(title = "Offsprings body weight by sex and parents obesity statuts", 
#        subtitle = paste("R model formula: ", as.character(paste(deparse(formula(mod_4), width.cutoff = 500), collapse=""))),
#        x="age [d]", y = "body weight [g]")
# 

# Snapshot environment) ----
sessionInfo()
save.image(file = here("scripts", "050_r_array_analysis_ah.RData"))
renv::snapshot()


# >>> Construction and old code code below: Volcano plot  ----

# AH code below - continue here after 31.05.2023 - also see in section 3 data read in - define Obesity variables from AH's diet variables

pVal <- resultTable$adj.P.Val   #adj. p-Value für cutoffs verwenden
CO1 <- 0.05                     #adj. p-Value cutoff
CO2 <- 0.5                      #FC cutoff
#vermutlich kannst du stringenter sein und als adj. P 0.01 nehmen, da wir sehr viele Hits haben
#(vielleicht mehr eingrenzen, um die Analysen zu vereinfachen)
#FC cutoff kannst du evt. auf 0.75 hochnehmen (bei 1 bleibt nicht mehr so viel übrig), aber p-Wert ist besser

rm(list=ls())
gc()


# use DGEL[[n]]]

resultTable <- eWAT_CD_HFD_VS_CD_CD #such dir die aus, die du plotten möchtest
para2 <- "eWAT_CD_HFD_VS_CD_CD"      #variable für path zum speichern des plots

resultTable$Col <- ifelse(pVal > CO1, "928a97", "grey40")
resultTable$Col <- ifelse(pVal < CO1 & resultTable$logFC > CO2, "#d72323", resultTable$Col)
resultTable$Col <- ifelse(pVal > CO1 & resultTable$logFC > CO2, "tomato", resultTable$Col)
resultTable$Col <- ifelse(pVal < CO1 & resultTable$logFC < -CO2, "blue", resultTable$Col)
resultTable$Col <- ifelse(pVal > CO1 & resultTable$logFC < -CO2, "dodgerblue", resultTable$Col)
  
#filter top genes 
sigGenes <- c(grep("#d72323",resultTable$Col),grep("blue",resultTable$Col))
topLFS <- resultTable[sigGenes,]
nosigGenes <- grep("dodgerblue",topLFS$Col)
if(length(nosigGenes) != 0){topLFS <- topLFS[-nosigGenes,]}
topLFS <- topLFS[ order( topLFS$adj.P.Val, decreasing=FALSE), ]    
topLFS <- head(topLFS, n=20)
resultTable$LABEL <- ifelse(match(resultTable$SYMBOL,topLFS$SYMBOL),as.character(resultTable$SYMBOL), NA)
  
#plot it
ggplot(data=resultTable, aes(x=logFC, y=-log10(pVal))) +
    geom_point(colour=resultTable$Col,size=2, alpha=0.7)+
    xlab("log2 fold change") + 
    ylab("-log10 adj. p-value") +     #change here
    geom_text_repel(aes(label=as.character(LABEL)),
                    direction="both",
                    max.overlaps = Inf,
                    min.segment.length = 0)+
    geom_hline(yintercept= -log10(CO1),  linetype="dotted", color="blue") + 
    geom_vline(xintercept= CO2, linetype="dotted", color="blue") + 
    geom_vline(xintercept= -CO2,  linetype="dotted",color="blue")+
    theme_classic()+
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15))
ggsave(paste0("Results_05.23/",para2,"_Volcano.pdf"),width = 6, height = 6)

# >>> Construction and old code code below: Upset plot  ----

 
# Upset plots ############################# 
#hier am Beispiel für Liver
para <- "Liver"

sub1 <- subset(Liv_CD_HFD_VS_CD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
colnames(sub1) <- c("RN","CD, HFD vs. CD, CD")
sub2 <- subset(Liv_HFD_CD_VS_CD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
colnames(sub2) <- c("RN","HFD, CD vs. CD, CD")
sub3 <- subset(Liv_HFD_HFD_VS_CD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
colnames(sub3) <- c("RN","HFD, HFD vs. CD, CD")
sub4 <- subset(Liv_HFD_CD_VS_CD_HFD, Col=="#d72323" | Col=="blue")[,c(1,8)]
colnames(sub4) <- c("RN","HFD, CD vs. CD, HFD")
sub5 <- subset(Liv_HFD_HFD_VS_CD_HFD, Col=="#d72323" | Col=="blue")[,c(1,8)]
colnames(sub5) <- c("RN","HFD, HFD vs. CD, HFD")
sub6 <- subset(Liv_HFD_HFD_VS_HFD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
colnames(sub6) <- c("RN","HFD, HFD vs. HFD, CD")

mTAB <- merge(merge(merge(merge(merge(
  sub1,
  sub2, all = TRUE,by="RN"),
  sub3, all = TRUE,by="RN"),
  sub4, all = TRUE,by="RN"),
  sub5, all = TRUE,by="RN"),
  sub6, all = TRUE,by="RN")

mTAB$`CD, HFD vs. CD, CD` <- ifelse(is.na(mTAB$`CD, HFD vs. CD, CD`),  0,1)
mTAB$`HFD, CD vs. CD, CD` <- ifelse(is.na(mTAB$`HFD, CD vs. CD, CD`),  0,1)
mTAB$`HFD, HFD vs. CD, CD` <- ifelse(is.na(mTAB$`HFD, HFD vs. CD, CD`),  0,1)
mTAB$`HFD, CD vs. CD, HFD` <- ifelse(is.na(mTAB$`HFD, CD vs. CD, HFD`),  0,1)
mTAB$`HFD, HFD vs. CD, HFD` <- ifelse(is.na(mTAB$`HFD, HFD vs. CD, HFD`),  0,1)
mTAB$`HFD, HFD vs. HFD, CD` <- ifelse(is.na(mTAB$`HFD, HFD vs. HFD, CD`),  0,1)

pdf(paste0("Results_05.23/",para,"_Upset.pdf"),width = 9,height = 5) 
  UpSetR::upset(mTAB,nsets = 6, decreasing=TRUE,order.by = "freq",group.by="degree",
              sets.bar.color=c("maroon","blue","orange","purple","steelblue","red"))
dev.off()

#für GO und KEGG nimmst du dann die DGE Tabellen und filterst nach den adj. P-Value und FC cutoffs


