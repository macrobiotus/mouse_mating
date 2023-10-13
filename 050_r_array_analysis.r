#' ---
#' title: "Mice Mating Study"
#' subtitle: "Re-analysis of array data"
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

#' ##  Packages

# _2.) Packages ----

library(here)
library(renv)

library(magrittr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggpubr)

library(factoextra)

library(pd.clariom.s.mouse) # get array sequences
library(Biobase)
library(BiocGenerics)
library(DESeq2)
library(affycoretools) # Functions useful for those doing repetitive analyses with Affymetrix GeneChips
library(org.Mm.eg.db)

library(swamp)

library(pheatmap)
library(EnhancedVolcano)

# Unused - for 3D PCA
# library(plotly)
# library(htmlwidgets)

# library(mgcv)
# library(parallel, lib.loc = "/Users/paul/Library/Caches/org.R-project.R/R/renv/sandbox/R-4.2/aarch64-apple-darwin20/84ba8b13")
# library(RVAideMemoire)
# library("GCSscore")
# library("oligoClasses") # to annotate array target
# library("AnnoProbe")  # Annotate the Gene Symbols for Probes in Expression Array
# # for GO analysis
# library("limma")
# library(gplots)
# library(RColorBrewer)
# library(KernSmooth)
# library(data.table)
# library(FactoMineR)
# library(UpSetR)
# library(usethis)

#' ## Increase memory if needed

# _3.) Increase meomeory if needed ----

# usethis::edit_r_environ()

#' ##  Functions

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
get_pca_plot = function(expr_data_pca, expr_data_raw, variable, legend_title, plot_title, percent_var) {
  require("factoextra")
  require("ggplot2")
  
  # correct legend labels (inset 9-Oct-2023)
  if(variable == "Tissue"){
  correct_lables <- case_when(pData(expr_data_raw)[[variable]] == "BRAT" ~ "BAT",
            pData(expr_data_raw)[[variable]] == "EVAT" ~ "EVAT",
            pData(expr_data_raw)[[variable]] == "LIAT" ~ "L",
            pData(expr_data_raw)[[variable]] == "SCAT" ~ "SCAT",
            TRUE ~ pData(expr_data_raw)[[variable]])
  } else if (variable != "Tissue") {
    correct_lables <- pData(expr_data_raw)[[variable]]
    
  }
  
  pca_plot <- fviz_pca_ind(
    expr_data_pca,
    label = "none",
    habillage = factor(correct_lables),
    pointsize = 2,
    palette = c("firebrick3", "purple", "steelblue3", "gold3"),
    legend.title = legend_title,
    invisible = "quali",
    addEllipses = FALSE,
    ellipse.level=0.95,
    title = plot_title
  ) +
    labs(
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
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

# Convenience function to define DGE among offsprings tissue types
get_dge_for_individal_tissues =  function(ExpSet){
  
  message("Convenience function to test DGE for each tissue compred to all other tissues, only makes sense for a data set with different tissue type. Analysis of factor \"Tissue\" is hard-coded.")
  
  # Building initial models
  message(paste0("\nUsing data set \"", deparse(substitute(ExpSet))), "\".")   
  
  
  message("\nConsider whether the sample size is appropriate for DGE, you should have at least six samples per group:")
  
  pData(ExpSet) %>% dplyr::select(ObeseParents) %>% table() %>% print()
  
  # Building initial models
  message("Building initial model.")   
  
  design_tissue_types <- model.matrix(~ ExpSet[["Tissue"]] - 1)
  colnames(design_tissue_types) <-c("BRAT", "EVAT", "LIAT", "SCAT") 
  fit_tissue_types <- lmFit(ExpSet, design_tissue_types)
  
  # Setting contrasts - likley not all will give results 
  contrast_names <- c(
    "BRAT vs EVAT & LIAT & SCAT",
    "EVAT vs BRAT & LIAT & SCAT",
    "LIAT vs EVAT & BRAT & SCAT",
    "SCAT vs EVAT & BRAT & LIAT")
  
  message("\nDefining contrasts: \"", paste0(contrast_names, collapse = "\", \""), "\".") 
  
  contrast_list <- vector(mode = "list", length = length(contrast_names))
  
  contrast_list[[1]]  <- makeContrasts("BRAT vs EVAT & LIAT & SCAT" =  BRAT - (EVAT + LIAT + SCAT)/3, levels = design_tissue_types)
  contrast_list[[2]]  <- makeContrasts("EVAT vs BRAT & LIAT & SCAT" =  EVAT - (BRAT + LIAT + SCAT)/3, levels = design_tissue_types)
  contrast_list[[3]]  <- makeContrasts("LIAT vs EVAT & BRAT & SCAT" =  LIAT - (EVAT + BRAT + SCAT)/3, levels = design_tissue_types)
  contrast_list[[4]]  <- makeContrasts("SCAT vs EVAT & BRAT & LIAT" =  SCAT - (EVAT + BRAT + LIAT)/3, levels = design_tissue_types)
  
  message("\nApplying contrasts.") 
  
  # Create list and fill it with contrats applied to the initial model
  contrast_model_list <- vector(mode = "list", length = length(contrast_names))
  contrast_model_list <- lapply(contrast_list, function (x) contrasts.fit(fit_tissue_types, x))
  contrast_model_list_eb <- lapply(contrast_model_list, function (x) eBayes(x))
  contrast_model_list_tp <- lapply(contrast_model_list_eb, function (x) topTable(x, p.value = 0.50, number = Inf))
  
  message("\nReturning results list - check list names for top table identification.") 
  
  # Setting names for identifying contrast among results
  names(contrast_model_list_tp) <- contrast_names
  return(contrast_model_list_tp)
}

# Convenience function to test offsprings' obesity status in all tissue samples
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
  
  return(topTable(fit1c, p.value = 0.50, number = Inf))
  
}

# Convenience function to test parents' obesity status in all tissue samples but LIAT
get_dge_for_parent_obesity = function(ExpSet){
  
  message("Convenience function to test parents' obesity status in all tissue samples - analysis of factor \"ObeseParents\" is hard-coded.")
  
  # Building initial models
  message(paste0("\nUsing data set \"", deparse(substitute(ExpSet))), "\".")   
  
  message("\nConsider whether the sample size is appropriate for DGE, you should have at least six samples per group:")
  
  pData(ExpSet) %>% dplyr::select(ObeseParents) %>% table() %>% print()
  
  # Building initial models
  message("Building initial model.")   
  
  design_parent_obese <- model.matrix(~ ExpSet[["ObeseParents"]] - 1)
  colnames(design_parent_obese) <-c("MotherFatherNotObese", "FatherObese", "MotherFatherObese", "MotherObese") 
  fit_parent_obese <- lmFit(ExpSet, design_parent_obese)
  
  # Setting contrasts - likely not all will give results 
  
  contrast_names <- c(
    "MotherFatherObese vs MotherObese",
    "MotherFatherObese vs FatherObese", 
    "MotherFatherObese vs FatherObese & MotherObese",
    "MotherFatherObese vs MotherFatherNotObese",
    "MotherFatherObese vs MotherFatherNotObese & MotherObese",
    "MotherFatherObese vs MotherFatherNotObese & FatherObese",
    "MotherFatherObese vs MotherFatherNotObese & FatherObese & MotherObese",
    "MotherFatherNotObese vs MotherObese", 
    "MotherFatherNotObese vs FatherObese", 
    "MotherFatherNotObese vs FatherObese & MotherObese", 
    "MotherFatherNotObese vs MotherFatherObese & MotherObese", 
    "MotherFatherNotObese vs MotherFatherObese & FatherObese", 
    "MotherFatherNotObese vs MotherFatherObese & FatherObese & MotherObese",
    "FatherObese vs MotherObese",
    "MotherObese vs MotherFatherNotObese & FatherObese",
    "FatherObese vs MotherFatherObese & MotherObese",
    "FatherObese vs MotherFatherNotObese & MotherFatherObese", 
    "MotherObese vs MotherFatherNotObese & MotherFatherObese")
  
  message("\nDefining contrasts: \"", paste0(contrast_names, collapse = "\", \""), "\".") 
  
  contrast_list <- vector(mode = "list", length = length(contrast_names))
  
  contrast_list[[1]]  <- makeContrasts("MotherFatherObese vs MotherObese" =  MotherFatherObese - MotherObese, levels = design_parent_obese)
  contrast_list[[2]]  <- makeContrasts("MotherFatherObese vs FatherObese" =  MotherFatherObese - FatherObese, levels = design_parent_obese)
  contrast_list[[3]]  <- makeContrasts("MotherFatherObese vs FatherObese & MotherObese" =  MotherFatherObese - (FatherObese + MotherObese)/2 , levels = design_parent_obese)
  contrast_list[[4]]  <- makeContrasts("MotherFatherObese vs MotherFatherNotObese" =  MotherFatherObese - MotherFatherNotObese, levels = design_parent_obese)
  contrast_list[[5]]  <- makeContrasts("MotherFatherObese vs MotherFatherNotObese & MotherObese" =  MotherFatherObese - (MotherFatherNotObese + MotherObese)/2, levels = design_parent_obese)
  contrast_list[[6]]  <- makeContrasts("MotherFatherObese vs MotherFatherNotObese & FatherObese" =  MotherFatherObese - (MotherFatherNotObese + FatherObese)/2, levels = design_parent_obese)
  contrast_list[[7]]  <- makeContrasts("MotherFatherObese vs MotherFatherNotObese & FatherObese & MotherObese" =  MotherFatherObese - (MotherFatherNotObese + FatherObese + MotherObese)/3, levels = design_parent_obese)
  contrast_list[[8]]  <- makeContrasts("MotherFatherNotObese vs MotherObese" =  MotherFatherNotObese - MotherObese, levels = design_parent_obese)
  contrast_list[[9]]  <- makeContrasts("MotherFatherNotObese vs FatherObese" =  MotherFatherNotObese - FatherObese, levels = design_parent_obese)
  contrast_list[[10]] <- makeContrasts("MotherFatherNotObese vs FatherObese & MotherObese" = MotherFatherNotObese - (FatherObese + MotherObese)/2, levels = design_parent_obese)
  contrast_list[[11]] <- makeContrasts("MotherFatherNotObese vs MotherFatherObese & MotherObese" = MotherFatherNotObese - (MotherFatherObese + MotherObese)/2, levels = design_parent_obese)
  contrast_list[[12]] <- makeContrasts("MotherFatherNotObese vs MotherFatherObese & FatherObese" = MotherFatherNotObese - (MotherFatherObese + FatherObese)/2, levels = design_parent_obese)
  contrast_list[[13]] <- makeContrasts("MotherFatherNotObese vs MotherFatherObese & FatherObese & MotherObese" = MotherFatherNotObese - (MotherFatherObese + FatherObese + MotherObese)/3, levels = design_parent_obese)
  contrast_list[[14]] <- makeContrasts("FatherObese vs MotherObese" = FatherObese - MotherObese, levels = design_parent_obese)
  contrast_list[[15]] <- makeContrasts("MotherObese vs MotherFatherNotObese & FatherObese" = MotherObese - (MotherFatherNotObese + FatherObese)/2, levels = design_parent_obese)
  contrast_list[[16]] <- makeContrasts("FatherObese vs MotherFatherObese & MotherObese" = FatherObese - (MotherFatherObese + MotherObese)/2, levels = design_parent_obese)
  contrast_list[[17]] <- makeContrasts("FatherObese vs MotherFatherNotObese & MotherFatherObese" = FatherObese - (MotherFatherNotObese + MotherFatherObese)/2, levels = design_parent_obese)
  contrast_list[[18]] <- makeContrasts("MotherObese vs MotherFatherNotObese & MotherFatherObese" = MotherObese - (MotherFatherNotObese + MotherFatherObese)/2, levels = design_parent_obese)
  
  message("\nApplying contrasts.") 
  
  # Create list and fill it with contrats applied to the initial model
  contrast_model_list <- vector(mode = "list", length = length(contrast_names))
  contrast_model_list <- lapply(contrast_list, function (x) contrasts.fit(fit_parent_obese, x))
  contrast_model_list_eb <- lapply(contrast_model_list, function (x) eBayes(x))
  contrast_model_list_tp <- lapply(contrast_model_list_eb, function (x) topTable(x, p.value = 0.50, number = Inf))
  
  message("\nReturning results list - check list names for top table identification.") 
  
  # Setting names for identifying contrast among results
  names(contrast_model_list_tp) <- contrast_names
  return(contrast_model_list_tp)
}

# Convenience function to test parents' obesity status in LIAT samples
get_some_dge_for_parent_obesity = function(ExpSet){
  
  message("Reduced contrast set for LIAT analysis") 
  message("Convenience function to test parents' obesity status in all tissue samples - analysis of factor \"ObeseParents\" is hard-coded.")
  
  # Building initial models
  message(paste0("\nUsing data set \"", deparse(substitute(ExpSet))), "\".")   
  
  message("\nConsider whether the sample size is appropriate for DGE, you should have at least six samples per group:")
  
  pData(ExpSet) %>% dplyr::select(ObeseParents) %>% table() %>% print()
  
  # Building initial models
  message("Building initial model.")   
  
  design_parent_obese <- model.matrix(~ ExpSet[["ObeseParents"]] - 1)
  colnames(design_parent_obese) <-c("MotherFatherNotObese", "FatherObese", "MotherFatherObese", "MotherObese") 
  fit_parent_obese <- lmFit(ExpSet, design_parent_obese)
  
  # Setting contrasts - likely not all will give results 
  
  contrast_names <- c(
    "MotherFatherObese vs FatherObese", 
    "MotherFatherObese vs MotherFatherNotObese",
    "MotherFatherObese vs MotherFatherNotObese & FatherObese",
    "MotherFatherNotObese vs FatherObese", 
    "MotherFatherNotObese vs MotherFatherObese & FatherObese", 
    "FatherObese vs  MotherFatherNotObese & MotherFatherObese")
  
  message("\nDefining contrasts: \"", paste0(contrast_names, collapse = "\", \""), "\".") 
  
  contrast_list <- vector(mode = "list", length = length(contrast_names))
  
  contrast_list[[1]]  <- makeContrasts("MotherFatherObese vs FatherObese" =  MotherFatherObese - FatherObese, levels = design_parent_obese)
  contrast_list[[2]]  <- makeContrasts("MotherFatherObese vs MotherFatherNotObese" =  MotherFatherObese - MotherFatherNotObese, levels = design_parent_obese)
  contrast_list[[3]]  <- makeContrasts("MotherFatherObese vs MotherFatherNotObese & FatherObese" =  MotherFatherObese - (MotherFatherNotObese + FatherObese)/2, levels = design_parent_obese)
  contrast_list[[4]]  <- makeContrasts("MotherFatherNotObese vs FatherObese" =  MotherFatherNotObese - FatherObese, levels = design_parent_obese)
  contrast_list[[5]] <- makeContrasts("MotherFatherNotObese vs MotherFatherObese & FatherObese" = MotherFatherNotObese - (MotherFatherObese + FatherObese)/2, levels = design_parent_obese)
  contrast_list[[6]] <- makeContrasts("FatherObese vs MotherFatherNotObese & MotherFatherObese" = FatherObese - (MotherFatherNotObese + MotherFatherObese)/2, levels = design_parent_obese)
  
  message("\nApplying contrasts.") 
  
  # Create list and fill it with contrats applied to the initial model
  contrast_model_list <- vector(mode = "list", length = length(contrast_names))
  contrast_model_list <- lapply(contrast_list, function (x) contrasts.fit(fit_parent_obese, x))
  contrast_model_list_eb <- lapply(contrast_model_list, function (x) eBayes(x))
  contrast_model_list_tp <- lapply(contrast_model_list_eb, function (x) topTable(x, p.value = 1.0, number = Inf))
  
  message("\nReturning results list - check list names for top table identification.") 
  
  # Setting names for identifying contrast among results
  names(contrast_model_list_tp) <- contrast_names
  return(contrast_model_list_tp)
}

# Get principal components of expression data
get_principal_components = function(ExpSet){
  
  # Building initial models
  message(paste0("\nGetting PCs of data set \"", deparse(substitute(ExpSet))), "\".")
  
  return( prcomp( t( exprs(ExpSet)), scale = TRUE, center = TRUE))
  
}

# Get percentage of variation of the PCs, for plotting
get_percent_variation = function(PcExpSet) {
  
  # Building initial models
  message(paste0("\nGetting percent variation of PC set \"", deparse(substitute(PcExpSet))), "\".")
  
  return(round(100*PcExpSet$sdev^2 / sum(PcExpSet$sdev^2), 1) )
  
}

# Isolate First 5 Principal Components
get_isolated_pcs = function(pca_ob){
  
  # Show which data is being processed
  message(paste0("\nGetting first 5 PCs of data set \"", deparse(substitute(pca_ob))), "\".")
  
  return(data.frame("PC1" = pca_ob[["x"]][ , "PC1"], "PC2" = pca_ob[["x"]][ , "PC2"],  "PC3" = pca_ob[["x"]][ , "PC3"], "PC4" = pca_ob[["x"]][ , "PC4"], "PC5" = pca_ob[["x"]][ , "PC5"]))
  
}

# Get data frames to model variation of first PC against chosen factors
get_model_data = function(ExprSet, PCsSet) {
  
  # require
  require(dplyr)
  require(tidyr)
  
  # Show which data is being processed
  message(paste0("For modeling uniting expression metadata  \"", deparse(substitute(ExprSet)), "\" and PC data \"", deparse(substitute(PCsSet)), "\"."))
  
  # Format input data as tiblles for joining
  PCsSet <- as_tibble(PCsSet, rownames = "Sample")
  ExpMData <- pData(ExprSet) %>% as_tibble %>% unite("Sample", "Animal", "Tissue", sep = "_", remove = FALSE)
  
  # Unite both data sets
  model_data <- left_join(PCsSet, ExpMData)
  
  return(model_data)  
  
}

# Get Volcano plots
get_one_volcanoplot <- function(TopTableListItem, TopTableListItemName){
  
  # transform input item for plotting function 
  top_tibble <- as_tibble(TopTableListItem)
  
  # diagnostic
  message(paste0("Creating plot for contrast: \"", TopTableListItemName, "\"", sep = ""))
  
  # get plot
  evplot <-
    EnhancedVolcano(
      top_tibble,
      pCutoff = 0.05,
      FCcutoff = 1,
      x = "logFC",
      y = "P.Value",
      pCutoffCol = "adj.P.Val",
      lab = top_tibble[["SYMBOL"]],
      title = TopTableListItemName,
      subtitle = NULL
    )
  
  # scales y axis
  evplot <-  evplot + ylim(0, 5)
  
  # return plot 
  return(evplot)
  
}

# Lookup Entrez IDs
get_entrez_ids = function(TopTable) {
  
  require("org.Mm.eg.db")
  require("clusterProfiler")
  require("biomaRt")
  
  message(paste0("Looking up Entrez IDs for data set."))  
  
  # create target object
  TopTableAppended <- TopTable
  
  # lookup Entrez ID from RefSeq ID and store both columns for separtaion
  translations <- bitr(geneID = TopTable[["ID"]], fromType = "REFSEQ", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = FALSE)
  
  # Move second column to traget obeject 
  TopTableAppended[["ENTREZ"]] <- translations[ , 2]
          
  return(TopTableAppended)
          
}

get_one_kegg_dotplot <- function(TopTableListItem, TopTableListItemName){
  
  # see https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
  
  require("clusterProfiler")
  require("enrichplot")
  
  # diagnostic
  message(paste0("Creating KEGG dot-plot for data set: \"", TopTableListItemName, "\".", sep = ""))
  
  # lookup kegg pathways
  kegg_result <- enrichKEGG(gene = TopTableListItem$ENTREZ, keyType = 'ncbi-geneid',  organism = 'mmu', minGSSize = 5)  # ncib-proteinid is not supported for mmu ...
  
  # check resulting object
  print(kegg_result)
  
  # get dipslay item
  kegg_plot <- enrichplot::dotplot(kegg_result, title =  paste0("KEGG pathways of data set: \"", TopTableListItemName, "\"", sep = ""), showCategory=15)
  
  # return plot 
  return(kegg_plot)
  
}

get_one_go_plot <- function(TopTableListItem, TopTableListItemName){
  
  require("clusterProfiler")
  require("enrichplot")
  
  # diagnostic
  message(paste0("Creating GO plot for data set: \"", TopTableListItemName, "\".", sep = ""))
  
  # lookup go pathways
  go_result <- enrichGO(gene = TopTableListItem$ENTREZ, keyType = "ENTREZID",  OrgDb = "org.Mm.eg.db", ont = "all")  # ncib-proteinid is not supported for mmu ...
  
  # get dipslay item
  go_plot <- enrichplot::dotplot(go_result, split="ONTOLOGY", title =  paste0("GO terms of data set: \"", TopTableListItemName, "\"", sep = ""), showCategory=15) + facet_grid(ONTOLOGY ~ ., scale="free")
  
  # return plot 
  return(go_plot)
  
}

save_kegg_plots <- function(ggplot_list_item, ggplot_list_item_name){
  
  require(ggplot2)
  
  filename <- gsub("[^[:alnum:][:space:]]","", ggplot_list_item_name)
  filename <- gsub("\\s+","_", filename)
  filename <- paste0("050_r_array_analysis__plot_kegg_", filename, ".pdf", collapse = "_")
  
  message(filename)
  try(ggsave(plot = ggplot_list_item, path = here("plots"), filename, scale = 0.75))
  
}

save_go_plots <- function(ggplot_list_item, ggplot_list_item_name){
  
  require(ggplot2)
  
  filename <- gsub("[^[:alnum:][:space:]]","", ggplot_list_item_name)
  filename <- gsub("\\s+","_", filename)
  filename <- paste0("050_r_array_analysis__plot_go_", filename, ".pdf", collapse = "_")
  
  message(filename)
  try(ggsave(plot = ggplot_list_item, path = here("plots"), filename, scale = 0.75, height = 800, width = 400, units = c("mm")))
  
}

#' ##  Color code for plotting

# _5.) Color code for plotting ----

hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(n=11, name="RdBu"))(50))

#' # Load and shape data

# Load and shape data ----

#' ## Full data of all 4 tissues

# _1.) Full data of all 4 tissues ----

#' ### Loading normalized data

# __a.) Loading normalized data ----

# Data are Clariom S mouse arrays
# I removed strong outliers: A285_liver clustert zu bAT, A339_liver weit weg vom Rest

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/allTissues_normData.RData") # only if you are interested to look into the normalized data

#' ### Annotate data

# __b.) Annotate data ----

normData  <- annotateEset(normData, pd.clariom.s.mouse, type = "probeset")

#' ### Copy, store, and discard data

# __c.) Copy, store, and discard data ----

# copy to stick to manuscript naming conventions
FLAT <- normData; rm(normData)

saveRDS(FLAT, file = here("rds_storage", "050_r_array_analysis__normalized_data.rds"))

#' ## Data normalized individually for each tissue

# _2.) Data normalized individually for each tissue ----

# Ich habe die normalisierung für jedes Gewebe getrennt für die DGE Analysen gemacht, 
# da dies genauer ist (Gewebe zu weit auseinander in PCA)

#' ### Loading normalized data

# __a.) Loading normalized data ----

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/normData4DGE.RData") #für jedes Gewebe die normalisierten Daten

#' ### Annotate data

# __b.) Annotate data ----

bAT_normData    <- annotateEset(bAT_normData, pd.clariom.s.mouse, type = "probeset")
ingWAT_normData <- annotateEset(ingWAT_normData, pd.clariom.s.mouse, type = "probeset")
Liv_normData    <- annotateEset(Liv_normData, pd.clariom.s.mouse, type = "probeset")
bAT_normData    <- annotateEset(bAT_normData, pd.clariom.s.mouse, type = "probeset")
eWAT_normData   <- annotateEset(eWAT_normData, pd.clariom.s.mouse, type = "probeset")

#' ### Copy, store, and discard data

# __c.) Copy, store, and discard data ----

# copy to stick to manuscript naming conventions - corrcted as per AH 25.05.2023
BRAT <- bAT_normData; rm(bAT_normData)        # brown adipose tissue
SCAT <- ingWAT_normData; rm(ingWAT_normData)  # subcutabneous adipose tissue 
LIAT <- Liv_normData; rm(Liv_normData)        # liver adipose tissue
EVAT <- eWAT_normData; rm(eWAT_normData)      # visceral adipose tissue

saveRDS(BRAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_bat.rds"))
saveRDS(SCAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_ewat.rds"))
saveRDS(LIAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_liv.rds" ))
saveRDS(EVAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_ingwat.rds"))

#' ## Loading metadata from modelling (obesity variables)

# _3.) Loading metadata from modelling (obesity variables) ----

mice_f1_modeled_data_with_rna_seq_data <- readRDS(file = here("rds_storage", "040_r_h3__mice_f1_modeled_data_with_rna_seq_data.rds"))

#' ## Adjust variable names and inspect data

# _4.) Adjust variable names and inspect data ----

FLAT <- adjust_array_data(FLAT, mice_f1_modeled_data_with_rna_seq_data)
BRAT <- adjust_array_data(BRAT, mice_f1_modeled_data_with_rna_seq_data)
SCAT <- adjust_array_data(SCAT, mice_f1_modeled_data_with_rna_seq_data)
LIAT <- adjust_array_data(LIAT, mice_f1_modeled_data_with_rna_seq_data)
EVAT <- adjust_array_data(EVAT, mice_f1_modeled_data_with_rna_seq_data)

# _5.) Save adjusted data ----

saveRDS(BRAT, file = here("rds_storage", "050_r_array_analysis__BRAT.rds"))
saveRDS(SCAT, file = here("rds_storage", "050_r_array_analysis__SCAT.rds"))
saveRDS(LIAT, file = here("rds_storage", "050_r_array_analysis__LIAT.rds"))
saveRDS(EVAT, file = here("rds_storage", "050_r_array_analysis__EVAT.rds"))

# check, if you like, using `pData(BRAT)` and `exprs(BRAT)` - check available metadata - LIAT does not have a lot

lapply(c(FLAT, BRAT, SCAT, LIAT, EVAT), pData)  

#' ## Loading AHs dietary DGE analysis results

# _6.) Loading AHs dietary DGE analysis results ----

#' ### Load data

# __a) Load data ----

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/DGELists.RData") #lädt alle DGE Tabellen

#' ### Define list with dietary variables ----

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

#' ### Clean environment and save data

# __d) Clean environment and save data ---- 

rm(bAT_CD_HFD_VS_CD_CD, bAT_HFD_CD_VS_CD_CD, bAT_HFD_CD_VS_CD_HFD, bAT_HFD_HFD_VS_CD_CD, bAT_HFD_HFD_VS_CD_HFD,    
   bAT_HFD_HFD_VS_HFD_CD, eWAT_CD_HFD_VS_CD_CD, eWAT_HFD_CD_VS_CD_CD, eWAT_HFD_CD_VS_CD_HFD, eWAT_HFD_HFD_VS_CD_CD, 
   eWAT_HFD_HFD_VS_CD_HFD, eWAT_HFD_HFD_VS_HFD_CD, ingWAT_CD_HFD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_HFD,
   ingWAT_HFD_HFD_VS_CD_CD, ingWAT_HFD_HFD_VS_CD_HFD, ingWAT_HFD_HFD_VS_HFD_CD, Liv_CD_HFD_VS_CD_CD, Liv_HFD_CD_VS_CD_CD,
   Liv_HFD_CD_VS_CD_HFD, Liv_HFD_HFD_VS_CD_CD, Liv_HFD_HFD_VS_CD_HFD, Liv_HFD_HFD_VS_HFD_CD)

saveRDS(DGEL_diet, file = here("rds_storage", "050_r_array_analysis__dge_lists_by_diet.rds"))

#' # Calculate Principal components

# Calculate Principal components ----

#' ## Get PCs for all subsequent analyses

# _1.) Get PCs for all subsequent analyses ----

PCA_FLAT <- get_principal_components(FLAT)
PCA_BRAT <- get_principal_components(BRAT)
PCA_SCAT <- get_principal_components(SCAT)
PCA_LIAT <- get_principal_components(LIAT)
PCA_EVAT <- get_principal_components(EVAT)

#' ## Get PC loadings for plots

# _2.) Get PC loadings for plots ----

FLAT_PV <- get_percent_variation(PCA_FLAT)
BRAT_PV <- get_percent_variation(PCA_BRAT)
SCAT_PV <- get_percent_variation(PCA_SCAT)
LIAT_PV <- get_percent_variation(PCA_LIAT)
EVAT_PV <- get_percent_variation(PCA_EVAT)

# #' Plot Principal Component Analyses

# Plot Principal Component Analyses ----

#' ## Plot FLAT PCs in 2D for several variables types

# _1.) Plot FLAT PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_flat_a <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a", percent_var = FLAT_PV)

# Overall expression differences and obesity status among f1 offspring
plot_pca_flat_b <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b", percent_var = FLAT_PV)

# Overall expression differences and obesity parental obesity status among f0
plot_pca_flat_c <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c", percent_var = FLAT_PV)

# ___ Combine and save plots ----

plot_pca_flat <- ggarrange(plot_pca_flat_a, plot_pca_flat_b, plot_pca_flat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_flat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_flat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_flat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_flat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# ___ Plot PCs in 3D (only for tissue types so far)  ----

# dataGG <- data.frame(PC1 = PCA_FLAT$x[,1], PC2 = PCA_FLAT$x[,2],PC3 = PCA_FLAT$x[,3])
# 
# plot_pca_sptial <- plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3) %>%
#   add_markers(color = ~pData(FLAT)$Tissue,colors = c("firebrick3","purple","steelblue3","gold3"),
#               symbol = ~pData(FLAT)$Tissue, symbols = c(19,19,19,19))
# 
# saveWidget(plot_pca_sptial, file = paste0(here("plots"),"/", "050_r_array_analysis__plot_pca_flat.html"), selfcontained = T, libdir = "lib")

#' ## Plot BRAT PCs in 2D for several variables types

# _2.) Plot BRAT PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_brat_a <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a", percent_var = BRAT_PV)

# Overall expression differences and obesity status among f1 offspring
plot_pca_brat_b <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b", percent_var = BRAT_PV)

# Overall expression differences and obesity parental obesity status among f0
plot_pca_brat_c <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c", percent_var = BRAT_PV)

# ___ Combine and save plots ----
plot_pca_brat <- ggarrange(plot_pca_brat_a, plot_pca_brat_b, plot_pca_brat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_brat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_brat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_brat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_brat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

#' ## Plot SCAT PCs in 2D for several variables types

# _3.) Plot SCAT PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_scat_a <- get_pca_plot(expr_data_pca = PCA_SCAT, expr_data_raw = SCAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a", percent_var = SCAT_PV)

# Overall expression differences and obesity status among f1 offspring
plot_pca_scat_b <- get_pca_plot(expr_data_pca = PCA_SCAT, expr_data_raw = SCAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b", percent_var = SCAT_PV)

# Overall expression differences and obesity parental obesity status among f0
plot_pca_scat_c <- get_pca_plot(expr_data_pca = PCA_SCAT, expr_data_raw = SCAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c", percent_var = SCAT_PV)

# ___ Combine and save plots ----

plot_pca_scat <- ggarrange(plot_pca_scat_a, plot_pca_scat_b, plot_pca_scat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_scat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_scat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_scat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_scat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

#' ## Plot EVAT PCs in 2D for several variables types

# _4.) Plot EVAT PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_evat_a <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a", percent_var = EVAT_PV)

# Overall expression differences and obesity status among f1 offspring
plot_pca_evat_b <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b", percent_var = EVAT_PV)

# Overall expression differences and obesity parental obesity status among f0
plot_pca_evat_c <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c", percent_var = EVAT_PV)

# ___ Combine and save plots ----

plot_pca_evat <- ggarrange(plot_pca_evat_a, plot_pca_evat_b, plot_pca_evat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_evat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_evat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_evat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_evat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

#' ## Plot LIAT PCs in 2D for several variables types

# _5.) Plot LIAT PCs in 2D for several variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_liat_a <- get_pca_plot(expr_data_pca = PCA_LIAT, expr_data_raw = LIAT , variable = "Tissue", legend_title = "Tissue", plot_title =  "a", percent_var = LIAT_PV)

# Overall expression differences and obesity status among f1 offspring
plot_pca_liat_b <- get_pca_plot(expr_data_pca = PCA_LIAT, expr_data_raw = LIAT , variable = "ObesityLgcl", legend_title = "F1 Obesity", plot_title =  "b", percent_var = LIAT_PV)

# Overall expression differences and obesity parental obesity status among f0
plot_pca_liat_c <- get_pca_plot(expr_data_pca = PCA_LIAT, expr_data_raw = LIAT , variable = "ObeseParents", legend_title = "F0 Obesity", plot_title =  "c", percent_var = LIAT_PV)

# ___ Combine and save plots ----

plot_pca_liat <- ggarrange(plot_pca_liat_a, plot_pca_liat_b, plot_pca_liat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_liat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_liat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_liat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_liat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

#' # Prepare getting numerical summaries of above PCAs

# Prepare getting numerical summaries of above PCAs ----

# use the PCA graphics and results above to get numerical summaties of hwat can be sees
# follow https://lauren-blake.github.io/Regulatory_Evol/analysis/gene_exp_corr.html
# in section "Correlations"

# use these results to inform DGE
# see DGE section to see aht models will be tested - try to match those modles here
# report those models used here and in DGE and write down the results
# possibly adjust DGE models

#' ## Isolate First 5 Pricipal Components

# _1.) Isolate First 5 Pricipal Components ----

PCs_FLAT <- get_isolated_pcs(PCA_FLAT)
PCs_BRAT <- get_isolated_pcs(PCA_BRAT)
PCs_EVAT <- get_isolated_pcs(PCA_EVAT)
PCs_LIAT <- get_isolated_pcs(PCA_LIAT)
PCs_SCAT <- get_isolated_pcs(PCA_SCAT)

#' ## Get data frames to model variation of first PC against chosen factors

# _2.) Get data frames to model variation of first PC against chosen factors ----

FLAT_md <- get_model_data(FLAT, PCs_FLAT)
BRAT_md <- get_model_data(BRAT, PCs_BRAT)
EVAT_md <- get_model_data(EVAT, PCs_EVAT)
LIAT_md <- get_model_data(LIAT, PCs_LIAT)
SCAT_md <- get_model_data(SCAT, PCs_SCAT)

#' # Getting numerical summaries of above PCAs

# Getting numerical summaries of above PCAs ----

#' ## Analyse FLAT's first PC

# _1.) Analyse FLAT's first PC ---- 

# opening connection to text with full path for Rmarkdown
# sink(file = paste0(here("plots/050_r_array_analysis__text_pca_FLAT.txt")))
# sink(file = "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_FLAT.txt")

#' ### Getting 0-model

# __a) Getting 0-model ----
model_intercept <- lm(PC1 ~  1, data = FLAT_md)

#' ### Testing effect of "Tissue"

# __b) Testing effect of "Tissue"  ----
model_tissue <- lm(PC1 ~  Tissue, data = FLAT_md)
summary(model_tissue) # all tissues distinct
anova(model_intercept, model_tissue) # ***Tissue highly significant structuring PC1*** 

#' ### Testing effect of "ObesityLgcl"

# __c) Testing effect of "ObesityLgcl" ----
model_obesity_off <- lm(PC1 ~  ObesityLgcl, data = FLAT_md)
summary(model_obesity_off) # no signal
anova(model_intercept, model_obesity_off) # offspring's obesity not significant structuring PC1 

#' ### Testing effect of "ObeseParents"

# __d) Testing effect of "ObeseParents" ----
model_obesity_par <- lm(PC1 ~ ObeseParents, data = FLAT_md)
summary(model_obesity_par) # no signal
anova(model_intercept, model_obesity_par) # parents's obesity not significant structuring PC1  

# sink()

#' ## Analyse BRAT's first PC

# _2.) Analyse BRAT's first PC ----  

# opening connection to text with full path for Rmarkdown
# sink(file = paste0(here("plots/050_r_array_analysis__text_pca_BRAT.txt")))
# sink(file = paste0(here("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_BRAT.txt")))

#' ### Getting 0-model

# __a) Getting 0-model ----

model_intercept <- lm(PC1 ~  1, data = BRAT_md)

#' ### Testing effect of "ObesityLgcl"

# __b) Testing effect of "ObesityLgcl" ----

model_obesity_off <- lm(PC1 ~  ObesityLgcl, data = BRAT_md)
summary(model_obesity_off) # no signal
anova(model_intercept, model_obesity_off) # Offspring's obesity not significant structuring PC1 

#' ### Testing effect of "ObeseParents"

# __c) Testing effect of "ObeseParents" ----

model_obesity_par <- lm(PC1 ~ ObeseParents, data = BRAT_md)
anova(model_intercept, model_obesity_par) # ***Parent's obesity significant structuring PC1, check in DGE***

#' ### Testing all reference levels of "ObeseParents"

# __d) Testing all reference levels of "ObeseParents" ----

summary(lm(PC1 ~ ObeseParents, data = BRAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 1))))
summary(lm(PC1 ~ ObeseParents, data = BRAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 2))))
summary(lm(PC1 ~ ObeseParents, data = BRAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 3))))
summary(lm(PC1 ~ ObeseParents, data = BRAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 4))))

# sink()

#' ## Analyse SCAT's first PC 

# _3.) Analyse SCAT's first PC ---- 

# opening connection to text with full path for Rmarkdown
# sink(file = paste0(here("plots/050_r_array_analysis__text_pca_SCAT.txt")))
# sink(file = "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_SCAT.txt")

#' ### Getting 0-model

# __a) Getting 0-model ----

model_intercept <- lm(PC1 ~  1, data = SCAT_md)

#' ### Testing effect of "ObesityLgcl"

# __b) Testing effect of "ObesityLgcl" ----

model_obesity_off <- lm(PC1 ~  ObesityLgcl, data = SCAT_md)
summary(model_obesity_off) # no signal
anova(model_intercept, model_obesity_off) # Offspring's obesity not significant structuring PC1 

#' ### Testing effect of "ObeseParents"

# __c) Testing effect of "ObeseParents" ----

model_obesity_par <- lm(PC1 ~ ObeseParents, data = SCAT_md)
anova(model_intercept, model_obesity_par) # ***Parent's obesity significant structuring PC1, check in DGE***

#' ### Testing all reference levels of "ObeseParents"

# __d) Testing all reference levels of "ObeseParents" ----

summary(lm(PC1 ~ ObeseParents, data = SCAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 1))))
summary(lm(PC1 ~ ObeseParents, data = SCAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 2))))
summary(lm(PC1 ~ ObeseParents, data = SCAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 3))))
summary(lm(PC1 ~ ObeseParents, data = SCAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 4))))

# sink()

#' ## Analyse LIAT's first PC 

# _4.) Analyse LIAT's first PC ---- 

# opening connection to text with full path for Rmarkdown
# sink(file = paste0(here("plots/050_r_array_analysis__text_pca_LIAT.txt")))
# sink(file = "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_LIAT.txt")

#' ### Getting 0-model

# __a) Getting 0-model ----

model_intercept <- lm(PC1 ~  1, data = LIAT_md)

#' ### Testing effect of "ObesityLgcl"

# __b) Testing effect of "ObesityLgcl" ----

model_obesity_off <- lm(PC1 ~  ObesityLgcl, data = LIAT_md)
summary(model_obesity_off) # no signal
anova(model_intercept, model_obesity_off) # Offspring's obesity not significant structuring PC1 

#' ### Testing effect of "ObeseParents"

# __c) Testing effect of "ObeseParents" ----

model_obesity_par <- lm(PC1 ~ ObeseParents, data = LIAT_md)
anova(model_intercept, model_obesity_par) # ***Parent's obesity significant structuring PC1, check in DGE***

#' ### Testing all reference levels of "ObeseParents"

# __d) Testing all reference levels of "ObeseParents" ----

summary(lm(PC1 ~ ObeseParents, data = LIAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 1))))
summary(lm(PC1 ~ ObeseParents, data = LIAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 2))))
summary(lm(PC1 ~ ObeseParents, data = LIAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 3))))
summary(lm(PC1 ~ ObeseParents, data = LIAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 4))))

# sink()

#' ## Analyse EVAT's first PC

# _5.) Analyse EVAT's first PC ----

# opening connection to text with full path for Rmarkdown
# sink(file = paste0(here("plots/050_r_array_analysis__text_pca_EVAT.txt")))
# sink(file = "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_EVAT.txt")

#' ### Getting 0-model

# __a) Getting 0-model ----

model_intercept <- lm(PC1 ~  1, data = EVAT_md)

#' ### Testing effect of "ObesityLgcl"

# __b) Testing effect of "ObesityLgcl" ----

model_obesity_off <- lm(PC1 ~  ObesityLgcl, data = EVAT_md)
summary(model_obesity_off) # no signal
anova(model_intercept, model_obesity_off) # Offspring's obesity not significant structuring PC1 

#' ### Testing effect of "ObeseParents"

# __c) Testing effect of "ObeseParents" ----

model_obesity_par <- lm(PC1 ~ ObeseParents, data = EVAT_md)
anova(model_intercept, model_obesity_par) # ***Parent's obesity significant structuring PC1, check in DGE***

#' ### Testing all reference levels of "ObeseParents"

# __d) Testing all reference levels of "ObeseParents" ----

summary(lm(PC1 ~ ObeseParents, data = EVAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 1))))
summary(lm(PC1 ~ ObeseParents, data = EVAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 2))))
summary(lm(PC1 ~ ObeseParents, data = EVAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 3))))
summary(lm(PC1 ~ ObeseParents, data = EVAT_md %>% mutate(ObeseParents = relevel(ObeseParents, 4))))

# sink()

#' # Re-implement analysis of array intensities

# Re-implement analysis of array intensities ----

#' ## Shape and check array intensity data

# _1.) Shape and check array intensity data ----

# see https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html

#' ###  Check input data formats

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

#' ### Covert expression set to data table for inspection

# __b) Covert expression set to data table for inspection ----

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

#' ### Inspect expression data raw intensities distribution

# __c) Inspect expression data raw intensities distribution (defunct)  ----

# to check that equal amounts of data are available for comparison - they are not

# *** plotting commented out until melting command above is adjusted ****

# ggplot(FLAT_DT.m1) +
#   # geom_density(aes(Intensity), stat = "bin", bins = 200) +
#   geom_density(aes(Intensity, colour = ObesityLgcl), stat = "bin", bins = 200) +
#   ylab("Number of measurments across all other variables") +
#   xlab("Intensity") +
#   ggtitle("Intensities' availibilty and distribution for offsprings obesity") +
#   facet_wrap(.~Tissue) + 
#   theme_bw()

# ggplot(FLAT_DT.m1) +
#   # geom_density(aes(Intensity), stat = "bin", bins = 200) +
#   geom_density(aes(Intensity, colour = ObeseParents), stat = "bin", bins = 200) +
#   ylab("Number of measurments across all other variables") +
#   xlab("Intensity") +
#   ggtitle("Intensities' availibilty and distribution for parents obesity") +
#   facet_wrap(.~Tissue) + 
#   theme_bw()

#' ### Inspect expression data raw intensities' density

# __d) Inspect expression data raw intensities' density  ----

# To check if distributions are different - hopefully they are a bit - yes perhaps in EVAT when Mother and Father are not Obese

# ggplot(FLAT_DT.m1) +
#   # geom_density(aes(Intensity), stat = "bin", bins = 200) +
#   geom_density(aes(Intensity, colour = ObesityLgcl), stat = "density") +
#   ylab("Density of measurments across all other variables") +
#   xlab("Intensity") +
#   ggtitle("Intensities' density for offsprings obesity") +
#   facet_wrap(.~Tissue) + 
#   theme_bw()

# ggplot(FLAT_DT.m1) +
#   # geom_density(aes(Intensity), stat = "bin", bins = 200) +
#   geom_density(aes(Intensity, colour = ObeseParents), stat = "density") +
#   ylab("Denisty of measurments across all other variables") +
#   xlab("Intensity") +
#   ggtitle("Intensities' density for parents obesity") +
#   facet_wrap(.~Tissue) + 
#   theme_bw()

#' ## DGE analysis using Limma

# _2.) DGE analysis using Limma ----

# see DESeq2 tutorial at https://colauttilab.github.io/RNA-Seq_Tutorial.html
# see GCSCore tutaorial at https://www.bioconductor.org/packages/release/bioc/vignettes/GCSscore/inst/doc/GCSscore.pdf
# see Limma slides at https://s3.amazonaws.com/assets.datacamp.com/production/course_6456/slides/chapter1.pdf
# **use this!** - Limma slides at  https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
# see also here for LFC shrinkage and genral procedure https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/08_practical_DE.pdf
# see here why shrinking LFCs is considered unenecssary by the limma authors - https://support.bioconductor.org/p/100804/

#' ###  Test for DGE for each tissue against all others, using FLAT

# __a) Test for DGE for each tissue against all others, using FLAT ----

# As justified per PCA results 

FLAT_Tissue_TopTableList <- get_dge_for_individal_tissues(FLAT)

#' ### Test for DGE among obese and non-obese offspring

# __b)  Test for DGE among obese and non-obese offspring ----

# No PCA signal (see above) nor DGE detected across all tissues or in any tissue based on offsprings' obesity status
# No further work necessary - but report!

# get_dge_for_offspring_obesity(FLAT)
# get_dge_for_offspring_obesity(BRAT)
# get_dge_for_offspring_obesity(LIAT)
# get_dge_for_offspring_obesity(SCAT)
# get_dge_for_offspring_obesity(EVAT)

#' ### Test for DGE among offspring based on parental obesity

# __c)  Test for DGE among offspring based on parental obesity ----

# Defining and applying contrasts: One of "MotherFatherNotObese", "FatherObese", "MotherFatherObese", or "MotherObese" 
#  against all remaining three levels.
#  Compare to PCA results, including coefficients - Expression is variable based on tissue, and within each tissue based on parental obesity

FLAT_TopTableList <- get_dge_for_parent_obesity(FLAT) # analysis off all tissues together doesn't really make doesn't really make sense
EVAT_TopTableList <- get_dge_for_parent_obesity(EVAT) # likely not needed - see manuscript results 05.07.2023
SCAT_TopTableList <- get_dge_for_parent_obesity(SCAT) # likely not needed - see manuscript results 05.07.2023
BRAT_TopTableList <- get_dge_for_parent_obesity(BRAT) 
LIAT_TopTableList <- get_some_dge_for_parent_obesity(LIAT) 

#' ###  Choose DGE results for further analyses

# __d) Choose DGE results for further analyses ----

# ___ FLAT - keeping all contrasts for all tissues ----

# see file "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_FLAT.txt"

names(FLAT_Tissue_TopTableList)
          
# ___ FLAT, BRAT, SCAT, LIAT, EVAT - not keeping any contrasts defined by offspring' obesity  ----

#  see PCA results: 
#  "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_EVAT.txt" # likely not needed - see manuscript results 05.07.2023
#  "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_SCAT.txt" # likely not needed - see manuscript results 05.07.2023
#  "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_BRAT.txt"
#  "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_LIAT.txt"

#  see DGE results - above

# ___ EVAT - keeping some contrasts defined by parents' obesity  ----

# no contrasts ar needed  
# - see PCA results: "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_EVAT.txt"
# - see manuscript results 05.07.2023

# ___ SCAT - keeping some contrasts defined by parents' obesity  ----

# no contrasts ar needed  
# - see PCA results: "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_SCAT.txt"
# - see manuscript results 05.07.2023

# ___ BRAT - keeping some contrasts defined by parents' obesity  ----

# some contrast are needed: 
# - see PCA results: "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_BRAT.txt"
# - see manuscript results 05.07.2023
# - "MotherFatherObese vs FatherObese"

names(BRAT_TopTableList) # check available slots 
names(BRAT_TopTableList[c(2)]) # select slots corresponding to contrasts readout above 

BRAT__Select_TopTableList <- BRAT_TopTableList[c(2)] # selecting individual and compound contrasts

# ___ LIAT - keeping some contrasts defined by parents' obesity  ----

# some contrast are needed: 
# - see PCA results: "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_LIAT.txt"
# - see manuscript results 05.07.2023
# - "MotherFatherNotObese vs MotherFatherObese"

names(LIAT_TopTableList)
names(LIAT_TopTableList[c(2)]) 

LIAT__Select_TopTableList <- LIAT_TopTableList[c(2)] 

#' ### Compile a well-labelled list with all DGE results

# __e)  Compile a well-labelled list with all DGE results  ----

# - these results are likley not needed
names(FLAT_Tissue_TopTableList) <- paste0("FLAT - ", names(FLAT_Tissue_TopTableList))

# - these results are undefined, as PCA didn't show structure and no DEGs were define (see 5.-Jul.2023)
#   names(EVAT__Select_TopTableList) <- paste0("EVAT - ", names(EVAT__Select_TopTableList))
#   names(SCAT__Select_TopTableList) <- paste0("SCAT - ", names(SCAT__Select_TopTableList))

names(BRAT__Select_TopTableList) <- paste0("BRAT - ", names(BRAT__Select_TopTableList))
names(LIAT__Select_TopTableList) <- paste0("LIAT - ", names(LIAT__Select_TopTableList))

FULL_TopTableList <- c(
  # FLAT_Tissue_TopTableList,
  # EVAT__Select_TopTableList, 
  BRAT__Select_TopTableList,
  LIAT__Select_TopTableList
  )

# receiving a table with 2 slots, each containing DGE results fo a specific tissue and statistically relavent contrasts
names(FULL_TopTableList)

# __f) Filter Toptables further if required ----   

# check tables
FULL_TopTableList[[1]]
FULL_TopTableList[[2]]

# check for NAs
FULL_TopTableList[[1]][rowSums(is.na(FULL_TopTableList[[1]])) > 0, ]
FULL_TopTableList[[2]][rowSums(is.na(FULL_TopTableList[[2]])) > 0, ]

# Incpect data in question - green tail on Volcano plot
FULL_TopTableList[[1]][ which(  log2(FULL_TopTableList[[1]]$logFC) > 1 ) , ]
FULL_TopTableList[[1]][ which(  log2(FULL_TopTableList[[1]]$logFC) > 1 ) , ]


#' ### Get, assort, arrange, and save Vulcano plots

# __g)  Get, assort, arrange, and save Vulcano plots  ----

FULL_VolcanoPlots <- mapply(get_one_volcanoplot, TopTableListItem = FULL_TopTableList, TopTableListItemName = names(FULL_TopTableList), SIMPLIFY = FALSE)

# the following three sets are undefined or serve no purpose (see 5-Jul-2023)
#  FLAT_VolcanoPlots <- FULL_VolcanoPlots[grep("FLAT", names(FULL_VolcanoPlots))]
#  SCAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^SCAT", names(FULL_VolcanoPlots))]
#  EVAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^EVAT", names(FULL_VolcanoPlots))]

BRAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^BRAT", names(FULL_VolcanoPlots))]
LIAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^LIAT", names(FULL_VolcanoPlots))]

# the following three sets are undefined or serve no purpose (see 5-Jul-2023)
#  FLAT_VolcanoPlotsComposite <- ggarrange(plotlist = FLAT_VolcanoPlots, ncol = 2, nrow = 2, labels = "auto")
#  SCAT_VolcanoPlotsComposite <- ggarrange(plotlist = SCAT_VolcanoPlots, ncol = 1, nrow = 1, labels = NULL)
#  EVAT_VolcanoPlotsComposite <- ggarrange(plotlist = EVAT_VolcanoPlots, ncol = 3, nrow = 1, labels = "auto")
BRAT_VolcanoPlotsComposite <- ggarrange(plotlist = BRAT_VolcanoPlots, ncol = 1, nrow = 1, labels = NULL)
LIAT_VolcanoPlotsComposite <- ggarrange(plotlist = LIAT_VolcanoPlots, ncol = 1, nrow = 1, labels = NULL)
BRLI_VolcanoPlotsComposite <- ggarrange(plotlist = c(BRAT_VolcanoPlots, LIAT_VolcanoPlots), ncol = 2, nrow = 1, labels = "auto")

# ggsave(plot = FLAT_VolcanoPlotsComposite, path = here("plots"), 
#        filename = "050_r_array_analysis__plot_volcano_flat.pdf",  
#        width = 180, height = 200, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
# 
# ggsave(plot = SCAT_VolcanoPlotsComposite, path = here("plots"), 
#        filename = "050_r_array_analysis__plot_volcano_scat.pdf",  
#        width = 180, height = 150, units = "mm", dpi = 150,  limitsize = TRUE, scale = 1.3)
# 
# ggsave(plot = EVAT_VolcanoPlotsComposite, path = here("plots"), 
#        filename = "050_r_array_analysis__plot_volcano_evat.pdf",  
#        width = 300, height = 125, units = "mm", dpi = 200,  limitsize = TRUE, scale = 3)

ggsave(plot = BRAT_VolcanoPlotsComposite, path = here("plots"), 
       filename = "050_r_array_analysis__plot_volcano_brat.pdf",  
       width = 180, height = 200, units = "mm", dpi = 300,  limitsize = TRUE, scale = 1)

ggsave(plot = LIAT_VolcanoPlotsComposite, path = here("plots"), 
       filename = "050_r_array_analysis__plot_volcano_liat.pdf",  
       width = 180, height = 200, units = "mm", dpi = 300,  limitsize = TRUE, scale = 1)

ggsave(plot = BRLI_VolcanoPlotsComposite, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_volcano_brli.pdf",  
       width = 160, height = 80, units = "mm", dpi = 300,   limitsize = TRUE, scale = 2.5)

#' ### Get, assort, arrange, and save heat maps (drafted)

# __h) Get, assort, arrange, and save heat maps (drafted) ----

# code below very dirty, could be cleaned out and made into a function
# sample annotation can be improved using a data frame as shown here
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

# ___ BRAT ----

# get expression data as matrix with probe ids, subset for speed 
foo <- as_tibble(exprs(BRAT), rownames = c("PROBEID")) %>% filter(PROBEID %in%  BRAT__Select_TopTableList[[1]][["PROBEID"]])
foo_mat <- base::as.matrix(foo %>% select(-PROBEID))
rownames(foo_mat) <- (foo %>% pull(PROBEID))

# isolate the DEG top table
bar <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]]

# join DEG top table and expression data to a new tibble  
foobar <- left_join(bar, foo)

# keep only relavent logFC, conforming with Volcano plots, p values do not need further adjustments
foobar <- foobar %>% filter(logFC < -1 | logFC > 1) %>% arrange(logFC)
foobar <- foobar %>% filter(adj.P.Val < 0.05 )%>% arrange(logFC)


# pheatmap needs a matrix - hence converting tibble to matrix
foo_mat <- base::as.matrix (foobar %>% select(contains("_BRAT")))
rownames(foo_mat) <- foobar[["SYMBOL"]]

# to add more sample annotations to heat map matrix colnames - isolating this data here
mdata_tibble <- tibble(pData(BRAT)) %>% mutate(SAMPLE = paste0(Animal,"_",Tissue))

# checking possibilty to extend matrix column names
which(colnames(foo_mat) %in% mdata_tibble$SAMPLE)

# extending matrix column names (samples) with relevant metadata
colnames(foo_mat) <- paste0(colnames(foo_mat), " | " , 
       mdata_tibble[which(colnames(foo_mat) %in% mdata_tibble$SAMPLE) , "ObeseParents"][["ObeseParents"]], " | " ,
       mdata_tibble[which(colnames(foo_mat) %in% mdata_tibble$SAMPLE) , "ObesityLgcl"][["ObesityLgcl"]]
       )

# adding stars to contrast
colnames(foo_mat)[grep(pattern = "MotherFatherObese|FatherObese", colnames(foo_mat))] <- paste(colnames(foo_mat)[grep(pattern = "MotherFatherObese|FatherObese", colnames(foo_mat))], "*")

# print heat map to script and file
pheatmap(foo_mat, scale = "row")
pheatmap(foo_mat, scale = "row", filename =  paste0(here("plots"),"/","050_r_array_analysis__plot_heatmap_brat.pdf"))

# ___ LIAT ----

# get expression data as matrix with probe ids, subset for speed 
foo <- as_tibble(exprs(LIAT), rownames = c("PROBEID")) %>% filter(PROBEID %in%  LIAT__Select_TopTableList[[1]][["PROBEID"]])
foo_mat <- base::as.matrix(foo %>% select(-PROBEID))
rownames(foo_mat) <- (foo %>% pull(PROBEID))

# isolate the DEG top table
bar <- LIAT__Select_TopTableList[["LIAT - MotherFatherObese vs MotherFatherNotObese"]]

# join DEG top table and expression data to a new tibble  
foobar <- left_join(bar, foo)

# keep only relavent logFC
foobar <- foobar %>% filter(logFC < -1 | logFC > 1) %>% arrange(logFC)
foobar <- foobar %>% filter(adj.P.Val < 0.05 ) %>% arrange(logFC)

# pheatmap needs a matrix - hence converting tibble to matrix
foo_mat <- base::as.matrix (foobar %>% select(contains("_LIAT")))
rownames(foo_mat) <- foobar[["SYMBOL"]]

# to add more sample annotations to heat map matrix colnames - isolating this data here
mdata_tibble <- tibble(pData(LIAT)) %>% mutate(SAMPLE = paste0(Animal,"_",Tissue))

# checking possibilty to extend matrix column names
which(colnames(foo_mat) %in% mdata_tibble$SAMPLE)

# extending matrix column names (samples) with relevant metadata
colnames(foo_mat) <- paste0(colnames(foo_mat), " | " , 
                            mdata_tibble[which(colnames(foo_mat) %in% mdata_tibble$SAMPLE) , "ObeseParents"][["ObeseParents"]], " | " ,
                            mdata_tibble[which(colnames(foo_mat) %in% mdata_tibble$SAMPLE) , "ObesityLgcl"][["ObesityLgcl"]]
)

# adding stars to contrast
colnames(foo_mat)[grep(pattern = "MotherFatherObese|MotherFatherNotObese", colnames(foo_mat))] <- paste(colnames(foo_mat)[grep(pattern = "MotherFatherObese|MotherFatherNotObese", colnames(foo_mat))], "*")

# print heat map to script and file
pheatmap(foo_mat, scale = "row")
pheatmap(foo_mat, scale = "row", filename =  paste0(here("plots"),"/","050_r_array_analysis__plot_heatmap_liat.pdf"))

#' ### Save DGE lists

# __i) Save DGE lists ----

BRAT_TTL_sign <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]] %>% filter(adj.P.Val < 0.05 )
BRAT_TTL_uprg <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]] %>% filter(logFC > 1)
BRAT_TTL_down <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]] %>% filter(logFC < -1)

LIAT_TTL_sign <- LIAT__Select_TopTableList[["LIAT - MotherFatherObese vs MotherFatherNotObese"]] %>% filter(adj.P.Val < 0.05 )
LIAT_TTL_uprg <- LIAT__Select_TopTableList[["LIAT - MotherFatherObese vs MotherFatherNotObese"]] %>% filter(logFC > 1)
LIAT_TTL_down <- LIAT__Select_TopTableList[["LIAT - MotherFatherObese vs MotherFatherNotObese"]] %>% filter(logFC < -1)

openxlsx::write.xlsx(BRAT_TTL_sign, paste0(here("tables"), "/050_r_array_analysis_", "BRAT_TTL_sign", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(BRAT_TTL_down, paste0(here("tables"), "/050_r_array_analysis_", "BRAT_TTL_down", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(BRAT_TTL_uprg, paste0(here("tables"), "/050_r_array_analysis_", "BRAT_TTL_uprg", ".xlsx"), asTable = TRUE, overwrite = TRUE)

openxlsx::write.xlsx(LIAT_TTL_sign, paste0(here("tables"), "/050_r_array_analysis_", "LIAT_TTL_sign", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(LIAT_TTL_down, paste0(here("tables"), "/050_r_array_analysis_", "LIAT_TTL_down", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(LIAT_TTL_uprg, paste0(here("tables"), "/050_r_array_analysis_", "LIAT_TTL_uprg", ".xlsx"), asTable = TRUE, overwrite = TRUE)

#' ## Gene Set Enrichment Analysis (GSEA)

# _3.) Gene Set Enrichment Analysis (GSEA) ----

# the following resources are susefule 
# step-by-step GSEA https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/08_practical_DE.pdf
# step-by-step GSEA https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# why shrinkage of LFCs is not necessary, or already done https://support.bioconductor.org/p/100804/

# check test data - possibly used tissue-specific data sets

#' ### Lookup Entrez ID for Clusterprofiler

# __a) Lookup Entrez ID for Clusterprofiler ----

FULL_TopTableListAppended <- lapply(FULL_TopTableList, get_entrez_ids) 

#' ### Get plots of KEEG pathways

# __b) Get plots of KEEG pathways ----

FULL_KeggPlots <- mapply(get_one_kegg_dotplot, TopTableListItem = FULL_TopTableListAppended, TopTableListItemName = names(FULL_TopTableListAppended), SIMPLIFY = FALSE)
names(FULL_KeggPlots)

# FULL_KeggPlots[1] - no result
FULL_KeggPlots[2]

#' ### Save plots of KEEG pathways

# __c) Save plots of KEEG pathways ----

mapply(save_kegg_plots, ggplot_list_item = FULL_KeggPlots, ggplot_list_item_name = names(FULL_KeggPlots), SIMPLIFY = FALSE)


#' ### Implement GO analysis

# __d) Implement GO analysis ----

FULL_GoPlots <- mapply(get_one_go_plot, TopTableListItem = FULL_TopTableListAppended, TopTableListItemName = names(FULL_TopTableListAppended), SIMPLIFY = FALSE)
names(FULL_GoPlots)

#' ### Show and save plots of GO pathways

# __e) Show and save plots of GO pathways ----

FULL_GoPlots[[1]]
FULL_GoPlots[[2]]

mapply(save_go_plots, ggplot_list_item = FULL_GoPlots, ggplot_list_item_name = names(FULL_KeggPlots), SIMPLIFY = FALSE)

ggsave(plot = ggarrange(plotlist =  FULL_GoPlots, ncol = 2, labels = "auto"), filename = "050_r_array_analysis__plot_go_both.pdf", path = here("../manuscript/display_items/"), width = 280, height = 420, unit = "mm", scale = 1)


# #' Experimental: DGE-analysis using GAMs - Investigate overall tissue specific expression differences based on obesity variables

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
model_data <- subset(FLAT_DT.m1, ArrayTarget %in% ArrayTargets) # subest to 100 of 20000 genes to speed up compuaterion
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

#' # Snapshot environment

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "050_r_array_analysis.RData"))
renv::snapshot()

# >>> Reference code from AH ----

# AH code below - UpSet plots ----

# # Upset plots ############################# 
# #hier am Beispiel für Liver
# para <- "Liver"
# 
# sub1 <- subset(Liv_CD_HFD_VS_CD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
# colnames(sub1) <- c("RN","CD, HFD vs. CD, CD")
# sub2 <- subset(Liv_HFD_CD_VS_CD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
# colnames(sub2) <- c("RN","HFD, CD vs. CD, CD")
# sub3 <- subset(Liv_HFD_HFD_VS_CD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
# colnames(sub3) <- c("RN","HFD, HFD vs. CD, CD")
# sub4 <- subset(Liv_HFD_CD_VS_CD_HFD, Col=="#d72323" | Col=="blue")[,c(1,8)]
# colnames(sub4) <- c("RN","HFD, CD vs. CD, HFD")
# sub5 <- subset(Liv_HFD_HFD_VS_CD_HFD, Col=="#d72323" | Col=="blue")[,c(1,8)]
# colnames(sub5) <- c("RN","HFD, HFD vs. CD, HFD")
# sub6 <- subset(Liv_HFD_HFD_VS_HFD_CD, Col=="#d72323" | Col=="blue")[,c(1,8)]
# colnames(sub6) <- c("RN","HFD, HFD vs. HFD, CD")
# 
# mTAB <- merge(merge(merge(merge(merge(
#   sub1,
#   sub2, all = TRUE,by="RN"),
#   sub3, all = TRUE,by="RN"),
#   sub4, all = TRUE,by="RN"),
#   sub5, all = TRUE,by="RN"),
#   sub6, all = TRUE,by="RN")
# 
# mTAB$`CD, HFD vs. CD, CD` <- ifelse(is.na(mTAB$`CD, HFD vs. CD, CD`),  0,1)
# mTAB$`HFD, CD vs. CD, CD` <- ifelse(is.na(mTAB$`HFD, CD vs. CD, CD`),  0,1)
# mTAB$`HFD, HFD vs. CD, CD` <- ifelse(is.na(mTAB$`HFD, HFD vs. CD, CD`),  0,1)
# mTAB$`HFD, CD vs. CD, HFD` <- ifelse(is.na(mTAB$`HFD, CD vs. CD, HFD`),  0,1)
# mTAB$`HFD, HFD vs. CD, HFD` <- ifelse(is.na(mTAB$`HFD, HFD vs. CD, HFD`),  0,1)
# mTAB$`HFD, HFD vs. HFD, CD` <- ifelse(is.na(mTAB$`HFD, HFD vs. HFD, CD`),  0,1)
# 
# pdf(paste0("Results_05.23/",para,"_Upset.pdf"),width = 9,height = 5) 
#   UpSetR::upset(mTAB,nsets = 6, decreasing=TRUE,order.by = "freq",group.by="degree",
#               sets.bar.color=c("maroon","blue","orange","purple","steelblue","red"))
# dev.off()

