#' ---
#' title: "Mice Mating Study"
#' subtitle: "Re-analysis of array data for first revision."
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

# Prepare environment ---- 


# _1.) Collect garbage ----

rm(list=ls())
gc()

# _2.) Packages ----

library(here)
library(renv)

library(magrittr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(tidyr)

library(data.table)
library(ggpubr)
library(writexl) # write excel sheets

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

# for GO analysis
library(limma)
library(gplots)
library(RColorBrewer)
library(KernSmooth)
library(FactoMineR)

# for cleaning, mostly unused

suppressMessages(require(tidyverse))
suppressMessages(require(Seurat))
suppressMessages(require(cowplot))
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(igraph))

library(sechm) # heatmap 



# _3.) Encode plotting colours  ----

hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(n=11, name="RdBu"))(50))

# _4.) Functions ----

# Updates metadata accompanying expression data (last updated 27-May-2024)
get_adjusted_array_data = function(expression_set, model_variables) {
  
  require(dplyr)
  require(Biobase)
  
  # stop("Remove function building code")
  # expression_set <- FLAT
  # model_variables <- mice_f1_modeled_data_with_rna_seq_data
  
  # adjust column and row names names in expression data
  colnames(expression_set) <- colnames(expression_set) %>%
    str_replace_all("bAT",    "BRAT") %>%
    str_replace_all("ingWAT", "IWAT") %>%
    str_replace_all("liver",  "LIVT") %>%
    str_replace_all("eWAT",   "EVAT")
  
  # adjust "Tissue" column values  in expression data
  pData(expression_set) %<>% mutate(
    Tissue = case_when(
      Tissue == "bAT"    ~ "BRAT",
      Tissue == "ingWAT" ~ "IWAT",
      Tissue == "liver"  ~ "LIVT",
      Tissue == "eWAT"   ~ "EVAT"
    )
  )
  
  # merge obesity variables from modelling  to metadata from array experiments
  pData(expression_set) <- left_join((pData(expression_set) %>% dplyr::select(
    -c(Sex, Parental_diet, Diet_mother, Diet_father, Group)
  )),
  (model_variables),
  by = c("Animal" = "AnimalId"))
  
  # convert from Expression Set to Summarized Experiment object
  warning("Converting from Expression Set object to Summarized Experiment object")
  rownames(pData(expression_set)) <- paste0(pData(expression_set)[["Animal"]] ,"_",  pData(expression_set)[["Tissue"]])
  colnames(assayData(expression_set)$exprs)
  stopifnot(identical(rownames(pData(expression_set)), colnames(assayData(expression_set)$exprs)))
  expression_set <- as(expression_set, "SummarizedExperiment")
  
  # Get uppercase symbols
  rowData(expression_set)[["SYMBOL"]]  <- toupper(rowData(expression_set)[["SYMBOL"]])
  
  return(expression_set)
  
}

# Some basic data checks in section "Clean out expression data" (last updated 28-May-2024)
plot_mean_expression_and_variance = function(se_object, plot_label){
  
  # Calculate gene mean across cell
  gene_mean <- rowMeans(assays(se_object)$exprs) 
  
  # Calculate gene variances
  gene_var  <- rowVars(assays(se_object)$exprs)  
  
  # get plotting data frame
  gene_stat_df <- tibble(gene_mean,gene_var)
  
  #ggplot plot
  plot <- ggplot(data = gene_stat_df, aes(x = log(gene_mean), y = log(gene_var))) + 
    geom_point(size=0.5) +
    theme_minimal() +
    xlab("log expression means") +
    ylab("log expression variances")
  
  return(plot)
}

# Some basic data checks in section "Clean out expression data" (last updated 28-May-2024)
plot_rare_genes = function(se_object, treshhold = 0.5){
  
  # for function testing only
  # se_object <- BRAT
  
  # Calculate gene mean across cell
  gene_mean <- rowMeans(assays(se_object)$exprs) 
  
  abundant_genes <- gene_mean >= treshhold # Remove Low abundance genes
  
  message(paste("Among", length(gene_mean), "genes", length(abundant_genes) , "are above a treshhold of", treshhold, "."))
  
  # plot low abundance gene filtering - histogram
  hist(log10(gene_mean), breaks= 100 , main="", col="grey80",
       xlab = expression(Log[10] ~ "average count"))
  # plot low abundance gene filtering - red line
  abline( v = log10(treshhold), col="red", lwd = 2, lty = 2)
  
}

# Get principal components of expression data (modified 28-May-2024 to use SummarizedExperiment)
get_principal_components = function(SummarizedExperiment){
  
  # Building initial models
  message(paste0("\nGetting PCs of data set \"", deparse(substitute(SummarizedExperiment))), "\".")
  
  return( prcomp( t( assays(SummarizedExperiment)$exprs), scale = TRUE, center = TRUE))
  
}

# Get percentage of variation of the principal components for plotting. (modified 28-May-2024 to use SummarizedExperiment)
get_percent_variation = function(prcomp_outout) {
  
  # Building initial models
  message(paste0("\nGetting percent variation of PC set \"", deparse(substitute(prcomp_outout))), "\".")
  
  return(round(100 * prcomp_outout$sdev^2 / sum(prcomp_outout$sdev^2), 1) )
  
}

# Get q-mode principal components of expression data (modified 28-May-2024 to use SummarizedExperiment)
get_q_mode_principal_components  = function(SummarizedExperiment){
  
  # Building initial models
  message(paste0("\nGetting PCs of data set \"", deparse(substitute(SummarizedExperiment))), "\".")
  
  return( prcomp(( assays(SummarizedExperiment)$exprs), scale = FALSE, center = TRUE))
  
}

# Get PCA plots
get_pca_plot = function(expr_data_pca, expr_data_raw, variable, legend_title, plot_title, percent_var) {
  
  require("factoextra")
  require("ggplot2")
  
  pca_plot <- fviz_pca_ind(
    expr_data_pca,
    label = "none",
    habillage = colData(expr_data_raw)[[variable]],
    pointsize = 2,
    legend.title = legend_title,
    invisible = "quali",
    addEllipses = TRUE,
    ellipse.level=0.95,
    title = plot_title
  ) +
    labs(
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    ) +
    theme_bw() +
    scale_shape_manual(values = c(19, 19, 19, 19, 19))
  
  return(pca_plot)
  
}

# Get q-mode PCA plots
get_q_mode_pca_plot = function(expr_data_pca, expr_data_raw, variable, plot_title, percent_var) {
  
  require("factoextra")
  require("ggplot2")
  
  # stop("remove function buiulding code")
  # expr_data_pca = qPCA_BRAT
  # expr_data_raw = BRAT
  # variable = "ParentalDietMoFa"
  # plot_title =  "a"
  # percent_var = "qPV_BRAT"
  # 
  # build plot as 
  pca_plot <- fviz_pca_var(expr_data_pca, col.var = "cos2", repel = TRUE, title = plot_title) + 
    theme_bw() + labs(
      x = paste0("PC1: ", percent_var[1], "% variance"),
      y = paste0("PC2: ", percent_var[2], "% variance")
    )
  
  # edit axis labels
  
  
  # - open plot for manual editing
  pca_plot_built <- ggplot_build(pca_plot)
  
  # - edit plot
  pca_plot_built$data[[1]]$label <- colData(expr_data_raw)[ which( rownames(colData(expr_data_raw)) %in% pca_plot_built$data[[1]]$label), ][[variable]]
  
  # - close plot after manual editing for plotting 
  pca_plot_edited <- plot_grid(ggplot_gtable(pca_plot_built))
  
  
  return(pca_plot_edited)
  
}

# All functions below are unrevised ----

# Convenience function to define DGE among offspring tissue types
get_dge_for_individal_tissues =  function(ExpSet){
  
  message("Convenience function to test DGE for each tissue compred to all other tissues, only makes sense for a data set with different tissue type. Analysis of factor \"Tissue\" is hard-coded.")
  
  # Building initial models
  message(paste0("\nUsing data set \"", deparse(substitute(ExpSet))), "\".")   
  
  
  message("\nConsider whether the sample size is appropriate for DGE, you should have at least six samples per group:")
  
  pData(ExpSet) %>% dplyr::select(ObeseParents) %>% table() %>% print()
  
  # Building initial models
  message("Building initial model.")   
  
  design_tissue_types <- model.matrix(~ ExpSet[["Tissue"]] - 1)
  colnames(design_tissue_types) <-c("BRAT", "EVAT", "LIVT", "IWAT") 
  fit_tissue_types <- lmFit(ExpSet, design_tissue_types)
  
  # Setting contrasts - likley not all will give results 
  contrast_names <- c(
    "BRAT vs EVAT & LIVT & IWAT",
    "EVAT vs BRAT & LIVT & IWAT",
    "LIVT vs EVAT & BRAT & IWAT",
    "IWAT vs EVAT & BRAT & LIVT")
  
  message("\nDefining contrasts: \"", paste0(contrast_names, collapse = "\", \""), "\".") 
  
  contrast_list <- vector(mode = "list", length = length(contrast_names))
  
  contrast_list[[1]]  <- makeContrasts("BRAT vs EVAT & LIVT & IWAT" =  BRAT - (EVAT + LIVT + IWAT)/3, levels = design_tissue_types)
  contrast_list[[2]]  <- makeContrasts("EVAT vs BRAT & LIVT & IWAT" =  EVAT - (BRAT + LIVT + IWAT)/3, levels = design_tissue_types)
  contrast_list[[3]]  <- makeContrasts("LIVT vs EVAT & BRAT & IWAT" =  LIVT - (EVAT + BRAT + IWAT)/3, levels = design_tissue_types)
  contrast_list[[4]]  <- makeContrasts("IWAT vs EVAT & BRAT & LIVT" =  IWAT - (EVAT + BRAT + LIVT)/3, levels = design_tissue_types)
  
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

# Convenience function to test parents' obesity status in all tissue samples but LIVT
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

# Convenience function to test parents' obesity status in LIVT samples
get_some_dge_for_parent_obesity = function(ExpSet){
  
  message("Reduced contrast set for LIVT analysis") 
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
  
  # stop("remove function development code")
  # TopTable <- FULL_TopTableList[[2]]
  # rm(TopTable)
  # rm(translations)
  # rm(TopTableAppended)

  message(paste0("Looking up Entrez IDs for data set, storing in column \"ENTREZID\"."))  
  
  # create target object
  TopTableAppended <- TopTable
  
  # lookup Entrez ID from RefSeq ID and store both columns for separtaion
  translations <- bitr(geneID = TopTable[["ID"]], fromType = "REFSEQ", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = FALSE)
  
  # message("Debugging code")
  # nrow(translations)
  # head(translations)
  # nrow(TopTableAppended)
  # head(TopTableAppended)
  
  # Move second column to traget obeject 
  # Doesn't always work: 
  #   TopTableAppended[["ENTREZ"]] <- translations[ , 2]
  # Using dplyr instead: 
  
  message("Joing with \"by = c(\"ID\" = \"REFSEQ\")\": Hardcoded \"REFSEQ\" may need to be chanaged if looking at transcripts.")
  
  TopTableAppended <- left_join(as_tibble(TopTableAppended), as_tibble(translations),  by = c("ID" = "REFSEQ"), keep = FALSE)
  
  return(TopTableAppended)
          
}

# KEGG plot creation
get_one_kegg_dotplot <- function(TopTableListItem, TopTableListItemName, save_to_disk = TRUE, table_path = NULL) {
    
  # see https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
  
  # packages
  require("clusterProfiler")
  require("enrichplot")
  require("stringr")
  
  # for function building only
  # stop("Remove function building code")
  # TopTableListItem = FULL_TopTableListAppended[[1]]
  # TopTableListItemName = names(FULL_TopTableListAppended)[[1]]
  # save_to_disk = TRUE
  # table_path = NULL
  # rm(list(TopTableListItem, TopTableListItemName))
  
  # diagnostic
  message(paste0("Creating KEGG plot for data set: \"", TopTableListItemName, "\".", sep = ""))
  
  # look- up KEGG pathways
  kegg_result <- enrichKEGG(gene = TopTableListItem$ENTREZID, keyType = 'ncbi-geneid', organism = 'mmu', minGSSize = 5)
    
  # erasing superflous descriptions
  kegg_result@result$Description <- gsub(" - Mus musculus (house mouse)", "", kegg_result@result$Description, fixed = TRUE)
  
  # save table to disk
  # stop("Coding of function is not finished yet")
  
  if (isTRUE(save_to_disk)) {
      
    message("Formatting results table")
      
    kegg_result_tibble <- as_tibble(kegg_result@result)
      
    geneID_names_list <- vector(mode = 'list', length = length(kegg_result_tibble[["geneID"]]))
      
    for (i in seq(length(kegg_result_tibble[["geneID"]]))){
      
      geneID_names_list[[i]] <- TopTableListItem[ which(TopTableListItem[["ENTREZID"]] %in% str_split(kegg_result_tibble[["geneID"]], pattern = "/")[[i]]), "SYMBOL"] %>% pull(SYMBOL)
      
      }
      
    kegg_result_tibble[["geneName"]] <- unlist(lapply(geneID_names_list,  function (x) paste(x, collapse = "/")))
      
    message("Saving table to disk.")
    
    # set table path if it isn't set
    if(is.null(table_path)){
      
      message("Table export path not provided, using hard-coded one.")
      
      table_path <- paste0(here("tables"), "/050_r_array_analysis_", "KEGG_terms_sign__", gsub(" ", "", TopTableListItemName)  , ".xlsx")
    }
      
    stopifnot("Please provide a full path for the xslx output." = !is.null(table_path))
      
    openxlsx::write.xlsx(data.frame(kegg_result_tibble), file = table_path, asTable = TRUE)
      
    }
    
    # check resulting object
    print(kegg_result)
    
    # get dipslay item
    kegg_plot <- enrichplot::dotplot(kegg_result, title =  paste0("KEGG pathways of data set: \"", TopTableListItemName,
          "\"", sep = ""
        ), showCategory = 5)
    
    # return plot
    return(kegg_plot)
    
  }

# GO plot creation
get_one_go_plot <- function(TopTableListItem, TopTableListItemName, save_to_disk = TRUE, table_path = NULL){
  
  # packages
  require("clusterProfiler")
  require("enrichplot")
  
  # for function building only
  # stop("Remove function building code")
  # TopTableListItem = FULL_TopTableListAppended[[2]]
  # TopTableListItemName = names(FULL_TopTableListAppended)[[2]]
  # save_to_disk = TRUE
  # table_path = NULL
  # rm(list(TopTableListItem, TopTableListItemName))
  
  # diagnostic
  message(paste0("Creating GO plot for data set: \"", TopTableListItemName, "\".", sep = ""))
  
 # look- up go pathways
  go_result <- enrichGO(gene = TopTableListItem$ENTREZID, keyType = "ENTREZID",  OrgDb = "org.Mm.eg.db", ont = "all")
  
  # save table to disk
  # stop("Coding of function is not finished yet")
  
  if (isTRUE(save_to_disk)) {
    message("Formatting results table")
    
    go_result_tibble <- as_tibble(go_result@result)
    
    geneID_names_list <- vector(mode = 'list', length = length(go_result_tibble[["geneID"]]))
    
    for (i in seq(length(go_result_tibble[["geneID"]]))){
      
      message(paste0("Looking up line ", i," of ", length(go_result_tibble[["geneID"]]) ," lines."))
      
      geneID_names_list[[i]] <- TopTableListItem[ which(TopTableListItem[["ENTREZID"]] %in% str_split(go_result_tibble[["geneID"]], pattern = "/")[[i]]), "SYMBOL"] %>% pull(SYMBOL)
      
    }
    
    go_result_tibble[["geneName"]] <- unlist(lapply(geneID_names_list,  function (x) paste(x, collapse = "/")))
    
    message("Saving table to disk.")
    
    # set table path if it isn't set
    if(is.null(table_path)){
      
      message("Table export path not provided, using hard-coded one.")
      
      table_path <- paste0(here("tables"), "/050_r_array_analysis_", "GO_terms_sign__", gsub(" ", "", TopTableListItemName)  , ".xlsx")
    }
    
    stopifnot("Please provide a full path for the xslx output." = !is.null(table_path))
    
    openxlsx::write.xlsx(go_result_tibble, file = table_path, asTable = TRUE)
    
  }
  
  # get display item
  go_plot <- enrichplot::dotplot(go_result, split="ONTOLOGY", title =  paste0("GO terms of data set: \"", TopTableListItemName, "\"", sep = ""), showCategory = 5) + facet_grid(ONTOLOGY ~ ., scale="free")
  
  # return plot 
  return(go_plot)
  
}

# save GO plots
save_kegg_plots <- function(ggplot_list_item, ggplot_list_item_name){
  
  require(ggplot2)
  
  filename <- gsub("[^[:alnum:][:space:]]","", ggplot_list_item_name)
  filename <- gsub("\\s+","_", filename)
  filename <- paste0("050_r_array_analysis__plot_kegg_", filename, ".pdf", collapse = "_")
  
  message(filename)
  try(ggsave(plot = ggplot_list_item, path = here("../manuscript/display_items"), filename, scale = 0.75, height = 297/2, width = 210, units = c("mm")))
  
}

# save KEGG plots
save_go_plots <- function(ggplot_list_item, ggplot_list_item_name){
  
  require(ggplot2)
  
  filename <- gsub("[^[:alnum:][:space:]]","", ggplot_list_item_name)
  filename <- gsub("\\s+","_", filename)
  filename <- paste0("050_r_array_analysis__plot_go_", filename, ".pdf", collapse = "_")
  
  message(filename)
  try(ggsave(plot = ggplot_list_item, path = here("../manuscript/display_items"), filename, scale = 1.75, height = 210, width = 210, units = c("mm")))
  
}


# Load, fomat, and shape data ----

# _1.) Full data of all 4 tissues ----

#  __a) Loading normalized data ----

# Data are Clariom S mouse arrays
# I removed strong outliers: A285_liver clustert zu bAT, A339_liver weit weg vom Rest

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/allTissues_normData.RData") # only if you are interested to look into the normalized data

# __b) Annotate data ----

normData  <- annotateEset(normData, pd.clariom.s.mouse, type = "probeset")

# __c) Copy, store, and discard data ----

# copy to stick to manuscript naming conventions
FLAT <- normData; rm(normData)

saveRDS(FLAT, file = here("rds_storage", "050_r_array_analysis__normalized_data.rds"))

# _2.) Data normalized individually for each tissue ----

# Ich habe die normalisierung für jedes Gewebe getrennt für die DGE Analysen gemacht, 
# da dies genauer ist (Gewebe zu weit auseinander in PCA)

# __a) Loading normalized data ----

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/normData4DGE.RData") #für jedes Gewebe die normalisierten Daten

# __b) Annotate data ----

bAT_normData    <- annotateEset(bAT_normData, pd.clariom.s.mouse, type = "probeset")
ingWAT_normData <- annotateEset(ingWAT_normData, pd.clariom.s.mouse, type = "probeset")
Liv_normData    <- annotateEset(Liv_normData, pd.clariom.s.mouse, type = "probeset")
bAT_normData    <- annotateEset(bAT_normData, pd.clariom.s.mouse, type = "probeset")
eWAT_normData   <- annotateEset(eWAT_normData, pd.clariom.s.mouse, type = "probeset")

# __c) Copy, store, and discard data ----

# copy to stick to manuscript naming conventions - corrected as per AH 25.05.2023, again on 07.05.2024
BRAT <- bAT_normData; rm(bAT_normData)        # brown adipose tissue
IWAT <- ingWAT_normData; rm(ingWAT_normData)  # subcutaneous adipose tissue 
LIVT <- Liv_normData; rm(Liv_normData)        # liver adipose tissue
EVAT <- eWAT_normData; rm(eWAT_normData)      # visceral adipose tissue

saveRDS(BRAT, file = here("rds_storage", "055_r_array_analysis__normalized_data_bat.rds"))
saveRDS(IWAT, file = here("rds_storage", "055_r_array_analysis__normalized_data_ewat.rds"))
saveRDS(LIVT, file = here("rds_storage", "055_r_array_analysis__normalized_data_liv.rds" ))
saveRDS(EVAT, file = here("rds_storage", "055_r_array_analysis__normalized_data_ingwat.rds"))

# _3.) Loading metadata from modelling (obesity variables) ----

# __a) Data from modelling prior to first revision ----

mice_f1_modeled_data_with_rna_seq_data <- readRDS(file = here("rds_storage", "040_r_h3__mice_f1_modeled_data_with_rna_seq_data.rds"))

# __b) Data from modelling for first submission, from {seamix} ----

mice_f0_slct_from_saemix <- readRDS(file = here("rds_storage", "mice_f0_slct_from_saemix.rds"))
mice_f1_slct_from_saemix <- readRDS(file = here("rds_storage", "mice_f1_slct_from_saemix.rds"))

# _4.) Adjust variable names and inspect data ----

# __a) Re-code experiment metadata to new nomenclature, merge in litter sizes from {saemix}-derived metadata.  --- 

mice_f1_modeled_data_with_rna_seq_data %<>% # re-code parental diet variable - slash will be buggy
  dplyr::mutate(ParentalDietMoFa = case_when( 
    ParentalDietMoFa == "chow / chow" ~ "CD CD",
    ParentalDietMoFa == "chow / HFD"  ~ "CD WD", 
    ParentalDietMoFa == "HFD / chow"  ~ "WD CD", 
    ParentalDietMoFa == "HFD / HFD"   ~ "WD WD",
    .default = as.factor(ParentalDietMoFa)
  )) %>% mutate(ParentalDietMoFa = as.factor(ParentalDietMoFa)) %>%  # re-code maternal diet variable 
  dplyr::mutate(MotherDiet = case_when( 
    MotherDiet == "HFD" ~ "WD",
    MotherDiet == "CD" ~ "CD",
    .default = as.factor(MotherDiet)
  )) %>% mutate(MotherDiet = as.factor(MotherDiet)) %>% # re-code paternal diet variable 
  dplyr::mutate(FatherDiet = case_when( 
    FatherDiet == "HFD" ~ "WD",
    FatherDiet == "CD" ~ "CD",
    .default = as.factor(FatherDiet)
  )) %>%  mutate(FatherDiet = as.factor(FatherDiet)) %>%  # re-code old analysis variable - to avoid unforeseen crashes 
  dplyr::mutate(ObeseParents = case_when( 
    ParentalDietMoFa == "CD CD" ~ "MotherFatherNotObese"      ,
    ParentalDietMoFa == "CD WD" ~ "FatherObese", 
    ParentalDietMoFa == "WD CD" ~ "MotherObese", 
    ParentalDietMoFa == "WD WD" ~ "MotherFatherObese",
    .default = as.factor(ParentalDietMoFa)
  )) %>% mutate(ObeseParents = as.factor(ObeseParents)) %>% 
  relocate(ParentalDietMoFa, .after = ObeseParents) %>% 
  dplyr::select(-Sex) %>% left_join( 
    {mice_f1_slct_from_saemix %>% dplyr::select(AnimalId, LitterSize) %>% distinct()}, by = "AnimalId"
  ) %>% arrange(DietGroup, AnimalId)

# __c) Get summary of sample sizes and treatments --- 

write_xlsx(mice_f1_modeled_data_with_rna_seq_data, 
           path =  here("../manuscript/display_items", "055_r_array_analysis_mice_f1_slct__mice_f1_modeled_data_with_rna_seq_data.xlsx")) 

# __d) Adjust array data ---

FLAT <- get_adjusted_array_data(FLAT, mice_f1_modeled_data_with_rna_seq_data)
BRAT <- get_adjusted_array_data(BRAT, mice_f1_modeled_data_with_rna_seq_data)
IWAT <- get_adjusted_array_data(IWAT, mice_f1_modeled_data_with_rna_seq_data)
LIVT <- get_adjusted_array_data(LIVT, mice_f1_modeled_data_with_rna_seq_data)
EVAT <- get_adjusted_array_data(EVAT, mice_f1_modeled_data_with_rna_seq_data)

# _5.) Save adjusted data ----

saveRDS(BRAT, file = here("rds_storage", "055_r_array_analysis__BRAT.rds"))
saveRDS(IWAT, file = here("rds_storage", "055_r_array_analysis__IWAT.rds"))
saveRDS(LIVT, file = here("rds_storage", "055_r_array_analysis__LIVT.rds"))
saveRDS(EVAT, file = here("rds_storage", "055_r_array_analysis__EVAT.rds"))

# check data, if needed, using `pData(BRAT)` and `exprs(BRAT)` - check available metadata - LIVT does not have a lot

lapply(list(FLAT, BRAT, IWAT, LIVT, EVAT), function(x) colData(x))  

# _6.) Loading AHs dietary DGE analysis results

# __a) Load data 

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/DGELists.RData") #lädt alle DGE Tabellen

# __b) Define list with dietary variables

DGEL_diet <- list(bAT_CD_HFD_VS_CD_CD, bAT_HFD_CD_VS_CD_CD, bAT_HFD_CD_VS_CD_HFD, bAT_HFD_HFD_VS_CD_CD, bAT_HFD_HFD_VS_CD_HFD,    
 bAT_HFD_HFD_VS_HFD_CD, eWAT_CD_HFD_VS_CD_CD, eWAT_HFD_CD_VS_CD_CD, eWAT_HFD_CD_VS_CD_HFD, eWAT_HFD_HFD_VS_CD_CD, 
 eWAT_HFD_HFD_VS_CD_HFD, eWAT_HFD_HFD_VS_HFD_CD, ingWAT_CD_HFD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_HFD,
 ingWAT_HFD_HFD_VS_CD_CD, ingWAT_HFD_HFD_VS_CD_HFD, ingWAT_HFD_HFD_VS_HFD_CD, Liv_CD_HFD_VS_CD_CD, Liv_HFD_CD_VS_CD_CD,
 Liv_HFD_CD_VS_CD_HFD, Liv_HFD_HFD_VS_CD_CD, Liv_HFD_HFD_VS_CD_HFD, Liv_HFD_HFD_VS_HFD_CD) 

names(DGEL_diet) <-  c("BRAT_CD_HFD_VS_CD_CD", "BRAT_HFD_CD_VS_CD_CD", "BRAT_HFD_CD_VS_CD_HFD", "BRAT_HFD_HFD_VS_CD_CD", "BRAT_HFD_HFD_VS_CD_HFD",     
 "BRAT_HFD_HFD_VS_HFD_CD", "EVAT_CD_HFD_VS_CD_CD", "EVAT_HFD_CD_VS_CD_CD", "EVAT_HFD_CD_VS_CD_HFD", "EVAT_HFD_HFD_VS_CD_CD",  
 "EVAT_HFD_HFD_VS_CD_HFD", "EVAT_HFD_HFD_VS_HFD_CD", "IWAT_CD_HFD_VS_CD_CD", "IWAT_HFD_CD_VS_CD_CD", "IWAT_HFD_CD_VS_CD_HFD", 
 "IWAT_HFD_HFD_VS_CD_CD", "IWAT_HFD_HFD_VS_CD_HFD", "IWAT_HFD_HFD_VS_HFD_CD", "LIVT_CD_HFD_VS_CD_CD", "LIVT_HFD_CD_VS_CD_CD", 
 "LIVT_HFD_CD_VS_CD_HFD", "LIVT_HFD_HFD_VS_CD_CD", "LIVT_HFD_HFD_VS_CD_HFD", "LIVT_HFD_HFD_VS_HFD_CD")

# __c) CONTINUE HERE IF RE-IMPLEMNETATION OF DGE FAILS: Define list with obesity variables 

DGEL_obes <- DGEL_diet

# As per 01.06.2023 in README commit `5fd8790e5024bce8e05f885d08f219b1c736ef58`:
# modify or exchange dietary variables to contain information regarding offsprings and/or parental obesity as per last scripts overview and manuscript tasks

# names(DGEL_obes) <-  c( "foo", "bar")

# __d) Clean environment and save data 

rm(bAT_CD_HFD_VS_CD_CD, bAT_HFD_CD_VS_CD_CD, bAT_HFD_CD_VS_CD_HFD, bAT_HFD_HFD_VS_CD_CD, bAT_HFD_HFD_VS_CD_HFD,    
   bAT_HFD_HFD_VS_HFD_CD, eWAT_CD_HFD_VS_CD_CD, eWAT_HFD_CD_VS_CD_CD, eWAT_HFD_CD_VS_CD_HFD, eWAT_HFD_HFD_VS_CD_CD, 
   eWAT_HFD_HFD_VS_CD_HFD, eWAT_HFD_HFD_VS_HFD_CD, ingWAT_CD_HFD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_CD, ingWAT_HFD_CD_VS_CD_HFD,
   ingWAT_HFD_HFD_VS_CD_CD, ingWAT_HFD_HFD_VS_CD_HFD, ingWAT_HFD_HFD_VS_HFD_CD, Liv_CD_HFD_VS_CD_CD, Liv_HFD_CD_VS_CD_CD,
   Liv_HFD_CD_VS_CD_HFD, Liv_HFD_HFD_VS_CD_CD, Liv_HFD_HFD_VS_CD_HFD, Liv_HFD_HFD_VS_HFD_CD)

saveRDS(DGEL_diet, file = here("rds_storage", "055_r_array_analysis__dge_lists_by_diet.rds"))

# Clean out expression data ----

# as described here https://nbisweden.github.io/excelerate-scRNAseq/session-clustering/Clustering.html


# _1.) Check data (abbreviate) ----

FLAT;BRAT;IWAT;LIVT;EVAT

table(colData(FLAT)$ParentalDietMoFa)
table(colData(BRAT)$ParentalDietMoFa)
table(colData(IWAT)$ParentalDietMoFa)
table(colData(LIVT)$ParentalDietMoFa)
table(colData(EVAT)$ParentalDietMoFa)

# _2.) Feature selection ----

# __a) Plot average expression and variance ----

# export plots in this section if needed
plot_mean_expression_and_variance(FLAT)
plot_mean_expression_and_variance(BRAT)
plot_mean_expression_and_variance(IWAT)
plot_mean_expression_and_variance(EVAT)

# __b) Plot rare genes ----

plot_rare_genes(FLAT)
plot_rare_genes(BRAT)
plot_rare_genes(IWAT)
plot_rare_genes(LIVT)
plot_rare_genes(EVAT)

# __c) Remove rare genes (not needed)

# omitted - code on web page

# __d) Remove genes only detected in few samples (not needed)

# omitted - code from web page
# se_object <- BRAT
# samples_per_gene <- nexprs(se_object, assay.type = "exprs",  byrow = TRUE) 
# samples_per_gene_subest <- samples_per_gene >= 3
# se_object <- se_object[samples_per_gene_subest,]
# dim(se_object)

# __e) Isolating highly variable genes (not needed)

# omitted - code from web page
# se_object <- EVAT
# dec <- modelGeneVar(se_object, assay.type = "exprs" )
# dec$HVG <- (dec$FDR < 0.00001)
# plot(dec$mean, dec$total, pch = 16, cex = 0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
# o <- order(dec$mean)
# lines(dec$mean[o], dec$tech[o], col="dodgerblue", lwd=2)
# points(dec$mean[dec$HVG], dec$total[dec$HVG], col="red", pch=16)

# Principal Component Analysis ----

# as described here https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
# see file:///Users/paul/Documents/HM_miscellaneous/231116_cf-statcon_mv_statistics/Exercises/Part_1/Exercise_3.html#Supplementary_individuals_and_weights

# _1.) Get PCs ----

# __a) R-mode ----

PCA_FLAT <- get_principal_components(FLAT)
PCA_BRAT <- get_principal_components(BRAT)
PCA_IWAT <- get_principal_components(IWAT)
PCA_LIVT <- get_principal_components(LIVT)
PCA_EVAT <- get_principal_components(EVAT)

# __b) Q-mode ----

# As per Camargo, Arley. 2022. "PCAtest: Testing the Statistical Significance of
# Principal Component Analysis in R." PeerJ 10 (February):e12967.
# https://doi.org/10.7717/peerj.12967. Get q-mode PCAs if explained variation on
# the PC1 and PC2 is low

qPCA_FLAT <- get_q_mode_principal_components(FLAT)
qPCA_BRAT <- get_q_mode_principal_components(BRAT)
qPCA_IWAT <- get_q_mode_principal_components(IWAT)
qPCA_LIVT <- get_q_mode_principal_components(LIVT)
qPCA_EVAT <- get_q_mode_principal_components(EVAT)

# _2.) Get PC loads ----

# __a) R-mode ----

PV_FLAT <- get_percent_variation(PCA_FLAT)
PV_BRAT <- get_percent_variation(PCA_BRAT)
PV_IWAT <- get_percent_variation(PCA_IWAT)
PV_LIVT <- get_percent_variation(PCA_LIVT)
PV_EVAT <- get_percent_variation(PCA_EVAT)

# __b) Q-mode ----

qPV_FLAT <- get_percent_variation(qPCA_FLAT)
qPV_BRAT <- get_percent_variation(qPCA_BRAT)
qPV_IWAT <- get_percent_variation(qPCA_IWAT)
qPV_LIVT <- get_percent_variation(qPCA_LIVT)
qPV_EVAT <- get_percent_variation(qPCA_EVAT)

# _3.) Plot Principal Component Analyses ----

# __a.) Plot FLAT PCs (R mode - enough variation)----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_flat_a <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "Tissue", legend_title = "F1 Tissue", plot_title =  "a", percent_var = PV_FLAT)

# Overall expression differences and obesity status among f1 offspring
plot_pca_flat_b <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "b", percent_var = PV_FLAT)

# Overall expression differences and obesity parental obesity status among f0
plot_pca_flat_c <- get_pca_plot(expr_data_pca = PCA_FLAT, expr_data_raw = FLAT , variable = "LitterSize", legend_title = "F1 litter size", plot_title =  "c", percent_var = PV_FLAT)

# Combine and save plots

# factor labeling needs to change - check above
plot_pca_flat <- ggarrange(plot_pca_flat_a, plot_pca_flat_b, plot_pca_flat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_flat, path = here("plots"), 
       filename = "055_r_array_analysis__plot_pca_flat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_flat, path = here("../manuscript/display_items"), 
       filename = "055_r_array_analysis__plot_pca_flat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# __b.) Plot BRAT PCs  ----

# Overall expression differences and obesity status among f1 offspring

# plot_pca_brat_a <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "a", percent_var = PV_BRAT)
# see `file:///Users/paul/Documents/HM_miscellaneous/231116_cf-statcon_mv_statistics/Exercises/Part_1/Exercise_3.html#Principle_Component_Analysis` for interpretation

plot_pca_brat_a <- get_q_mode_pca_plot(expr_data_pca = qPCA_BRAT, expr_data_raw = BRAT, variable = "ParentalDietMoFa", plot_title =  "a", percent_var = qPV_FLAT)

# Overall expression differences and obesity parental obesity status among f0

# plot_pca_brat_b <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT, variable = "LitterSize", legend_title = "F1 litter size", plot_title =  "b", percent_var = PV_BRAT)
plot_pca_brat_b <- get_q_mode_pca_plot(expr_data_pca = qPCA_BRAT, expr_data_raw = BRAT, variable = "LitterSize", plot_title =  "b", percent_var = qPV_FLAT)

# Combine and save plots
plot_pca_brat <- ggarrange(plot_pca_brat_a, plot_pca_brat_b, nrow = 1, ncol = 2)

ggsave(plot = plot_pca_brat, path = here("plots"), 
       filename = "055_r_array_analysis__plot_pca_brat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_brat, path = here("../manuscript/display_items"), 
       filename = "055_r_array_analysis__plot_pca_brat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# __c.) Plot IWAT PCs----

# Overall expression differences and obesity status among f1 offspring

# plot_pca_scat_a <- get_pca_plot(expr_data_pca = PCA_IWAT, expr_data_raw = IWAT , variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "a", percent_var = PV_IWAT)
plot_pca_scat_a <- get_q_mode_pca_plot(expr_data_pca = qPCA_IWAT, expr_data_raw = IWAT, variable = "ParentalDietMoFa", plot_title =  "a", percent_var = qPV_IWAT)

# Overall expression differences and obesity parental obesity status among f0

# plot_pca_scat_b <- get_pca_plot(expr_data_pca = PCA_IWAT, expr_data_raw = IWAT , variable = "LitterSize", legend_title = "F1 litter size", plot_title =  "b", percent_var = PV_IWAT)
plot_pca_scat_b <- get_q_mode_pca_plot(expr_data_pca = qPCA_IWAT, expr_data_raw = IWAT, variable = "LitterSize", plot_title =  "b", percent_var = qPV_IWAT)

# Combine and save plots

plot_pca_scat <- ggarrange(plot_pca_scat_a, plot_pca_scat_b, nrow = 1, ncol = 2)

ggsave(plot = plot_pca_scat, path = here("plots"), 
       filename = "055_r_array_analysis__plot_pca_scat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_scat, path = here("../manuscript/display_items"), 
       filename = "055_r_array_analysis__plot_pca_scat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# __d.) Plot EVAT PCs ----

# Overall expression differences and obesity status among f1 offspring

# plot_pca_evat_a <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "a", percent_var = PV_EVAT)
plot_pca_evat_a <- get_q_mode_pca_plot(expr_data_pca = qPCA_EVAT, expr_data_raw = EVAT , variable = "ParentalDietMoFa", percent_var = qPV_EVAT, plot_title =  "a")

# Overall expression differences and obesity parental obesity status among f0

# plot_pca_evat_b <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT , variable = "LitterSize", legend_title = "F1 litter size", plot_title =  "b", percent_var = PV_EVAT)
plot_pca_evat_b <- get_q_mode_pca_plot(expr_data_pca = qPCA_EVAT, expr_data_raw = EVAT , variable = "LitterSize", percent_var = qPV_EVAT, plot_title =  "b")

# Combine and save plots
plot_pca_evat <- ggarrange(plot_pca_evat_a, plot_pca_evat_b, nrow = 1, ncol = 2)

ggsave(plot = plot_pca_evat, path = here("plots"), 
       filename = "055_r_array_analysis__plot_pca_evat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_evat, path = here("../manuscript/display_items"), 
       filename = "055_r_array_analysis__plot_pca_evat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# __e.) Plot LIVT PCs ----

# Overall expression differences and obesity status among f1 offspring

# plot_pca_liat_a <- get_pca_plot(expr_data_pca = PCA_LIVT, expr_data_raw = LIVT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "a", percent_var = PV_LIVT)
plot_pca_liat_a <- get_q_mode_pca_plot(expr_data_pca = qPCA_LIVT, expr_data_raw = LIVT, variable = "ParentalDietMoFa", plot_title =  "a", percent_var = qPV_LIVT)

# Overall expression differences and obesity parental obesity status among f0

# plot_pca_liat_b <- get_pca_plot(expr_data_pca = PCA_LIVT, expr_data_raw = LIVT, variable = "LitterSize", legend_title = "F1 litter size", plot_title =  "b", percent_var = PV_LIVT)
plot_pca_liat_b <- get_q_mode_pca_plot(expr_data_pca = qPCA_LIVT, expr_data_raw = LIVT, variable = "LitterSize", plot_title =  "a", percent_var = qPV_LIVT)

# Combine and save plots

plot_pca_liat <- ggarrange(plot_pca_liat_a, plot_pca_liat_b, nrow = 1, ncol = 2)

ggsave(plot = plot_pca_liat, path = here("plots"), 
       filename = "055_r_array_analysis__plot_pca_liat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_liat, path = here("../manuscript/display_items"), 
       filename = "055_r_array_analysis__plot_pca_liat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# __f.) Plot for talks (in R mode) ----

plot_pca_brat_one <- get_pca_plot(expr_data_pca = PCA_BRAT, expr_data_raw = BRAT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "BAT: interscapular brown AT ", percent_var = PV_BRAT)
plot_pca_evat_two <- get_pca_plot(expr_data_pca = PCA_EVAT, expr_data_raw = EVAT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "EVAT: epigonal visceral AT", percent_var = PV_EVAT)
plot_pca_livr_thr <- get_pca_plot(expr_data_pca = PCA_LIVT, expr_data_raw = LIVT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "L: liver T", percent_var = PV_LIVT)
plot_pca_scat_for <- get_pca_plot(expr_data_pca = PCA_IWAT, expr_data_raw = IWAT, variable = "ParentalDietMoFa", legend_title = "F1 parental diet\n(mother / father)", plot_title =  "IWAT: inguinal subcutaneous AT", percent_var = PV_IWAT)

plot_pca_summ <- ggarrange(plot_pca_brat_one,
                           plot_pca_evat_two,
                           plot_pca_livr_thr,
                           plot_pca_scat_for, nrow = 2, ncol = 2)
ggsave(plot = plot_pca_summ, path = here("plots"), 
       filename = "055_r_array_analysis__plot_pca_summ.pdf",  
       width = 500, height = 275, units = "mm", dpi = 300,  limitsize = TRUE, scale = 0.75)

# >>> Code construction in progress - implementing new analysis flow here. ----

warning("Code construction in progress.")

# Show that obesity-related genes are in the four tissues ----

# Will be  used to warrant DEG search, PCA analysis will move to supplement. AS
# per 10.3390/ijms231911005 using: Leptin (LEP), the leptin receptor (LEPR),
# proopiomelanocortin (POMC), prohormone convertase 1 (PCSK1), the melanocortin
# 4 receptor (MC4R), single-minded homolog 1 (SIM1), brain-derived neurotrophic
# factor (BDNF), and the neurotrophic tyrosine kinase receptor type 2 gene
# (NTRK2)

# 1.) Compile data ----

SE_list <- list(BRAT, IWAT, LIVT, EVAT)
names(SE_list) <- c("BRAT", "IWAT", "LIVT", "EVAT")
obesity_genes <- c("LEP", "LEPR", "POMC", "PCSK1", "MC4R", "SIM1", "BDNF", "NTRK2") # from 10.3390/ijms231911005

# 2.) Filter data to obesity genes and check ----

SE_list_og <- lapply(SE_list, function(seob) seob[   which( rowData(seob)[["SYMBOL"]] %in% obesity_genes ) ] )
lapply(SE_list_og,  function(seob) toupper(rowData(seob)[["SYMBOL"]]))

# 3.) Show obesity genes ----

# https://github.com/plger/sechm
# https://www.bioconductor.org/packages/release/bioc/vignettes/sechm/inst/doc/sechm.html

get_one_heatmap = function(se_ob, gaps_at = "ParentalDietMoFa", mark = c("LEP", "LEPR", "POMC", "PCSK1", "MC4R", "SIM1", "BDNF", "NTRK2")){
  
  # for function building
  se_ob <-  SE_list_og[[1]]
  
  # build_heatmap
  sechm_heatmap <- sechm(
    se_ob,
    features = rownames(se_ob),
    do.scale = TRUE,
    gaps_at = gaps_at,
    mark = mark,
    show_rownames = TRUE)
  # modify row lables
  
  rownames(foo@matrix)
  
  
  
  
  
}







# Get Fig. 2

# Look for DEGs for all contrasts in each of the 4 tissues. ---- 

# Get Suppl. Tables 1-6.  Look for DEGs for 3 (and all) contrasts (CD CD against WD CD, CD WD, and WD WD) in each of the 4 tissues (results in 12 lists).

# Decsribe  besity related genes in intersections of full DEG lists ----

# Get Fig. 3: Upset plot of intersections of full list for all contrasts

# Subest DEG  lists to set differences among all contrasts in each of the 4 tissues. ----

# Get Suppl. Tables 6-12.  Look for DEGs for 3 (and all) contrasts (CD CD against WD CD, CD WD, and WD WD) in each of the 4 tissues (results in 12 lists).

# Decsribe  besity related genes in intersections of set different DEG lists ----

#  Get Fig. 4: Upset plot of intersections of full list for all contrasts

# GSEA of set differencess with unique obesity genes

# Get Suppl. Tab 12-n 

# Save Environment ----



stop("Unadjusted old analysis code below - possibly integrate into code above. ")


# >>> Unadjusted old analysis code below - possibly integrate into code above. ----





# Re-implement analysis of array intensities ----

#' ## Shape and check array intensity data

# _1.) Shape and check array intensity data ----

# see https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html

#' ###  Check input data formats

# __a) Check input data formats ----

# see available expression data 
FLAT; BRAT; IWAT; LIVT; EVAT

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

FLAT_DT.m1 <- melt(FLAT_DT,  id.vars = c("Sample", "Animal", "Tissue", "AnimalSex", "ObesityLgcl", "ObeseParents", "MotherDiet", "FatherDiet", "AnimalSex", "ParentalDietMoFa", "DietGroup"), 
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
warning("get_dge_for_offspring_obesity() has not been adjusted for 1st revisions because it is likely unenecssary.")

# get_dge_for_offspring_obesity(FLAT)
# get_dge_for_offspring_obesity(BRAT)
# get_dge_for_offspring_obesity(LIVT)
# get_dge_for_offspring_obesity(IWAT)
# get_dge_for_offspring_obesity(EVAT)

#' ### Test for DGE among offspring based on parental obesity

# __c)  Test for DGE among offspring based on parental obesity ----

stop("Update code based on methods in main text, on or after 22.05.2024, then update results")
stop("LFC needs tp be set to 2 for most contarsts, see main text.")

# Defining and applying contrasts: One of "MotherFatherNotObese", "FatherObese", "MotherFatherObese", or "MotherObese" 
#  against all remaining three levels.
#  Compare to PCA results, including coefficients - Expression is variable based on tissue, and within each tissue based on parental obesity

FLAT_TopTableList <- get_dge_for_parent_obesity(FLAT) # analysis off all tissues together doesn't really make doesn't really make sense
EVAT_TopTableList <- get_dge_for_parent_obesity(EVAT) # likely not needed - see manuscript results 05.07.2023
SCAT_TopTableList <- get_dge_for_parent_obesity(IWAT) # likely not needed - see manuscript results 05.07.2023
BRAT_TopTableList <- get_dge_for_parent_obesity(BRAT) 
LIAT_TopTableList <- get_some_dge_for_parent_obesity(LIVT) 

#' ###  Choose DGE results for further analyses

# __d) Choose DGE results for further analyses ----

# ___ FLAT - keeping all contrasts for all tissues ----

# see file "/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__text_pca_FLAT.txt"

names(FLAT_Tissue_TopTableList)
          
# ___ FLAT, BRAT, IWAT, LIVT, EVAT - not keeping any contrasts defined by offspring' obesity  ----

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

# ___ IWAT - keeping some contrasts defined by parents' obesity  ----

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

# ___ LIVT - keeping some contrasts defined by parents' obesity  ----

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
#   names(SCAT__Select_TopTableList) <- paste0("IWAT - ", names(SCAT__Select_TopTableList))

names(BRAT__Select_TopTableList) <- paste0("BRAT - ", names(BRAT__Select_TopTableList))
names(LIAT__Select_TopTableList) <- paste0("LIVT - ", names(LIAT__Select_TopTableList))

FULL_TopTableList <- c(
  # FLAT_Tissue_TopTableList,
  # EVAT__Select_TopTableList, 
  BRAT__Select_TopTableList,
  LIAT__Select_TopTableList
  )

# receiving a table with 2 slots, each containing DGE results fo a specific tissue and statistically relavent contrasts
names(FULL_TopTableList)

# __f) Filter top tables further if required ----   

# check tables
FULL_TopTableList[[1]]
FULL_TopTableList[[2]]

# check for NAs
FULL_TopTableList[[1]][rowSums(is.na(FULL_TopTableList[[1]])) > 0, ]
FULL_TopTableList[[2]][rowSums(is.na(FULL_TopTableList[[2]])) > 0, ]

# Inspect data in question - green tail on Volcano plot
FULL_TopTableList[[1]][ which(  log2(FULL_TopTableList[[1]]$logFC) > 1 ) , ]
FULL_TopTableList[[1]][ which(  log2(FULL_TopTableList[[1]]$logFC) > 1 ) , ]


#' ### Get, assort, arrange, and save Vulcano plots

# __g)  Get, assort, arrange, and save Vulcano plots  ----

FULL_VolcanoPlots <- mapply(get_one_volcanoplot, TopTableListItem = FULL_TopTableList, TopTableListItemName = gsub("LIVT", "L", gsub("BRAT", "BAT", names(FULL_TopTableList),  fixed = TRUE), fixed = TRUE), SIMPLIFY = FALSE)

# the following three sets are undefined or serve no purpose (see 5-Jul-2023)
#  FLAT_VolcanoPlots <- FULL_VolcanoPlots[grep("FLAT", names(FULL_VolcanoPlots))]
#  SCAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^IWAT", names(FULL_VolcanoPlots))]
#  EVAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^EVAT", names(FULL_VolcanoPlots))]

BRAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^BRAT", names(FULL_VolcanoPlots))]
LIAT_VolcanoPlots <- FULL_VolcanoPlots[grep("^LIVT", names(FULL_VolcanoPlots))]

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

ggsave(plot = BRAT_VolcanoPlotsComposite, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_volcano_brat.pdf",  
       width = 180, height = 200, units = "mm", dpi = 300,  limitsize = TRUE, scale = 1)

ggsave(plot = LIAT_VolcanoPlotsComposite, path = here("../manuscript/display_items"), 
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
foo_mat <- base::as.matrix(foo %>% dplyr::select(-PROBEID))
rownames(foo_mat) <- (foo %>% pull(PROBEID))

# isolate the DEG top table
bar <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]]

# join DEG top table and expression data to a new tibble  
foobar <- left_join(bar, foo)

# keep only relavent logFC, conforming with Volcano plots, p values do not need further adjustments
foobar <- foobar %>% filter(logFC < -1 | logFC > 1) %>% arrange(logFC)
foobar <- foobar %>% filter(adj.P.Val < 0.05 )%>% arrange(logFC)

# pheatmap needs a matrix - hence converting tibble to matrix
foo_mat <- base::as.matrix (foobar %>% dplyr::select(contains("_BRAT")))
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

# renaming tissues
colnames(foo_mat) <- gsub("LIVT", "L", gsub("BRAT", "BAT", colnames(foo_mat),  fixed = TRUE), fixed = TRUE)

# print heat map to script and file
pheatmap(foo_mat, scale = "row")
pheatmap(foo_mat, scale = "row", filename =  paste0(here("../manuscript/display_items"),"/","050_r_array_analysis__plot_heatmap_brat.pdf"))

# ___ LIVT ----

# get expression data as matrix with probe ids, subset for speed 
foo <- as_tibble(exprs(LIVT), rownames = c("PROBEID")) %>% filter(PROBEID %in%  LIAT__Select_TopTableList[[1]][["PROBEID"]])
foo_mat <- base::as.matrix(foo %>% dplyr::select(-PROBEID))
rownames(foo_mat) <- (foo %>% pull(PROBEID))

# isolate the DEG top table
bar <- LIAT__Select_TopTableList[["LIVT - MotherFatherObese vs MotherFatherNotObese"]]

# join DEG top table and expression data to a new tibble  
foobar <- left_join(bar, foo)

# keep only relavent logFC and pValues
foobar <- foobar %>% filter(logFC < -1 | logFC > 1) %>% arrange(logFC)
foobar <- foobar %>% filter(adj.P.Val < 0.05 ) %>% arrange(logFC)

# pheatmap needs a matrix - hence converting tibble to matrix
foo_mat <- base::as.matrix (foobar %>% dplyr::select(contains("_LIAT")))
rownames(foo_mat) <- foobar[["SYMBOL"]]

# to add more sample annotations to heat map matrix colnames - isolating this data here
mdata_tibble <- tibble(pData(LIVT)) %>% mutate(SAMPLE = paste0(Animal,"_",Tissue))

# checking possibilty to extend matrix column names
which(colnames(foo_mat) %in% mdata_tibble$SAMPLE)

# extending matrix column names (samples) with relevant metadata
colnames(foo_mat) <- paste0(colnames(foo_mat), " | " , 
                            mdata_tibble[which(colnames(foo_mat) %in% mdata_tibble$SAMPLE) , "ObeseParents"][["ObeseParents"]], " | " ,
                            mdata_tibble[which(colnames(foo_mat) %in% mdata_tibble$SAMPLE) , "ObesityLgcl"][["ObesityLgcl"]]
)

# adding stars to contrast
colnames(foo_mat)[grep(pattern = "MotherFatherObese|MotherFatherNotObese", colnames(foo_mat))] <- paste(colnames(foo_mat)[grep(pattern = "MotherFatherObese|MotherFatherNotObese", colnames(foo_mat))], "*")

# renaming tissues
colnames(foo_mat) <- gsub("LIVT", "L", gsub("BRAT", "BAT", colnames(foo_mat),  fixed = TRUE), fixed = TRUE)

# print heat map to script and file
pheatmap(foo_mat, scale = "row")
pheatmap(foo_mat, scale = "row", filename =  paste0(here("../manuscript/display_items"),"/","050_r_array_analysis__plot_heatmap_liat.pdf"))

#' ### Save DGE lists

# __i) Save DGE lists ----

BRAT_TTL_sign <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]] %>% filter(adj.P.Val < 0.05 )
BRAT_TTL_uprg <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]] %>% filter(logFC > 1)
BRAT_TTL_down <- BRAT__Select_TopTableList[["BRAT - MotherFatherObese vs FatherObese"]] %>% filter(logFC < -1)

LIAT_TTL_sign <- LIAT__Select_TopTableList[["LIVT - MotherFatherObese vs MotherFatherNotObese"]] %>% filter(adj.P.Val < 0.05 )
LIAT_TTL_uprg <- LIAT__Select_TopTableList[["LIVT - MotherFatherObese vs MotherFatherNotObese"]] %>% filter(logFC > 1)
LIAT_TTL_down <- LIAT__Select_TopTableList[["LIVT - MotherFatherObese vs MotherFatherNotObese"]] %>% filter(logFC < -1)

openxlsx::write.xlsx(BRAT_TTL_sign, paste0(here("tables"), "/050_r_array_analysis_", "BAT_DEGs_sign", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(BRAT_TTL_down, paste0(here("tables"), "/050_r_array_analysis_", "BAT_DEGs_down", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(BRAT_TTL_uprg, paste0(here("tables"), "/050_r_array_analysis_", "BAT_DEGs_uprg", ".xlsx"), asTable = TRUE, overwrite = TRUE)

openxlsx::write.xlsx(LIAT_TTL_sign, paste0(here("tables"), "/050_r_array_analysis_", "L_DEGs_sign", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(LIAT_TTL_down, paste0(here("tables"), "/050_r_array_analysis_", "L_DEGs_down", ".xlsx"), asTable = TRUE, overwrite = TRUE)
openxlsx::write.xlsx(LIAT_TTL_uprg, paste0(here("tables"), "/050_r_array_analysis_", "L_DEGs_uprg", ".xlsx"), asTable = TRUE, overwrite = TRUE)

#' ## Gene Set Enrichment Analysis (GSEA)

# _3.) Gene Set Enrichment Analysis (GSEA) ----

# the following resources are useful 
# step-by-step GSEA https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/08_practical_DE.pdf
# step-by-step GSEA https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# why shrinkage of LFCs is not necessary, or already done https://support.bioconductor.org/p/100804/

# check test data - possibly used tissue-specific data sets

#' ### Lookup Entrez ID for Clusterprofiler

# __a) Lookup Entrez ID for Clusterprofiler ----

FULL_TopTableListAppended <- lapply(FULL_TopTableList, get_entrez_ids) 

# __b) Adjust data for KEGG and GO analysis ----

# Keep only relavent logFC and pValues Clusterprofiler ----
#  not needed - default values for `enrichGo()` and `enrichKegg()` are `pvalueCutoff = 0.05, pAdjustMethod = "BH"` # "holm"

# FULL_TopTableListAppended <- lapply(FULL_TopTableListAppended, function (li) li %<>% filter(logFC < -1 | logFC > 1) %>% arrange(logFC))
# FULL_TopTableListAppended <- lapply(FULL_TopTableListAppended, function (li) li %<>% filter(adj.P.Val < 0.05 ) %>% arrange(logFC))

# Adjust names
names(FULL_TopTableListAppended) <- gsub("BRAT", "BAT", names(FULL_TopTableListAppended), fixed = TRUE )
names(FULL_TopTableListAppended) <- gsub("LIVT", "L", names(FULL_TopTableListAppended), fixed = TRUE )

#' ### Get plots of KEEG pathways

# __c) Get plots of KEEG pathways ----

FULL_KeggPlots <- mapply(get_one_kegg_dotplot, TopTableListItem = FULL_TopTableListAppended, TopTableListItemName = names(FULL_TopTableListAppended), SIMPLIFY = FALSE)
names(FULL_KeggPlots)

#' ### Save plots of KEEG pathways

# __d) Save plots of KEEG pathways ----

FULL_KeggPlots[1]
FULL_KeggPlots[2]

mapply(save_kegg_plots, ggplot_list_item = FULL_KeggPlots, ggplot_list_item_name = names(FULL_KeggPlots), SIMPLIFY = FALSE)

ggsave(plot = ggarrange(plotlist =  FULL_KeggPlots, ncol = 2, labels = "auto"), filename = "050_r_array_analysis__plot_kegg_both.pdf", path = here("../manuscript/display_items/"), width = 280, height = 420, unit = "mm", scale = 1)

#' ### Implement GO analysis

# __e) Implement GO analysis ----

FULL_GoPlots <- mapply(get_one_go_plot, TopTableListItem = FULL_TopTableListAppended, TopTableListItemName = names(FULL_TopTableListAppended), SIMPLIFY = FALSE)
names(FULL_GoPlots)

#' ### Show and save plots of GO pathways

# __f) Show and save plots of GO pathways ----

FULL_GoPlots[[1]]
FULL_GoPlots[[2]]

mapply(save_go_plots, ggplot_list_item = FULL_GoPlots, ggplot_list_item_name = names(FULL_KeggPlots), SIMPLIFY = FALSE)

ggsave(plot = ggarrange(plotlist =  FULL_GoPlots, ncol = 2, labels = "auto"), filename = "050_r_array_analysis__plot_go_both.pdf", path = here("../manuscript/display_items/"), width = 280, height = 420, unit = "mm", scale = 1)

# #' Experimental: DGE-analysis using GAMs - Investigate overall tissue specific expression differences based on obesity variables

# Experimental: DGE-analysis using GAMs - Investigate overall tissue specific expression differences based on obesity variables ----

# from inspection 
# - check EVAT of MotherFatherNotObese across all genes
# - check LIVT and MotherFatherObese across all genes
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

