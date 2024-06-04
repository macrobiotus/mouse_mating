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

# Get a "complex heatmap" of a reasonably small summarized experiment object. 
get_one_heatmap = function(se_ob, gaps_at = "ParentalDietMoFa", mark = c("LEP", "LEPR", "POMC", "PCSK1", "MC4R", "SIM1", "BDNF", "NTRK2")){
  
  # for function building
  # se_ob <-  SE_list_og[[1]]
  # gaps_at = "ParentalDietMoFa"
  # mark = c("LEP", "LEPR", "POMC", "PCSK1", "MC4R", "SIM1", "BDNF", "NTRK2")
  
  # modify SE object 
  # - setting rwonames which can be understood
  warning("Overwriting gene names, do not use resulting obejcts outside function.") 
  rownames(se_ob) <- rowData(se_ob)[["SYMBOL"]][which(rownames(rowData(se_ob)) %in% rownames(se_ob))] 
  # - setting colour which can be understood
  metadata(se_ob)$hmcols <- c("#0072B2","white","#D55E00")
  
  # build_heatmap
  sechm_heatmap <- sechm(
    se_ob,
    features = rownames(se_ob),
    do.scale = TRUE,
    gaps_at = gaps_at,
    # mark = mark,
    show_rownames = TRUE)
  
  return(sechm_heatmap)
  
}

# Calculate DEGs - factor is currently hardcoded.
get_deg_lists = function(se_ob, peval = 0.05, logfc = 2.0){
  
  require("limma")
  require("SummarizedExperiment")
  
  # initialize for function development
  # stop("Remove function building code")
  # se_ob <- SE_list[[1]]
  # treat <- "ParentalDietMoFa"
  # peval <- 0.05
  # logfc <- 2.0
  
  treat = "ParentalDietMoFa"
  
  # Abort if wrong object is passed in and if treatment isn't available
  stopifnot(class(se_ob) == "SummarizedExperiment")
  stopifnot(!is.null(colData(se_ob)[[treat]]))
  
  # Diagnostic messages
  message(paste0("\nUsing data set \"", deparse(substitute(se_ob))), "\".")
  message(paste0("Analysis of factor \"", treat, "\" is hard-coded, to not soft-code contrast definition."))
  
  # Creating design matrix
  model_matrix <- model.matrix(~ se_ob[[treat]] - 1)
  colnames(model_matrix) <- make.names(levels(se_ob[[treat]]))
  
  # Fit linear model for each gene given a series of arrays
  fit_limma <- lmFit(assay(se_ob), model_matrix)
  
  # Defining contrasts
  contrast_list <- vector(mode = "list", length = 6)
  names(contrast_list) <- c("CD CD - WD WD", "CD CD - CD WD", "CD CD - WD CD", "WD WD - CD WD", "WD WD - WD CD", "WD CD - CD WD")
  contrast_list[[1]]  <- makeContrasts("CD CD - WD WD" =  CD.CD - WD.WD, levels = model_matrix)
  contrast_list[[2]]  <- makeContrasts("CD CD - CD WD" =  CD.CD - CD.WD, levels = model_matrix)
  contrast_list[[3]]  <- makeContrasts("CD CD - WD CD" =  CD.CD - WD.CD, levels = model_matrix)
  contrast_list[[4]]  <- makeContrasts("WD WD - CD WD" =  WD.WD - CD.WD, levels = model_matrix)
  contrast_list[[5]]  <- makeContrasts("WD WD - WD CD" =  WD.WD - WD.CD, levels = model_matrix)
  contrast_list[[6]]  <- makeContrasts("WD CD - CD WD" =  WD.CD - CD.WD, levels = model_matrix)
  
  # Applying contrasts, and empirical Bayes moderation
  fit_limma_contrast_list <- lapply(contrast_list, function (x) contrasts.fit(fit_limma, x))
  fit_limma_contrast_list_ebayes <- lapply(fit_limma_contrast_list, function (x) eBayes(x))
  fit_limma_contrast_list_ebayes_toptable <- lapply(fit_limma_contrast_list_ebayes, function (x) topTable(x, p.value = peval, lfc = logfc, number = Inf))
  names(fit_limma_contrast_list_ebayes_toptable) <- names(contrast_list)
  
  # Adding gene symbols to nested lists
  fit_limma_contrast_list_ebayes_toptable <- lapply(fit_limma_contrast_list_ebayes_toptable, function(top_table)
    as.data.frame( merge(top_table, rowData(se_ob)[c("SYMBOL", "GENENAME")][which(rownames(rowData(se_ob)) %in% rownames(top_table)), ], by = "row.names")))
  
  # Storing list in in Summarized experiment objects 
  metadata(se_ob)$toptable_list <- fit_limma_contrast_list_ebayes_toptable  
  
  message(paste0("Returning list of top tables with log fold-chnage ", logfc , " and adjusted p-value below ", peval , " in Summarized Experiment object"))
  return(se_ob)
}

# Flatten DEG lists
get_flat_deg_lists = function(se_ob){
  
  # stop("Remove function building code")
  # se_ob <- SE_list[[1]]
  
  # isolate top table list
  toptable_list <- metadata(se_ob)$toptable_list
  # isolate tissue name
  tissue <- unique(colData(se_ob)[["Tissue"]])
  
  # basic input tests
  # pluck depth should be one - otherwise function will crash
  stopifnot(purrr:::pluck_depth(toptable_list) == 3 )
  # there should only be one tissue tipe in the input object
  stopifnot(length(tissue) == 1)
  
  # build tibble from tiblle
  flat_list <- as_tibble(bind_rows(toptable_list, .id = "Contrasts"))
  flat_list[["Tissue"]] = unique(colData(se_ob)[["Tissue"]])
  
  # format tiblle
  names(flat_list) <- base::toupper(names(flat_list)) 
  
  flat_list %<>% relocate(TISSUE, .before = "CONTRASTS") %>% dplyr::rename(PROBEID = ROW.NAMES) %>% relocate(SYMBOL, .before = PROBEID)
  
  return(flat_list)
  
}

# All functions below are unrevised ----

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

#  __a) Loading normalized data

# Data are Clariom S mouse arrays
# I removed strong outliers: A285_liver clustert zu bAT, A339_liver weit weg vom Rest

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/allTissues_normData.RData") # only if you are interested to look into the normalized data

# __b) Annotate data

normData  <- annotateEset(normData, pd.clariom.s.mouse, type = "probeset")

# __c) Copy, store, and discard data

# copy to stick to manuscript naming conventions
FLAT <- normData; rm(normData)

saveRDS(FLAT, file = here("rds_storage", "050_r_array_analysis__normalized_data.rds"))

# _2.) Data normalized individually for each tissue ----

# Ich habe die normalisierung für jedes Gewebe getrennt für die DGE Analysen gemacht, 
# da dies genauer ist (Gewebe zu weit auseinander in PCA)

# __a) Loading normalized data

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/normData4DGE.RData") #für jedes Gewebe die normalisierten Daten

# __b) Annotate data

bAT_normData    <- annotateEset(bAT_normData, pd.clariom.s.mouse, type = "probeset")
ingWAT_normData <- annotateEset(ingWAT_normData, pd.clariom.s.mouse, type = "probeset")
Liv_normData    <- annotateEset(Liv_normData, pd.clariom.s.mouse, type = "probeset")
bAT_normData    <- annotateEset(bAT_normData, pd.clariom.s.mouse, type = "probeset")
eWAT_normData   <- annotateEset(eWAT_normData, pd.clariom.s.mouse, type = "probeset")

# __c) Copy, store, and discard data

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

# __a) Data from modelling prior to first revision

mice_f1_modeled_data_with_rna_seq_data <- readRDS(file = here("rds_storage", "040_r_h3__mice_f1_modeled_data_with_rna_seq_data.rds"))

# __b) Data from modelling for first submission, from {seamix}

mice_f0_slct_from_saemix <- readRDS(file = here("rds_storage", "mice_f0_slct_from_saemix.rds"))
mice_f1_slct_from_saemix <- readRDS(file = here("rds_storage", "mice_f1_slct_from_saemix.rds"))

# _4.) Adjust variable names and inspect data ----

# __a) Re-code experiment metadata to new nomenclature, merge in litter sizes from {saemix}-derived metadata.

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

# __c) Get summary of sample sizes and treatments

write_xlsx(mice_f1_modeled_data_with_rna_seq_data, 
           path =  here("../manuscript/display_items", "055_r_array_analysis_mice_f1_slct__mice_f1_modeled_data_with_rna_seq_data.xlsx")) 

# __d) Adjust array data

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

# __a) Plot average expression and variance

# export plots in this section if needed
plot_mean_expression_and_variance(FLAT)
plot_mean_expression_and_variance(BRAT)
plot_mean_expression_and_variance(IWAT)
plot_mean_expression_and_variance(EVAT)

# __b) Plot rare genes

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

# __a) R-mode

PCA_FLAT <- get_principal_components(FLAT)
PCA_BRAT <- get_principal_components(BRAT)
PCA_IWAT <- get_principal_components(IWAT)
PCA_LIVT <- get_principal_components(LIVT)
PCA_EVAT <- get_principal_components(EVAT)

# __b) Q-mode

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

# __a) R-mode

PV_FLAT <- get_percent_variation(PCA_FLAT)
PV_BRAT <- get_percent_variation(PCA_BRAT)
PV_IWAT <- get_percent_variation(PCA_IWAT)
PV_LIVT <- get_percent_variation(PCA_LIVT)
PV_EVAT <- get_percent_variation(PCA_EVAT)

# __b) Q-mode

qPV_FLAT <- get_percent_variation(qPCA_FLAT)
qPV_BRAT <- get_percent_variation(qPCA_BRAT)
qPV_IWAT <- get_percent_variation(qPCA_IWAT)
qPV_LIVT <- get_percent_variation(qPCA_LIVT)
qPV_EVAT <- get_percent_variation(qPCA_EVAT)

# _3.) Plot Principal Component Analyses ----

# __a.) Plot FLAT PCs (R mode - enough variation)

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

# __b.) Plot BRAT PCs

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

# __c.) Plot IWAT PCs

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

# __d.) Plot EVAT PCs

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

# __e.) Plot LIVT PCs

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

# __f.) Plot for talks (in R mode)

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

# Inspect obesity-related genes among tissues ----

# Will be  used to warrant DEG search, PCA analysis will move to supplement. As
# per 10.3390/ijms231911005 using: Leptin (LEP), the leptin receptor (LEPR),
# proopiomelanocortin (POMC), prohormone convertase 1 (PCSK1), the melanocortin
# 4 receptor (MC4R), single-minded homolog 1 (SIM1), brain-derived neurotrophic
# factor (BDNF), and the neurotrophic tyrosine kinase receptor type 2 gene
# (NTRK2), NEED TO CONFIRM THAT THESE GENES ARE RELAVNT FOR FAT

# _1.) Compile data ----

# __a) SE objects

SE_all_tissues_all_genes <- list(BRAT, IWAT, LIVT, EVAT)
names(SE_all_tissues_all_genes) <- c("BRAT", "IWAT", "LIVT", "EVAT")

# __b) Obesity genes 

# genes from Mahmoud, Ranim, Virginia Kimonis, and Merlin G. Butler. 2022. “Genetics of Obesity in Humans: A Clinical Review.” International Journal of Molecular Sciences 23 (19): 11005. https://doi.org/10.3390/ijms231911005.
obesity_genes <- c("LEP", "LEPR", "POMC", "PCSK1", "MC4R", "SIM1", "BDNF", "NTRK2") # from 10.3390/ijms231911005

# genes from Hua, Yuchen, Danyingzhu Xie, Yugang Zhang, Ming Wang, Weiheng Wen, and Jia Sun. 2023. “Identification and Analysis of Key Genes in Adipose Tissue for Human Obesity Based on Bioinformatics.” Gene 888 (December):147755. https://doi.org/10.1016/j.gene.2023.147755.
obesity_genes <- c("EGR2", "GREM1", "NPY1R", obesity_genes)

# genes from Dahlman, Ingrid, and Peter Arner. 2010. “Chapter 3 - Genetics of Adipose Tissue Biology.” In Progress in Molecular Biology and Translational Science, edited by Claude Bouchard, 94:39–74. Genes and Obesity. Academic Press. https://doi.org/10.1016/B978-0-12-375003-7.00003-0.
obesity_genes <- c("ADRB2","GPR74", "GPR74", "PPARG", "SREBP1", "ADIPOQ","ARL15", obesity_genes)

warning(c("Confirm that these genes are relevant for fat: ", paste(obesity_genes,  sep=" ", collapse=" ")))

# _2.) Filter data to obesity genes and check ----

SE_all_tissues_obs_genes <- lapply(SE_all_tissues_all_genes, function(seob) seob[   which( rowData(seob)[["SYMBOL"]] %in% obesity_genes ) ] )
lapply(SE_all_tissues_obs_genes,  function(seob) toupper(rowData(seob)[["SYMBOL"]]))

# _3.) Show and plot obesity genes ----

# https://github.com/plger/sechm
# https://www.bioconductor.org/packages/release/bioc/vignettes/sechm/inst/doc/sechm.html

# get a compound plot
plot_heatmap_all_tissues_obesity_genes <- ggarrange(
  grid.grabExpr(draw(get_one_heatmap(SE_all_tissues_obs_genes[[1]]))),
  grid.grabExpr(draw(get_one_heatmap(SE_all_tissues_obs_genes[[2]]))),
  grid.grabExpr(draw(get_one_heatmap(SE_all_tissues_obs_genes[[3]]))),
  grid.grabExpr(draw(get_one_heatmap(SE_all_tissues_obs_genes[[4]]))),
  labels = list("a: BRAT", "b: IWAT", "c: LIVT", "d: EVAT"),
  font.label = list(size = 14, face = "bold"),
  ncol = 2, nrow = 2, hjust = "0", vjust = "0")

# show plot
plot_heatmap_all_tissues_obesity_genes

# save plot
ggsave(plot = plot_heatmap_all_tissues_obesity_genes, path = here("plots"), 
       filename = "055_r_array_analysis__plot_expr_obesity_flat.pdf",  
       width = 180, height = 85, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_heatmap_all_tissues_obesity_genes, path = here("../manuscript/display_items"), 
       filename = "055_r_array_analysis__plot_expr_obesity_flat_unassigned.pdf",  
       width = 180, height = 85, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# Inspect DEGs  ---- 

# Get Suppl. Tables 1-6.  Look for DEGs for 3 (and all) contrasts (CD CD against WD CD, CD WD, and WD WD) in each of the 4 tissues (results in 12 lists).

# _1.) Compile data ----

# __a) Use these full objects ----

SE_all_tissues_all_genes
names(SE_all_tissues_all_genes) 

# __b) Use these obesity genes if needed ----

obesity_genes 
warning(c("Confirm that these genes are relvant for fat: ", paste(obesity_genes,  sep=" ", collapse=" ")))

# __c) Use these objects which only contain obesity-relevant genes  ----

SE_all_tissues_obs_genes
names(SE_all_tissues_obs_genes) 

# _2.) Isolate DEGs  ----

# __a) For all genes ----

SE_all_tissues_all_genes <- lapply(SE_all_tissues_all_genes, function(se_ob) get_deg_lists(se_ob, peval = 0.05, logfc = 2))

# __b) For genes of interest related to obesity  ----

SE_all_tissues_obs_genes <- lapply(SE_all_tissues_obs_genes, function(se_ob) get_deg_lists(se_ob, peval = 0.05, logfc = 2))

# __c) Flatten DEG lists to tibbles ----

warning("List of tissues is flatetned here.")
DEGs_all_tissues_all_genes <- bind_rows(lapply(SE_all_tissues_all_genes, get_flat_deg_lists))
DEGs_all_tissues_obs_genes <- bind_rows(lapply(SE_all_tissues_obs_genes, get_flat_deg_lists))

# _3.) Check DEG lists ----

# __a) Look at lists ----

DEGs_all_tissues_all_genes %<>% arrange(TISSUE, CONTRASTS, abs(LOGFC)) 
DEGs_all_tissues_obs_genes %<>% arrange(TISSUE, CONTRASTS, abs(LOGFC)) 

DEGs_all_tissues_all_genes %>% print(n = Inf) 
DEGs_all_tissues_obs_genes %>% print(n = Inf)

# for nkb look at numerical values

# __b) Establish whether there is an intersection between full DEG lists and DEG list derived from obesity - only gene lists ----

if (identical(obesity_genes[which(obesity_genes %in% DEGs_all_tissues_all_genes[["SYMBOL"]])], character(0))){
  message("No obesity relavent genes found among full DEG list, consider providing more target genes.")
  } else {
  message("Found ", paste(obesity_genes[which(obesity_genes %in% DEGs_all_tissues_all_genes[["SYMBOL"]])], sep=" ", collapse=" " ) )
  
}

if (identical(obesity_genes[which(obesity_genes %in% DEGs_all_tissues_obs_genes[["SYMBOL"]])], character(0))){
  message("No obesity relavent genes found among full DEG list, consider providing more target genes.")
} else {
  message("Found ", paste(obesity_genes[which(obesity_genes %in% DEGs_all_tissues_obs_genes[["SYMBOL"]])], sep=" ", collapse=" " ), " in obesity specific lists." )
  
}

# __c) Export results as Excel table ----

write_xlsx(DEGs_all_tissues_all_genes, path = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/055_r_array_analysis__degs_all_tissues_all_genes.xlsx", col_names = TRUE, format_headers = TRUE)
write_xlsx(DEGs_all_tissues_obs_genes, path = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/055_r_array_analysis__degs_all_tissues_obs_genes.xlsx", col_names = TRUE, format_headers = TRUE)

get_raw_summaries <- function(se) {
  
  require("SummarizedExperiment")
  
  # se <- SE_all_tissues_obs_genes[[1]]
  
  rownames(assays(se, withDimnames=FALSE)[["exprs"]])  <- rowData(se)[which(rownames(assay(se)) %in% rownames(rowData(se))) ,"SYMBOL"]
  colnames(assays(se, withDimnames=FALSE)[["exprs"]])  <- colData(se)[which(colnames(assay(se)) %in% rownames(colData(se))) ,"ParentalDietMoFa"]
  
  return(assays(se, withDimnames=FALSE)[["exprs"]])
}

# extra export for NKB only
raw_expression__all_tissues_obs_genes <- lapply(SE_all_tissues_obs_genes, get_raw_summaries)
capture.output(raw_expression__all_tissues_obs_genes, file = "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/055_r_array_analysis__raw_summary__all_tissues_obs_genes.txt")

# __d) Get Upset Plot ----

# >>> Code below needs outlining again ----
stop("No obesity genes found differentially expressed - reoutline code below.")



# [continue here]

# Get Fig. 3: Upset plot of intersections of full list for all contrasts

# Subset DEG  lists to set differences among all contrasts in each of the 4 tissues. ----

# Get Suppl. Tables 6-12. Look for DEGs for 3 (and all) contrasts (CD CD against WD CD, CD WD, and WD WD) in each of the 4 tissues (results in 12 lists).

# Describe obesity related genes in intersections of set different DEG lists ----

#  Get Fig. 4: Upset plot of intersections of full list for all contrasts

# GSEA of set differences with unique obesity genes

# Get Suppl. Tab 12-n 

# Save Environment ----

stop("Unadjusted old analysis code below - possibly integrate into code above. ")

# >>> Unadjusted old analysis code below - possibly integrate into code above. ----

# _2.) DGE analysis using Limma ----

# see DESeq2 tutorial at https://colauttilab.github.io/RNA-Seq_Tutorial.html
# see GCSCore tutaorial at https://www.bioconductor.org/packages/release/bioc/vignettes/GCSscore/inst/doc/GCSscore.pdf
# see Limma slides at https://s3.amazonaws.com/assets.datacamp.com/production/course_6456/slides/chapter1.pdf
# **use this!** - Limma slides at  https://kasperdanielhansen.github.io/genbioconductor/html/limma.html
# see also here for LFC shrinkage and genral procedure https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/08_practical_DE.pdf
# see here why shrinking LFCs is considered unenecssary by the limma authors - https://support.bioconductor.org/p/100804/



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

