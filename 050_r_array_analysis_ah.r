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

library(ggplot2)
library(stringr)
library(tidyr)
library(Biobase)
library(BiocGenerics)
library(gplots)
library(dplyr)
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

# _3.) Functions ----

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
  pData(expression_set) <- left_join((pData(expression_set) %>% select(
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

# _4.) Color code for plotting ----

hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(n=11, name="RdBu"))(50))

# Load and shape data ----

# _1.) Loading normalized data as ExpressionSet type ----

# Data are Clariom S mouse arrays
# I removed strong outliers: A285_liver clustert zu bAT, A339_liver weit weg vom Rest

base::load("/Users/paul/Documents/HM_MouseMating/analysis_ah/allTissues_normData.RData") # only if you are interrested to look into the normalized data

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
saveRDS(LIAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_ingwat.rds"))
saveRDS(EVAT, file = here("rds_storage", "050_r_array_analysis__normalized_data_liv.rds"))

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






# AH code below - continue here after 30.05.2023

# Volcano plot ############################# 
pVal <- resultTable$adj.P.Val   #adj. p-Value für cutoffs verwenden
CO1 <- 0.05                     #adj. p-Value cutoff
CO2 <- 0.5                      #FC cutoff
#vermutlich kannst du stringenter sein und als adj. P 0.01 nehmen, da wir sehr viele Hits haben
#(vielleicht mehr eingrenzen, um die Analysen zu vereinfachen)
#FC cutoff kannst du evt. auf 0.75 hochnehmen (bei 1 bleibt nicht mehr so viel übrig), aber p-Wert ist besser

load("Results_05.23/DGELists.RData") #lädt alle DGE Tabellen
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


