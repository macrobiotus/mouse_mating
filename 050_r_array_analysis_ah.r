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

# rewrite metadata in exprssion data sets
adjust_array_data = function(expression_set, model_variables){
  
  require(dplyr)
  require(Biobase)
  
  # adjust column and row names names in expression data 
  colnames(expression_set) <- colnames(expression_set) %>% 
    str_replace_all("bAT",    "BRAT") %>%
    str_replace_all("ingWAT", "SCAT") %>% 
    str_replace_all("liver",  "LIAT") %>% 
    str_replace_all("eWAT",   "EVAT")
  
  # adjust "Tissue" column values  in expression data 
  pData(expression_set) %<>% mutate(Tissue = case_when(
    Tissue == "bAT"    ~ "BRAT",
    Tissue == "ingWAT" ~ "SCAT",
    Tissue == "liver"  ~ "LIAT",
    Tissue == "eWAT"   ~ "EVAT"))
  
  # merge obesity variables from modelling  to metadata from array experiments
  pData(expression_set) <- left_join( 
    (pData(expression_set) %>% select(-c(Sex, Parental_diet, Diet_mother, Diet_father, Group))),
    (model_variables), by = c( "Animal" = "AnimalId"))
  
  return(expression_set)
  
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


# Principal Component Analysis ----

# _1.) Get PC for all tissues

PCA_FLAT <- prcomp( t (exprs(FLAT)), scale. = FALSE)
percentVar <- round(100*PCA_FLAT$sdev^2 / sum(PCA_norm$sdev^2), 1)

# _2.) Plot PCs in 2D for sevral variables types  ----

# "Overall expression differences between analysed tissues among f1 offspring"
plot_pca_flat_a  <- fviz_pca_ind(PCA_FLAT, label="none", habillage = factor(pData(FLAT)$Tissue),
                   pointsize = 2,
                   palette = c("firebrick3","purple","steelblue3","gold3"),
                   legend.title = "Tissue",
                   invisible = "none",
                   addEllipses = FALSE, 
                   title = "a") +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance"))+
  theme_bw() +
  scale_shape_manual(values=c(19,19,19,19))

# Overall expression differences and obesity status among f1 offspring
plot_pca_flat_b <- fviz_pca_ind(PCA_FLAT, label="none", habillage = factor(pData(FLAT)$ObesityLgcl),
                  pointsize = 2,
                  palette = c("firebrick3","purple","steelblue3","gold3"),
                  legend.title = "F1 Obesity",
                  invisible = "none",
                  addEllipses = FALSE, 
                  title = "b") +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance"))+
  theme_bw() +
  scale_shape_manual(values=c(19,19,19,19))

# Overall expression differences and obesity parental obesity status among f0
plot_pca_flat_c <- fviz_pca_ind(PCA_FLAT, label="none", habillage = factor(pData(FLAT)$ObeseParents),
                  pointsize = 2,
                  palette = c("firebrick3","purple","steelblue3","gold3"),
                  legend.title = "F0 Obesity",
                  invisible = "none",
                  addEllipses = FALSE, 
                  title = "c") +
  labs(x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance"))+
  theme_bw() +
  scale_shape_manual(values=c(19,19,19,19))

plot_pca_flat <- ggarrange(plot_pca_flat_a, plot_pca_flat_b, plot_pca_flat_c, nrow = 1, ncol = 3)

ggsave(plot = plot_pca_flat, path = here("plots"), 
       filename = "050_r_array_analysis__plot_pca_flat.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)
ggsave(plot = plot_pca_flat, path = here("../manuscript/display_items"), 
       filename = "050_r_array_analysis__plot_pca_flat_unassigned.pdf",  
       width = 180, height = 65, units = "mm", dpi = 300,  limitsize = TRUE, scale = 2)

# continue here after 26.05.2023 - commit first -- older code below ----


# _3.) PCA of expression values ----

PCA_norm <- prcomp(t(exprs(normData)),scale. = TRUE)
dataGG <- data.frame(PC1 = PCA_norm$x[,1], PC2 = PCA_norm$x[,2],PC3 = PCA_norm$x[,3])

# _4.) Plot PCs in 3D for tissue types  ----

p <- plot_ly(dataGG, x = ~PC1, y = ~PC2, z = ~PC3) %>%
  add_markers(color = ~pData(normData)$Tissue,colors = c("firebrick3","purple","steelblue3","gold3"),
              symbol = ~pData(normData)$Tissue, symbols = c(19,19,19,19))

saveWidget(p, file = paste0(here("plots"),"/", "050_r_array_analysis__PCA-3D.html"), selfcontained = T, libdir = "lib")

# Old code below - continue here after 24.05.2023 ---- 
# - env should be up-to-date and re-run is painless if needed

# PCAs for single tissue ############################# 

load("Results_05.23/normData4DGE.RData") #für jedes Gewebe die normalisierten Daten

#Beispiel Analyse für bAT (gleiche für die anderen Gewebe)
normData <- bAT_normData #choose one of the 4 tissues: bAT_normData, Liv_normData,ingWAT_normData, eWAT_normData
para <- "bAT" #choose: bAT, Liver, ingWAT, eWAT. (für Path beim Speichern der Plots)

group <- pData(normData)$Parental_diet
group <- gsub("_", ", ", group)

#2D PCA
PCA_norm <- prcomp(t(exprs(normData)),scale. = FALSE)
percentVar <- round(100*PCA_norm$sdev^2/sum(PCA_norm$sdev^2),1)

fviz_pca_ind(PCA_norm, label="none", habillage=factor(group),
             pointsize = 2,
             title = "",
             palette = c("firebrick3","purple","steelblue3","gold3"),
             legend.title = "Parental diet (m/f)",
             invisible="quali")+
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), 
       y = paste0("PC2: ",percentVar[2],"% variance"))+
  theme_test()+
  scale_shape_manual(values=c(19,19,19,19))
ggsave(paste0("Results_05.23/",para,"_PCA.pdf"),width = 5,height = 4)

#Variance PCA plots
phenos <- pData(normData)[,c(4:6)]
phenos$Parental_diet <- as.factor(phenos$Parental_diet)
phenos$Diet_mother <- as.factor(phenos$Diet_mother)
phenos$Diet_father <- as.factor(phenos$Diet_father)
colnames(phenos) <- c("Parental diet","Diet mother", "Diet father")
res1 <- prince(exprs(normData),phenos,top=5)

pdf(paste0("Results_05.23/",para,"_PCA_Variance.pdf"),width = 5,height = 5)
  prince.plot(prince=res1,note=TRUE,key = FALSE)
dev.off()

## Volcano plot ############################# 
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


