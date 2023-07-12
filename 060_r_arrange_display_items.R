# Get Display Items for Mouse Mating manuscript  

# see https://stackoverflow.com/questions/57020571/reading-pdf-plots-arranging-them-on-a-grid-save-in-one-page-pdf-using-r
# see https://cran.r-project.org/web/packages/magick/vignettes/intro.html

# Prepare environment  ---- 

# _1.) Collect garbage ----

rm(list=ls())

gc()

# _2.) Packages ----

library(magick) # see https://cran.r-project.org/web/packages/magick/vignettes/intro.html
library(dplyr)

# Read in data ----

# _1.) Fig. 2 ----

# created by 050_r_array_analysis

path_volcano_brat <- c("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__plot_volcano_brat.pdf")
path_volcano_liat <- c("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__plot_volcano_liat.pdf")

path_heatmap_brat <- c("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__plot_heatmap_brat.pdf")
path_heatmap_liat <- c("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__plot_heatmap_liat.pdf")

# Create text figures ----

# _1.) Fig. 2 ----

panel.ul <- image_trim(image_read_pdf(path_volcano_brat)) %>% image_annotate(" a", size = 16, gravity = "northwest", color = "black")
panel.ur <- image_trim(image_read_pdf(path_volcano_liat)) %>% image_annotate(" b", size = 16, gravity = "northwest", color = "black")
panel.ll <- image_trim(image_read_pdf(path_heatmap_brat)) %>% image_annotate(" c", size = 16, gravity = "northwest", color = "black")
panel.lr <- image_trim(image_read_pdf(path_heatmap_liat)) %>% image_annotate(" d", size = 16, gravity = "northwest", color = "black")

fig_1 <- image_append(
  c(
    image_append(c(panel.ul, panel.ll), stack = TRUE), 
    image_append(c(panel.ur, panel.lr), stack = TRUE)
    ))

image_write(fig_1, format = "png", density = 200, path =  "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/060_r_arrange_display_items__fig_2.png")

# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "060_r_arrange_display_items.RData"))
renv::snapshot()
