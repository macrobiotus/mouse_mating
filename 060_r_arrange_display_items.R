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

# _1.) Fig. 3 ----

# created by 050_r_array_analysis

path_volcano_brat <- c("/Users/paul/Documents/HM_MouseMating/manuscript/display_items/050_r_array_analysis__plot_volcano_brat.pdf")
path_volcano_liat <- c("/Users/paul/Documents/HM_MouseMating/manuscript/display_items//050_r_array_analysis__plot_volcano_liat.pdf")

path_heatmap_brat <- c("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__plot_heatmap_brat.pdf")
path_heatmap_liat <- c("/Users/paul/Documents/HM_MouseMating/analysis/plots/050_r_array_analysis__plot_heatmap_liat.pdf")

# _2.) Fig. 4 ----

path_kegg_liat <- c("/Users/paul/Documents/HM_MouseMating/manuscript/display_items/050_r_array_analysis__plot_kegg_L_MotherFatherObese_vs_MotherFatherNotObese.pdf")
path_kegg_brat <- c("/Users/paul/Documents/HM_MouseMating/manuscript/display_items/050_r_array_analysis__plot_kegg_BAT_MotherFatherObese_vs_FatherObese.pdf")

# _3.) Fig. 5 ----

path_go_liat <- c("/Users/paul/Documents/HM_MouseMating/manuscript/display_items/050_r_array_analysis__plot_go_L_MotherFatherObese_vs_MotherFatherNotObese.pdf")
path_go_brat <- c("/Users/paul/Documents/HM_MouseMating/manuscript/display_items/050_r_array_analysis__plot_go_BAT_MotherFatherObese_vs_FatherObese.pdf")

# Create text figures ----

# _1.) Fig. 3 ----

panel.ul <- image_trim(image_read_pdf(path_volcano_brat)) %>% image_annotate(" a", size = 16, gravity = "northwest", color = "black")
panel.ur <- image_trim(image_read_pdf(path_volcano_liat)) %>% image_annotate(" b", size = 16, gravity = "northwest", color = "black")
panel.ll <- image_trim(image_read_pdf(path_heatmap_brat)) %>% image_annotate(" c", size = 16, gravity = "northwest", color = "black")
panel.lr <- image_trim(image_read_pdf(path_heatmap_liat)) %>% image_annotate(" d", size = 16, gravity = "northwest", color = "black")

fig_1 <- image_append(
  c(
    image_append(c(panel.ul, panel.ll), stack = TRUE), 
    image_append(c(panel.ur, panel.lr), stack = TRUE)
    ))

image_write(fig_1, format = "png", density = 300, path =  "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/060_r_arrange_display_items__fig_2.png")


# _2.) Fig. 4 ----

fig_4 <- image_append(
  c(
    image_trim(image_read_pdf(path_kegg_brat)) %>% image_annotate(" a", size = 16, gravity = "northwest", color = "black"),
    image_trim(image_read_pdf(path_kegg_liat)) %>% image_annotate(" b", size = 16, gravity = "northwest", color = "black")
  ))

image_write(fig_4, format = "png", density = 300, path =  "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/060_r_arrange_display_items__fig_4.png")



# _3.) Fig. 5 ----

fig_5 <- image_append(
  c(
    image_trim(image_read_pdf(path_go_brat)) %>% image_annotate(" a", size = 16, gravity = "northwest", color = "black"),
    image_trim(image_read_pdf(path_go_liat)) %>% image_annotate(" b", size = 16, gravity = "northwest", color = "black")
  ))

image_write(fig_5, format = "png", density = 300, path =  "/Users/paul/Documents/HM_MouseMating/manuscript/display_items/060_r_arrange_display_items__fig_5.png")


# Snapshot environment ----

sessionInfo()
save.image(file = here("scripts", "060_r_arrange_display_items.RData"))
renv::snapshot()
