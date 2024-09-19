# Load packages ----

suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))



## Clinical ----



# meta_r <- readRDS("./data/meta_r.RDS")
# meta_p <- readRDS("./data/meta_p.RDS")

obese_df_RNA <- readRDS("./data/obese_df_RNA.RDS")

obese_df_RNA <- obese_df_RNA[!is.na(obese_df_RNA$SMD),]

obese_df_protein <- readRDS("./data/obese_df_protein.RDS")

obese_df_protein <- obese_df_protein[!is.na(obese_df_protein$SMD),]

portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")