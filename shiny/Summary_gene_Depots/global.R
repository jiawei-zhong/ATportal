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

depots_df_RNA <- readRDS("./data/depots_df_RNA.RDS")

depots_df_RNA <- depots_df_RNA[!is.na(depots_df_RNA$SMD),]

depots_df_protein <- readRDS("./data/depots_df_protein.RDS")

depots_df_protein <- depots_df_protein[!is.na(depots_df_protein$SMD),]

portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")
