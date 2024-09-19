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


# Proteomicseset <- readRDS("./data/Proteomicseset.rds")

sex_df_RNA <- readRDS("./data/sex_df_RNA.RDS")

sex_df_RNA <- sex_df_RNA[!is.na(sex_df_RNA$SMD),]

sex_df_protein <- readRDS("./data/sex_df_protein.RDS")

sex_df_protein <- sex_df_protein[!is.na(sex_df_protein$SMD),]

portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")