# Load packages ----

suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(scCustomize))
suppressMessages(library(ggplot2))
suppressMessages(library(extrafont))


## Clinical ----

discrete_color <- c("#1a4659", "#e2c744", "#4f736a", "#edf2f2", "#b2dfda", "#7dbfb3", "#adbec4", "#ba5c28", "#f0cbbf", "#f2ecde", "#f59b7c", "#939a64", "#8ca5a5", "#c3b6b6")

discrete_color <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#edf2f2", "#f2ecde", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")

discrete_color <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")

discrete_color <- c("#8ca5a5", "#7dbfb3", "#4f736a", "#1a4659", "#f59b7c", "#ba5c28", "#e2c744")


discrete_color <- c("T, NK & NKT" = "#8ca5a5", "B" = "#7dbfb3", "mono. & macro." = "#4f736a", "mast" = "#1a4659", "FAPs" = "#f59b7c", "adipocytes" = "#ba5c28", "vascular" = "#e2c744")

discrete_color_gradient <- colorRampPalette(discrete_color)


META_all <- readRDS("./data/META_all.sub.RDS")

temp <- as.character(Idents(META_all))

temp[temp=="Lymphoid"] <- "T, NK & NKT"
temp[temp=="B_cell"] <- "B"
temp[temp=="Myeloid"] <- "mono. & macro."
temp[temp=="Mast_cell"] <- "mast"
temp[temp=="Adipocyte"] <- "adipocytes"
temp[temp=="Vascular"] <- "vascular"

temp <- factor(temp,levels = c("T, NK & NKT", "B", "mono. & macro.", "mast", "FAPs", "Mesothelial", "adipocytes", "vascular"))

Idents(META_all) <- temp

rm(temp)

p1 <- DimPlot_scCustom(seurat_object = META_all,label = T,raster=FALSE)[[1]] +
  theme(aspect.ratio=1,
        text = element_text(family = "Red Hat Display"),
        legend.text = element_text(size = 12, family = "Red Hat Display"),
        legend.position = "bottom") +
  # xlim(-10.79425,14.05299) +
  # ylim(-10.60139,12.84636) +
  # scale_color_manual(values = discrete_color) +
  scale_color_manual(values = discrete_color_gradient(length(levels(META_all)))) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size=4)))


