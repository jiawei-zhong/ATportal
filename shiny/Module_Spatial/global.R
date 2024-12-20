# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
#read data
# ## Spatial ----

Jesper_et_al_baseline <- readRDS("./data/Jesper_et_al_baseline.RDS")
Jesper_et_al_insulin <- readRDS("./data/Jesper_et_al_insulin.RDS")

Jesper_et_al_baseline_marker <- readRDS("./data/Jesper_et_al_baseline_marker.RDS")
Jesper_et_al_baseline_marker <- Jesper_et_al_baseline_marker[Jesper_et_al_baseline_marker$avg_log2FC>0.25 & (Jesper_et_al_baseline_marker$pct.1>0.1 | Jesper_et_al_baseline_marker$pct.2>0.1) & Jesper_et_al_baseline_marker$p_val_adj<0.05,]
Jesper_et_al_insulin_marker <- readRDS("./data/Jesper_et_al_insulin_marker.RDS")
Jesper_et_al_insulin_marker <- Jesper_et_al_insulin_marker[Jesper_et_al_insulin_marker$avg_log2FC>0.25 & (Jesper_et_al_insulin_marker$pct.1>0.1 | Jesper_et_al_insulin_marker$pct.2>0.1) & Jesper_et_al_insulin_marker$p_val_adj<0.05,]




for (i in grep("marker$",ls(),value = T)) {
  eval(parse(text=paste0(i,'$p_val <- signif(',i,'$p_val,4)')))
  eval(parse(text=paste0(i,'$p_val_adj <- signif(',i,'$p_val_adj,4)')))
  eval(parse(text=paste0(i,'$avg_log2FC <- signif(',i,'$avg_log2FC,4)')))
}

# Jesper_et_al_baseline_marker$p_val <- as.character(Jesper_et_al_baseline_marker$p_val)


gene_all <- sort(rownames(Jesper_et_al_baseline))

gene_all <- c(gene_all[grep("^AC\\d{6}|^AD\\d{6}|^AE\\d{6}|^AF\\d{6}|^AJ\\d{6}|^AL\\d{6}|^AP\\d{6}|^CATG\\d{6}",gene_all,invert = T)],
              gene_all[grep("^AC\\d{6}|^AD\\d{6}|^AE\\d{6}|^AF\\d{6}|^AJ\\d{6}|^AL\\d{6}|^AP\\d{6}|^CATG\\d{6}",gene_all,invert = F)])

names(gene_all) <- gene_all



cell_type_all <- colnames(Jesper_et_al_baseline@meta.data)[24:85]

names(cell_type_all) <- cell_type_all
#colors
#heatmap green to purple
PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#e2c744", "#1a4659", "#ba5c28", "#4f736a", "#f0cbbf", "#edf2f2" ,"#f2ecde", "#adbec4", "#c3b6b6" , "#8ca5a5", "#b2dfda", "#f59b7c", "#7dbfb3" ,"#939a64")
portal_col_ex <-colorRampPalette(portalcol2)
farben <- list("default"= portalcol,"PRGn" = PRGn, "RdBu" = RdBu, "viridis"=virid)

mode_in <- setNames(c("raw","normalized"),c("mean abundance", "normalized mean abundance")) 


discrete_color <- c("T, NK & NKT" = "#8ca5a5", "B" = "#7dbfb3", "mono. & macro." = "#4f736a", "mast" = "#1a4659", "FAPs" = "#f59b7c", "adipocytes" = "#ba5c28", "vascular" = "#e2c744")

discrete_color_gradient <- colorRampPalette(discrete_color)