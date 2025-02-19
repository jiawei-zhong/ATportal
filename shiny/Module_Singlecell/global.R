# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(ggplotify))
suppressMessages(library(cowplot))
suppressMessages(library(scCustomize))
suppressMessages(library(ggplot2))
suppressMessages(library(ggplot2))

#read data

#colors

discrete_color <- c("T, NK & NKT" = "#8ca5a5", "B" = "#7dbfb3", "mono. & macro." = "#4f736a", "mast" = "#1a4659", "FAPs" = "#f59b7c", "adipocytes" = "#ba5c28", "vascular" = "#e2c744")

discrete_color_gradient <- colorRampPalette(discrete_color)

#heatmap green to purple
PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#e2c744", "#1a4659", "#ba5c28", "#4f736a", "#f0cbbf", "#edf2f2" ,"#f2ecde", "#adbec4", "#c3b6b6" , "#8ca5a5", "#b2dfda", "#f59b7c", "#7dbfb3" ,"#939a64")
portal_col_ex <-colorRampPalette(portalcol2)
farben <- list("default"= portalcol,"PRGn" = PRGn, "RdBu" = RdBu, "viridis"=virid)

mode_in <- setNames(c("raw","normalized"),c("mean abundance", "normalized mean abundance"))



Massier_et_al_all <- readRDS("./data/META_all.sub.RDS")
Massier_et_al_Lymphoid <- readRDS("./data/WAT_all_Lymphoid.RDS")
Massier_et_al_Myeloid <- readRDS("./data/WAT_all_Myeloid.RDS")
Massier_et_al_Vascular <- readRDS("./data/WAT_all_Vascular.RDS")
Massier_et_al_B <- readRDS("./data/WAT_all_B.RDS")
Massier_et_al_FAPs_sc <- readRDS("./data/WAT_sc_FAPs.RDS")
Massier_et_al_FAPs_om <- readRDS("./data/WAT_om_FAPs.RDS")
Massier_et_al_FAPs_pvat <- readRDS("./data/WAT_pvat_FAPs.RDS")



temp <- as.character(Idents(Massier_et_al_all))
temp[temp=="Lymphoid"] <- "T, NK & NKT"
temp[temp=="B_cell"] <- "B"
temp[temp=="Myeloid"] <- "mono. & macro."
temp[temp=="Mast_cell"] <- "mast"
temp[temp=="Adipocyte"] <- "adipocytes"
temp[temp=="Vascular"] <- "vascular"
temp <- factor(temp,levels = c("T, NK & NKT", "B", "mono. & macro.", "mast", "FAPs", "Mesothelial", "adipocytes", "vascular"))
Idents(Massier_et_al_all) <- temp
rm(temp)


Massier_et_al_all_marker <- readRDS("./data/META_all_marker.RDS")
Massier_et_al_Lymphoid_marker <- readRDS("./data/WAT_all_Lymphoid_marker.RDS")
Massier_et_al_Myeloid_marker <- readRDS("./data/WAT_all_Myeloid_marker.RDS")
Massier_et_al_Vascular_marker <- readRDS("./data/WAT_all_Vascular_marker.RDS")
Massier_et_al_B_marker <- readRDS("./data/WAT_all_B_marker.RDS")
Massier_et_al_FAPs_sc_marker <- readRDS("./data/WAT_sc_FAPs_marker.RDS")
Massier_et_al_FAPs_om_marker <- readRDS("./data/WAT_om_FAPs_marker.RDS")
Massier_et_al_FAPs_pvat_marker <- readRDS("./data/WAT_pvat_FAPs_marker.RDS")


Hinte_et_al_scAT_LTSS <- readRDS("./data/scAT_LTSS.rds")
Hinte_et_al_scAT_NEFA <- readRDS("./data/scAT_NEFA.rds")
Hinte_et_al_omAT_LTSS <- readRDS("./data/omAT_LTSS.rds")
Hinte_et_al_omAT_MTSS <- readRDS("./data/omAT_MTSS.rds")



Reinisch_et_al_sc <- readRDS("./data/scAT_MHUO.rds")
Reinisch_et_al_sc_Adipo <- readRDS("./data/scAdipo_MHUO.rds")
Reinisch_et_al_sc_APCs <- readRDS("./data/scAPCs_MHUO.rds")
Reinisch_et_al_sc_Immune <- readRDS("./data/scImmuneCs_MHUO.rds")



Reinisch_et_al_vis <- readRDS("./data/visAT_MHUO.rds")
Reinisch_et_al_vis_Adipo <- readRDS("./data/visAdipo_MHUO.rds")
Reinisch_et_al_vis_APCs <- readRDS("./data/visAPCs_MHUO.rds")
Reinisch_et_al_vis_Immune <- readRDS("./data/visImmuneCs_MHUO.rds")
Reinisch_et_al_vis_Meso <- readRDS("./data/visMesoCs_MHUO.rds")



# Reinisch_et_al_sc <- readRDS("./data/scAT_MHUO.rds")
# Reinisch_et_al_sc_Adipo <- readRDS("./data/scAdipo_MHUO.rds")
# Reinisch_et_al_sc_APCs <- readRDS("./data/scAPCs_MHUO.rds")
# Reinisch_et_al_sc_Immune <- readRDS("./data/scImmuneCs_MHUO.rds")
# Reinisch_et_al_vis <- readRDS("./data/visAT_MHUO.rds")
# Reinisch_et_al_vis_Adipo <- readRDS("./data/visAdipo_MHUO.rds")
# Reinisch_et_al_vis_APCs <- readRDS("./data/visAPCs_MHUO.rds")
# Reinisch_et_al_vis_Immune <- readRDS("./data/visImmuneCs_MHUO.rds")
# Reinisch_et_al_vis_Meso <- readRDS("./data/visMesoCs_MHUO.rds")
# 
# 
# Reinisch_et_al_sc[['SCT']] <- NULL
# Reinisch_et_al_sc_Adipo[['SCT']] <- NULL
# Reinisch_et_al_sc_APCs[['SCT']] <- NULL
# Reinisch_et_al_sc_Immune[['SCT']] <- NULL
# Reinisch_et_al_vis[['SCT']] <- NULL
# Reinisch_et_al_vis_Adipo[['SCT']] <- NULL
# Reinisch_et_al_vis_APCs[['SCT']] <- NULL
# Reinisch_et_al_vis_Immune[['SCT']] <- NULL
# Reinisch_et_al_vis_Meso[['SCT']] <- NULL
# 
# 
# Reinisch_et_al_sc[['integrated']] <- NULL
# Reinisch_et_al_sc_Adipo[['integrated']] <- NULL
# Reinisch_et_al_sc_APCs[['integrated']] <- NULL
# Reinisch_et_al_sc_Immune[['integrated']] <- NULL
# Reinisch_et_al_vis[['integrated']] <- NULL
# Reinisch_et_al_vis_Adipo[['integrated']] <- NULL
# Reinisch_et_al_vis_APCs[['integrated']] <- NULL
# Reinisch_et_al_vis_Immune[['integrated']] <- NULL
# Reinisch_et_al_vis_Meso[['integrated']] <- NULL
# 
# 
# saveRDS(Reinisch_et_al_sc,"./data/scAT_MHUO.rds")
# saveRDS(Reinisch_et_al_sc_Adipo,"./data/scAdipo_MHUO.rds")
# saveRDS(Reinisch_et_al_sc_APCs,"./data/scAPCs_MHUO.rds")
# saveRDS(Reinisch_et_al_sc_Immune,"./data/scImmuneCs_MHUO.rds")
# saveRDS(Reinisch_et_al_vis,"./data/visAT_MHUO.rds")
# saveRDS(Reinisch_et_al_vis_Adipo,"./data/visAdipo_MHUO.rds")
# saveRDS(Reinisch_et_al_vis_APCs,"./data/visAPCs_MHUO.rds")
# saveRDS(Reinisch_et_al_vis_Immune,"./data/visImmuneCs_MHUO.rds")
# saveRDS(Reinisch_et_al_vis_Meso,"./data/visMesoCs_MHUO.rds")





gene_all <- sort(rownames(Massier_et_al_all))
gene_all <- c(gene_all[grep("^AC\\d{6}|^AD\\d{6}|^AE\\d{6}|^AF\\d{6}|^AJ\\d{6}|^AL\\d{6}|^AP\\d{6}|^CATG\\d{6}",gene_all,invert = T)],
              gene_all[grep("^AC\\d{6}|^AD\\d{6}|^AE\\d{6}|^AF\\d{6}|^AJ\\d{6}|^AL\\d{6}|^AP\\d{6}|^CATG\\d{6}",gene_all,invert = F)])


citation <- readRDS("./data/citation.RDS")


# omAT_LTSS <- readRDS("omAT_LTSS.rds")
# omAT_MTSS <- readRDS("omAT_MTSS.rds")
# scAT_LTSS <- readRDS("scAT_LTSS.rds")
# scAT_NEFA <- readRDS("scAT_NEFA.rds")

# omAT_LTSS[['SCT']] <- NULL
# omAT_MTSS[['SCT']] <- NULL
# scAT_LTSS[['SCT']] <- NULL
# scAT_NEFA[['SCT']] <- NULL

# omAT_LTSS[['integrated']] <- NULL
# omAT_MTSS[['integrated']] <- NULL
# scAT_LTSS[['integrated']] <- NULL
# scAT_NEFA[['integrated']] <- NULL


# saveRDS(omAT_LTSS,"omAT_LTSS.rds")
# saveRDS(omAT_MTSS,"omAT_MTSS.rds")
# saveRDS(scAT_LTSS,"scAT_LTSS.rds")
# saveRDS(scAT_NEFA,"scAT_NEFA.rds")


### save legend (以防输出空白的画板，导致在vscode运行的时候报错)


Cohort_bl__Subcluster_bl <- list(Massier_et_al=c("all" = "all","T, NK & NKT" = "Lymphoid", "mono. & macro." = "Myeloid", "vascular" = "Vascular", "B" = "B", "FAPs subcutaneous"="FAPs_sc", "FAPs omental"="FAPs_om", "FAPs perivascular"="FAPs_pvat"),
                                 Hinte_et_al=c("omAT_LTSS" = "omAT_LTSS", "omAT_MTSS" = "omAT_MTSS", "scAT_LTSS" = "scAT_LTSS", "scAT_NEFA" = "scAT_NEFA"),
                                 Reinisch_et_al=c("sc" = "sc", "sc_Adipo" = "sc_Adipo", "sc_APCs" = "sc_APCs", "sc_Immune" = "sc_Immune", "vis" = "vis", "vis_Adipo" = "vis_Adipo", "vis_APCs" = "vis_APCs", "vis_Immune" = "vis_Immune", "vis_Meso" = "vis_Meso"))

# for (cohort in names(Cohort_bl__Subcluster_bl)) {
#   for (subcluster in Cohort_bl__Subcluster_bl[[cohort]]) {
#     legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(paste0(cohort,"_",subcluster)), repel = T)[[1]] +
#                                      theme(aspect.ratio=1,
#                                            text = element_text(family = "Red Hat Display"),
#                                            legend.position = "bottom",
#                                            legend.text = element_text(size = 10, family = "Red Hat Display")) +
#                                      scale_color_manual(values = discrete_color_gradient(length(levels(get(paste0(cohort,"_",subcluster)))))) +
#                                      NoAxes() +
#                                      guides(color = guide_legend(nrow = 3, override.aes = list(size=4)))))
#     assign(x = paste0("legend_",cohort,"_",subcluster), value = legend)
#     saveRDS(legend, paste0("./data/legend_",cohort,"_",subcluster,'.RDS'))
#   }
# }


for (cohort in names(Cohort_bl__Subcluster_bl)) {
  for (subcluster in Cohort_bl__Subcluster_bl[[cohort]]) {
    assign(x = paste0("legend_",cohort,"_",subcluster),value = readRDS(paste0("./data/legend_",cohort,"_",subcluster,'.RDS')))
  }
}


