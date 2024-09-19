library(shiny)
library(viridis)
library(RColorBrewer)
library(reshape)
library(dplyr)
# library(tidyverse)

data_line <- read.delim("./data/TimeCourse.txt")
# data_niga_heat <- read.delim("./data/TimeCourse_wide.txt", row.names = 1)
data_SVF <- read.delim("./data/TimeCourseSVF.txt")
data_SVF$timepoint <- factor(data_SVF$timepoint, levels = c("day4"   ,    "day8"   ,    "day12"  ,    "adipocytes"))
# data_SVF_heat <- read.delim("./data/TimeCourse_wideSVF.txt", row.names = 1)
# data_novo <- read.delim("./data/novoarray.txt")
# data_novo_heat <- read.delim("./data/novo_heat.txt", row.names = 1)
data_fantom <- read.delim("./data/fantom.txt", row.names = 1)
data_prot_heat <- read.delim("./data/NIGA_heat_proteome.txt", row.names = 1)
data_prot_raw <- read.delim("./data/NIGA_line_proteome_raw.txt")
data_prot_z <- read.delim("./data/NIGA_line_proteome_norm.txt")
data_SGBS_p_heat <- read.delim("./data/SGBS_heat_proteome.txt", row.names = 1)
data_SGBS_p_raw <- read.delim("./data/SGBS_line_proteome_raw.txt")
data_SGBS_p_z <- read.delim("./data/SGBS_line_proteome_norm.txt")
data_SGBS_t_heat <- read.delim("./data/SGBS_heat.txt", row.names = 1)
data_SGBS_t_raw <- read.delim("./data/SGBS_line_raw.txt")
data_SGBS_t_z <- read.delim("./data/SGBS_line_norm.txt")


# Jiawei add

FANTOM <- readRDS("./data/FANTOM.RDS")
FANTOM$classfication[FANTOM$ID %in% c("adipose, donor1", "adipose, donor2", "adipose, donor3", "adipose, donor4")] <- "other tissue"
FANTOM$classfication[FANTOM$classfication %in% c("muture adipocyte", "NiGa adipocyte")] <- "adipocyte"

adipogenesis <- readRDS("./data/adipogenesis.RDS")
adipogenesis <- adipogenesis[rownames(exprs(adipogenesis))[rowSums(exprs(adipogenesis))>0],1:48]

SVF <- readRDS("./data/SVF.RDS")

novo <- readRDS("./data/novoarray.RDS")
novo$Celltype[novo$Celltype=="Progenitor"] <- "progenitor"
novo$Celltype[novo$Celltype=="CD14+ Myeloid"] <- "total macrophage"
novo$Celltype[novo$Celltype=="Total T-cells"] <- "total T cell"
novo$Celltype[novo$Celltype=="CD4+ T-cells"] <- "CD4+ T cell"
novo$Celltype[novo$Celltype=="CD8+ T-cells"] <- "CD8+ T cell"
novo$Celltype[novo$Celltype=="M1 Macrophage"] <- "M1 macrophage"
novo$Celltype[novo$Celltype=="M2 Macrophage"] <- "M2 macrophage"
novo$Celltype[novo$Celltype=="Mature adipocytes"] <- "mature adipocyte"


data_niga_heat_temp <- exprs(adipogenesis)
data_niga_heat <- data.frame()
for (i in unique(gsub("_biol_rep[0-9]+$", "", colnames(data_niga_heat_temp)))) {
  # grep(i, colnames(data_niga_heat_temp), value = TRUE)
  if (length(data_niga_heat)==0) {
    data_niga_heat <- data.frame(rowMeans(data_niga_heat_temp[,grep(i, colnames(data_niga_heat_temp), value = TRUE)]))
  } else {
    data_niga_heat <- cbind(data_niga_heat,data.frame(rowMeans(data_niga_heat_temp[,grep(i, colnames(data_niga_heat_temp), value = TRUE)])))
  }
}
data_niga_heat <- as.matrix(data_niga_heat)
colnames(data_niga_heat) <- colnames(read.delim("./data/TimeCourse_wide.txt", row.names = 1)) %>% gsub(pattern = "X",replacement = "")
rm(data_niga_heat_temp)


data_SVF_heat_temp <- exprs(SVF)
data_SVF_heat <- data.frame()
for (i in unique(gsub("_donor[0-9]+$", "", colnames(data_SVF_heat_temp)))) {
  # grep(i, colnames(data_SVF_heat_temp), value = TRUE)
  if (length(data_SVF_heat)==0) {
    data_SVF_heat <- data.frame(rowMeans(data_SVF_heat_temp[,grep(i, colnames(data_SVF_heat_temp), value = TRUE)]))
  } else {
    data_SVF_heat <- cbind(data_SVF_heat,data.frame(rowMeans(data_SVF_heat_temp[,grep(i, colnames(data_SVF_heat_temp), value = TRUE)])))
  }
}
data_SVF_heat <- as.matrix(data_SVF_heat)
colnames(data_SVF_heat) <- colnames(read.delim("./data/TimeCourse_wideSVF.txt", row.names = 1)) %>% gsub(pattern = "X",replacement = "")
rm(data_SVF_heat_temp)


data_novo_heat <- data.frame()
for (i in unique(novo$Celltype)) {
  if (length(data_novo_heat)==0) {
    data_novo_heat <- data.frame(rowMeans(exprs(novo)[,i==novo$Celltype]))
  } else {
    data_novo_heat <- cbind(data_novo_heat, data.frame(rowMeans(exprs(novo)[,i==novo$Celltype])))
  }
}
# data_novo_heat$gene <- rownames(data_novo_heat)
data_novo_heat <- as.matrix(data_novo_heat)
colnames(data_novo_heat) <- unique(novo$Celltype)


data_WATfraction_prot_heat <- readRDS("./data/WAT_fraction_protromics.RDS")  
data_WATfraction_prot_heat_normalization <- readRDS("./data/WAT_fraction_protromics_normalization.RDS")  

#colors
#heatmap green to purple
PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#e2c744", "#1a4659", "#ba5c28", "#4f736a", "#f0cbbf", "#edf2f2" ,"#f2ecde", "#adbec4", "#c3b6b6" , "#8ca5a5", "#b2dfda", "#f59b7c", "#7dbfb3" ,"#939a64")
portal_col_ex <- colorRampPalette(portalcol2)
farben <- list("default"= portalcol,"PRGn" = PRGn, "RdBu" = RdBu, "viridis"=virid)

mode_in <- setNames(c("raw","normalized"),c("mean abundance", "normalized mean abundance")) 
