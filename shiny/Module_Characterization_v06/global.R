library(shiny)
library(viridis)
library(RColorBrewer)
library(reshape)
data_line <- read.delim("./data/TimeCourse.txt")
data_niga_heat <- read.delim("./data/TimeCourse_wide.txt", row.names = 1)
data_SVF <- read.delim("./data/TimeCourseSVF.txt")
data_SVF$timepoint <- factor(data_SVF$timepoint, levels = c("day4"   ,    "day8"   ,    "day12"  ,    "adipocytes"))
data_SVF_heat <- read.delim("./data/TimeCourse_wideSVF.txt", row.names = 1)
data_novo <- read.delim("./data/novoarray.txt")
data_novo_heat <- read.delim("./data/novo_heat.txt", row.names = 1)
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

SVF <- readRDS("./data/SVF.RDS")


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
