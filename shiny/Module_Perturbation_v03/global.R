library(shiny)
library(RColorBrewer)
library(reshape)
library(viridis)
library(shinycustomloader)
data_hyp <- readRDS("./data/hypoxia_eset.rds")
data_hyp_tt <- read.delim("./data/hypoxia_toptable.txt", row.names=1)
data_hyp <- exprs(data_hyp)


data_tnf_tt <- read.delim("./data/tnf_tt_portal.txt", row.names = 1)
data_tnf <- read.delim("./data/TNF_mat_portal.txt", row.names = 1)
data_gls_tt <- read.delim("./data/gls_tt_portal.txt", row.names = 1)
data_gls <- read.delim("./data/GLS_mat_portal.txt", row.names = 1)


#colors
#heatmap green to purple
PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
farben <- list("PRGn" = PRGn, "RdBu" = RdBu, "viridis"=virid)
