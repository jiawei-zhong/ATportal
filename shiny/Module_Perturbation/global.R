library(shiny)
library(viridis)
library(RColorBrewer)
library(reshape)
#cell stress
data_hyp <- read.delim("./data/hypoxia_toptable.txt", row.names = 1)
data_hyp_heat <- read.delim("./data/hypoxia_heat.txt", row.names = 1)
data_tnf <- read.delim("./data/tnf_tt_portal.txt", row.names = 1)
data_tnf_heat <- read.delim("./data/TNF_mat_portal.txt", row.names = 1)
#insulin cage
data_ins <- read.delim("./data/NEFA_0y_ctrls_pairedDESeq_2mioThreshold_noMale_noOB.txt", sep="\t")
data_ins_heat <- read.delim("./data/NEFA_0y_rlogTable_BaseMeanFilt_pairedDESeq_2mioT_noMale_noOB.txt", sep="\t", row.names = 1)
#knock down/ knock out
data_gls <- read.delim("./data/gls_tt_portal.txt")
data_gls_heat <- read.delim("./data/GLS_mat_portal.txt", row.names = 1)
data_cd248 <- read.delim("./data/CD248_toptable.txt")
data_cd248_heat <-read.delim("./data/CD248_heat.txt", row.names = 1)
data_aqp7 <- read.delim("./data/AQP7_toptable.txt")
data_aqp7_heat <- read.delim("./data/AQP7_heat.txt", row.names = 1)
data_c14 <- read.delim("./data/C14_toptable.txt")
data_c14_heat <- read.delim("./data/C14_heat.txt", row.names = 1)
ko_inputs <- list("GLS", "CD248","C14orf180","AQP7")
ko_inputs <-  ko_inputs[order(names(setNames(ko_inputs, ko_inputs)))]
#stimulation SGBS
data_stim <- readRDS("./data/stimuli_sgbs_tts.rds")
data_stim_heat <- read.delim("./data/stimuli_sgbs_matrix.txt", row.names = 1)
data_stim_group <- read.delim("./data/stimuli_sgbs_pheno.txt", row.names = 1)
#inflamed AT 
data_infl <- readRDS("./data/inflamed_tts.rds")
data_infl_heat <- read.delim("./data/inflamed_matrix.txt", row.names = 1)
data_infl_group <- read.delim("./data/inflamed_pheno.txt", row.names = 1)

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
