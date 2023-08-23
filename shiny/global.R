
# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))



## Clinical ----

clinical_measure_vector <- c("BMI" = "BMI", "HOMA" = "HOMA", "Age" = "Age", "WHR" = "WHR", "Waist" = "Waist", "Hip" = "Hip", "Glucose" = "Glucose", "Insulin" = "Insulin", "LEP Protein" = "LEP_protein", "TG" = "TG", "Chol" = "Chol", 
                             "HDL" = "HDL", "LDL" = "LDL", "CRP" = "CRP", "Hba1c" = "Hba1c", "TNF Protein" = "TNF_protein", "MCP1 Protein" = "MCP1_protein", "cell_vol" = "cell_vol", "basal_TG" = "basal_TG", "iso_TG" = "iso_TG", "iso_basal" = "iso_basal")

cohort_forest <- c("DEOSH"="DEOSHeset_Baseline", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF sc"="EMIFeset_sc", "Krieg et al. sc"="Krieg_et_al_sc", "Keller et al. sc"="Keller_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset","METSIM 770"="METSIM_770eset","METSIM 434"="METSIM_434eset","METSIM 200"="METSIM_200eset","EMIF om"="EMIFeset_om", "Krieg et al. om"="Krieg_et_al_om", "Keller et al. om"="Keller_et_al_om")

## Depots ----

data_GO_NES <- readRDS("./data/Depots/NES_mat.rds")
data_GO_P <- readRDS("./data/Depots/P_mat.rds")
data_GO_Size <- readRDS("./data/Depots/Size_mat.rds")


## Characterization ----

data_line <- read.delim("./data/Characterization/TimeCourse.txt")
data_niga_heat <- read.delim("./data/Characterization/TimeCourse_wide.txt", row.names = 1)
data_SVF <- read.delim("./data/Characterization/TimeCourseSVF.txt")
data_SVF$timepoint <- factor(data_SVF$timepoint, levels = c("day4"   ,    "day8"   ,    "day12"  ,    "adipocytes"))
data_SVF_heat <- read.delim("./data/Characterization/TimeCourse_wideSVF.txt", row.names = 1)
data_novo <- read.delim("./data/Characterization/novoarray.txt")
# data_novo_heat <- read.delim("./data/Characterization/novo_heat.txt", row.names = 1)
data_novo_heat <- read.delim("./data/Characterization/novo_heat_unique.txt", row.names = 1) # zjw edited
data_fantom <- read.delim("./data/Characterization/fantom.txt", row.names = 1)

# snRNAseq ----

META_all <- readRDS("./data/Singlecell/META_all.sub.RDS")
META_all_sc <- readRDS("./data/Singlecell/META_all.s.sub.RDS")
META_all_om <- readRDS("./data/Singlecell/META_all.o.sub.RDS")
META_all_pvat <- readRDS("./data/Singlecell/META_all.p.sub.RDS")

WAT_all_Lymphoid <- readRDS("./data/Singlecell/WAT_all_Lymphoid.RDS")
WAT_all_Myeloid <- readRDS( "./data/Singlecell/WAT_all_Myeloid.RDS")
WAT_all_Vascular <- readRDS("./data/Singlecell/WAT_all_Vascular.RDS")
WAT_all_B <- readRDS("./data/Singlecell/WAT_all_B.RDS")

WAT_sc_Lymphoid <- readRDS("./data/Singlecell/WAT_sc_Lymphoid.RDS")
WAT_sc_Myeloid <- readRDS( "./data/Singlecell/WAT_sc_Myeloid.RDS")
WAT_sc_FAPs <- readRDS(    "./data/Singlecell/WAT_sc_FAPs.RDS")
WAT_sc_Vascular <- readRDS("./data/Singlecell/WAT_sc_Vascular.RDS")
WAT_sc_B <- readRDS("./data/Singlecell/WAT_sc_B.RDS")

WAT_om_Lymphoid <- readRDS("./data/Singlecell/WAT_om_Lymphoid.RDS")
WAT_om_Myeloid <- readRDS( "./data/Singlecell/WAT_om_Myeloid.RDS")
WAT_om_FAPs <- readRDS(    "./data/Singlecell/WAT_om_FAPs.RDS")
WAT_om_Vascular <- readRDS("./data/Singlecell/WAT_om_Vascular.RDS")
WAT_om_B <- readRDS("./data/Singlecell/WAT_om_B.RDS")

WAT_pvat_Lymphoid <- readRDS("./data/Singlecell/WAT_pvat_Lymphoid.RDS")
WAT_pvat_Myeloid <- readRDS( "./data/Singlecell/WAT_pvat_Myeloid.RDS")
WAT_pvat_FAPs <- readRDS(    "./data/Singlecell/WAT_pvat_FAPs.RDS")
WAT_pvat_Vascular <- readRDS("./data/Singlecell/WAT_pvat_Vascular.RDS")

META_all_marker <- readRDS("./data/Singlecell/META_all_marker.RDS")
META_all_marker <- META_all_marker[META_all_marker$avg_log2FC>0.25 & (META_all_marker$pct.1>0.1 | META_all_marker$pct.2>0.1) & META_all_marker$p_val_adj<0.05,]
META_all_sc_marker <- readRDS("./data/Singlecell/META_all.s_marker.RDS")
META_all_sc_marker <- META_all_sc_marker[META_all_sc_marker$avg_log2FC>0.25 & (META_all_sc_marker$pct.1>0.1 | META_all_sc_marker$pct.2>0.1) & META_all_sc_marker$p_val_adj<0.05,]
META_all_om_marker <- readRDS("./data/Singlecell/META_all.o_marker.RDS")
META_all_om_marker <- META_all_om_marker[META_all_om_marker$avg_log2FC>0.25 & (META_all_om_marker$pct.1>0.1 | META_all_om_marker$pct.2>0.1) & META_all_om_marker$p_val_adj<0.05,]
META_all_pvat_marker <- readRDS("./data/Singlecell/META_all.p_marker.RDS")
META_all_pvat_marker <- META_all_pvat_marker[META_all_pvat_marker$avg_log2FC>0.25 & (META_all_pvat_marker$pct.1>0.1 | META_all_pvat_marker$pct.2>0.1) & META_all_pvat_marker$p_val_adj<0.05,]

WAT_all_B_marker <- readRDS("./data/Singlecell/WAT_all_B_marker.RDS")
WAT_all_B_marker <- WAT_all_B_marker[WAT_all_B_marker$avg_log2FC>0.25 & (WAT_all_B_marker$pct.1>0.1 | WAT_all_B_marker$pct.2>0.1) & WAT_all_B_marker$p_val_adj<0.05,]
WAT_all_Lymphoid_marker <- readRDS("./data/Singlecell/WAT_all_Lymphoid_marker.RDS")
WAT_all_Lymphoid_marker <- WAT_all_Lymphoid_marker[WAT_all_Lymphoid_marker$avg_log2FC>0.25 & (WAT_all_Lymphoid_marker$pct.1>0.1 | WAT_all_Lymphoid_marker$pct.2>0.1) & WAT_all_Lymphoid_marker$p_val_adj<0.05,]
WAT_all_Myeloid_marker <- readRDS("./data/Singlecell/WAT_all_Myeloid_marker.RDS")
WAT_all_Myeloid_marker <- WAT_all_Myeloid_marker[WAT_all_Myeloid_marker$avg_log2FC>0.25 & (WAT_all_Myeloid_marker$pct.1>0.1 | WAT_all_Myeloid_marker$pct.2>0.1) & WAT_all_Myeloid_marker$p_val_adj<0.05,]
WAT_all_Vascular_marker <- readRDS("./data/Singlecell/WAT_all_Vascular_marker.RDS")
WAT_all_Vascular_marker <- WAT_all_Vascular_marker[WAT_all_Vascular_marker$avg_log2FC>0.25 & (WAT_all_Vascular_marker$pct.1>0.1 | WAT_all_Vascular_marker$pct.2>0.1) & WAT_all_Vascular_marker$p_val_adj<0.05,]

WAT_sc_B_marker <- readRDS("./data/Singlecell/WAT_sc_B_marker.RDS")
WAT_sc_B_marker <- WAT_sc_B_marker[WAT_sc_B_marker$avg_log2FC>0.25 & (WAT_sc_B_marker$pct.1>0.1 | WAT_sc_B_marker$pct.2>0.1) & WAT_sc_B_marker$p_val_adj<0.05,]
WAT_sc_FAPs_marker <- readRDS("./data/Singlecell/WAT_sc_FAPs_marker.RDS")
WAT_sc_FAPs_marker <- WAT_sc_FAPs_marker[WAT_sc_FAPs_marker$avg_log2FC>0.25 & (WAT_sc_FAPs_marker$pct.1>0.1 | WAT_sc_FAPs_marker$pct.2>0.1) & WAT_sc_FAPs_marker$p_val_adj<0.05,]
WAT_sc_Lymphoid_marker <- readRDS("./data/Singlecell/WAT_sc_Lymphoid_marker.RDS")
WAT_sc_Lymphoid_marker <- WAT_sc_Lymphoid_marker[WAT_sc_Lymphoid_marker$avg_log2FC>0.25 & (WAT_sc_Lymphoid_marker$pct.1>0.1 | WAT_sc_Lymphoid_marker$pct.2>0.1) & WAT_sc_Lymphoid_marker$p_val_adj<0.05,]
WAT_sc_Myeloid_marker <- readRDS("./data/Singlecell/WAT_sc_Myeloid_marker.RDS")
WAT_sc_Myeloid_marker <- WAT_sc_Myeloid_marker[WAT_sc_Myeloid_marker$avg_log2FC>0.25 & (WAT_sc_Myeloid_marker$pct.1>0.1 | WAT_sc_Myeloid_marker$pct.2>0.1) & WAT_sc_Myeloid_marker$p_val_adj<0.05,]
WAT_sc_Vascular_marker <- readRDS("./data/Singlecell/WAT_sc_Vascular_marker.RDS")
WAT_sc_Vascular_marker <- WAT_sc_Vascular_marker[WAT_sc_Vascular_marker$avg_log2FC>0.25 & (WAT_sc_Vascular_marker$pct.1>0.1 | WAT_sc_Vascular_marker$pct.2>0.1) & WAT_sc_Vascular_marker$p_val_adj<0.05,]

WAT_om_B_marker <- readRDS("./data/Singlecell/WAT_om_B_marker.RDS")
WAT_om_B_marker <- WAT_om_B_marker[WAT_om_B_marker$avg_log2FC>0.25 & (WAT_om_B_marker$pct.1>0.1 | WAT_om_B_marker$pct.2>0.1) & WAT_om_B_marker$p_val_adj<0.05,]
WAT_om_FAPs_marker <- readRDS("./data/Singlecell/WAT_om_FAPs_marker.RDS")
WAT_om_FAPs_marker <- WAT_om_FAPs_marker[WAT_om_FAPs_marker$avg_log2FC>0.25 & (WAT_om_FAPs_marker$pct.1>0.1 | WAT_om_FAPs_marker$pct.2>0.1) & WAT_om_FAPs_marker$p_val_adj<0.05,]
WAT_om_Lymphoid_marker <- readRDS("./data/Singlecell/WAT_om_Lymphoid_marker.RDS")
WAT_om_Lymphoid_marker <- WAT_om_Lymphoid_marker[WAT_om_Lymphoid_marker$avg_log2FC>0.25 & (WAT_om_Lymphoid_marker$pct.1>0.1 | WAT_om_Lymphoid_marker$pct.2>0.1) & WAT_om_Lymphoid_marker$p_val_adj<0.05,]
WAT_om_Myeloid_marker <- readRDS("./data/Singlecell/WAT_om_Myeloid_marker.RDS")
WAT_om_Myeloid_marker <- WAT_om_Myeloid_marker[WAT_om_Myeloid_marker$avg_log2FC>0.25 & (WAT_om_Myeloid_marker$pct.1>0.1 | WAT_om_Myeloid_marker$pct.2>0.1) & WAT_om_Myeloid_marker$p_val_adj<0.05,]
WAT_om_Vascular_marker <- readRDS("./data/Singlecell/WAT_om_Vascular_marker.RDS")
WAT_om_Vascular_marker <- WAT_om_Vascular_marker[WAT_om_Vascular_marker$avg_log2FC>0.25 & (WAT_om_Vascular_marker$pct.1>0.1 | WAT_om_Vascular_marker$pct.2>0.1) & WAT_om_Vascular_marker$p_val_adj<0.05,]

WAT_pvat_FAPs_marker <- readRDS("./data/Singlecell/WAT_pvat_FAPs_marker.RDS")
WAT_pvat_FAPs_marker <- WAT_pvat_FAPs_marker[WAT_pvat_FAPs_marker$avg_log2FC>0.25 & (WAT_pvat_FAPs_marker$pct.1>0.1 | WAT_pvat_FAPs_marker$pct.2>0.1) & WAT_pvat_FAPs_marker$p_val_adj<0.05,]
WAT_pvat_Lymphoid_marker <- readRDS("./data/Singlecell/WAT_pvat_Lymphoid_marker.RDS")
WAT_pvat_Lymphoid_marker <- WAT_pvat_Lymphoid_marker[WAT_pvat_Lymphoid_marker$avg_log2FC>0.25 & (WAT_pvat_Lymphoid_marker$pct.1>0.1 | WAT_pvat_Lymphoid_marker$pct.2>0.1) & WAT_pvat_Lymphoid_marker$p_val_adj<0.05,]
WAT_pvat_Myeloid_marker <- readRDS("./data/Singlecell/WAT_pvat_Myeloid_marker.RDS")
WAT_pvat_Myeloid_marker <- WAT_pvat_Myeloid_marker[WAT_pvat_Myeloid_marker$avg_log2FC>0.25 & (WAT_pvat_Myeloid_marker$pct.1>0.1 | WAT_pvat_Myeloid_marker$pct.2>0.1) & WAT_pvat_Myeloid_marker$p_val_adj<0.05,]
WAT_pvat_Vascular_marker <- readRDS("./data/Singlecell/WAT_pvat_Vascular_marker.RDS")
WAT_pvat_Vascular_marker <- WAT_pvat_Vascular_marker[WAT_pvat_Vascular_marker$avg_log2FC>0.25 & (WAT_pvat_Vascular_marker$pct.1>0.1 | WAT_pvat_Vascular_marker$pct.2>0.1) & WAT_pvat_Vascular_marker$p_val_adj<0.05,]



## Spatial ----

Jesper_et_al_baseline <- readRDS("./data/Spatial/Jesper_et_al_baseline.RDS")
Jesper_et_al_insulin <- readRDS("./data/Spatial/Jesper_et_al_insulin.RDS")

Jesper_et_al_baseline_marker <- readRDS("./data/Spatial/Jesper_et_al_baseline_marker.RDS")
Jesper_et_al_baseline_marker <- Jesper_et_al_baseline_marker[Jesper_et_al_baseline_marker$avg_log2FC>0.25 & (Jesper_et_al_baseline_marker$pct.1>0.1 | Jesper_et_al_baseline_marker$pct.2>0.1) & Jesper_et_al_baseline_marker$p_val_adj<0.05,]
Jesper_et_al_insulin_marker <- readRDS("./data/Spatial/Jesper_et_al_insulin_marker.RDS")
Jesper_et_al_insulin_marker <- Jesper_et_al_insulin_marker[Jesper_et_al_insulin_marker$avg_log2FC>0.25 & (Jesper_et_al_insulin_marker$pct.1>0.1 | Jesper_et_al_insulin_marker$pct.2>0.1) & Jesper_et_al_insulin_marker$p_val_adj<0.05,]




for (i in grep("marker$",ls(),value = T)) {
  eval(parse(text=paste0(i,'$p_val <- signif(',i,'$p_val,4)')))
  eval(parse(text=paste0(i,'$p_val_adj <- signif(',i,'$p_val_adj,4)')))
  eval(parse(text=paste0(i,'$avg_log2FC <- signif(',i,'$avg_log2FC,4)')))
}



## Color ----
#heatmap green to purple
PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
farben <- list("PRGn" = PRGn, "RdBu" = RdBu, "viridis"=virid)

