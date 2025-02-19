# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(dplyr))


trait_vector <- c("BMI" = "BMI", "HOMA" = "HOMA", "Age" = "Age", "WHR" = "WHR", "Waist" = "Waist", "Hip" = "Hip", "Glucose" = "Glucose", "Insulin" = "Insulin", "LEP Protein" = "LEP_protein", "TG" = "TG", "Chol" = "Chol", 
                  "HDL" = "HDL", "LDL" = "LDL", "CRP" = "CRP", "Hba1c" = "Hba1c", "TNF Protein" = "TNF_protein", "MCP1 Protein" = "MCP1_protein", "cell_vol" = "cell_vol", "basal_TG" = "basal_TG", "iso_TG" = "iso_TG", "iso_basal" = "iso_basal")


# cohort_forest <- c("DEOSH"="DEOSHeset_Baseline", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF sc"="EMIFeset_sc", "Krieg et al. sc"="Krieg_et_al_sc", "Keller et al. sc"="Keller_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset","METSIM 770"="METSIM_770eset","METSIM 434"="METSIM_434eset","METSIM 200"="METSIM_200eset","EMIF om"="EMIFeset_om", "Krieg et al. om"="Krieg_et_al_om", "Keller et al. om"="Keller_et_al_om")
cohort_name <- c("Keller, M. (2017)"="Keller_et_al",
                 "Krieg, L. (2021)"="Krieg_et_al",
                 "Arner, P. (2016)"="EMIFeset",
                 "Schleinitz, D. (2020)"="Schleinitz_et_al",
                 "Mazaki-Tovi, S. (2016)"="Mazaki_et_al",
                 "Hoggard, N. (2012)"="Hoggard_et_al",
                 "MacLaren, RE. (2010)"="GSE15524eset",
                 "Salcedo-Tacuma, D. (2022)"="GSE188799eset",
                 "Hardy, OT. (2011)"="GSE20950eset",
                 "Du Plessis, J. (2015)"="GSE58979eset",
                 "Hruska, P. (2021)"="Hruska_2021",
                 "Hruska, P. (2023)"="Hruska_2023",
                 "Krieg, L. (2021)"="Krieg_et_al_proteomics"
                 )

citation <- readRDS("./data/citation.RDS")

## Depot Comparisons ----

for (i in cohort_name) {
  temp <- readRDS(paste0("./data/",i,".rds"))
  temp$Tissue[temp$Tissue=="sc"] <- "SAT Abdomen"
  temp$Tissue[temp$Tissue=="om"] <- "VAT Omentum"
  
  if (i=="Schleinitz_et_al") {
    rownames(temp) <- temp@featureData@data[["Genes"]]
    
    temp[,colnames(temp)[temp$Tissue %in% c("SAT Abdomen", "VAT Omentum")]]
    
    temp$ID <- gsub("_.*|-.*","",gsub("DS_|DS ","",temp$ID))
    # temp$ID <- temp$ID %>% strsplit("_") %>% lapply( function(x) x[1]) %>% unlist()
  }
  
  temp <- temp[,colnames(temp)[temp$ID %in% names(table(temp$ID)[table(temp$ID)==length(unique(temp$Tissue))])]]
  
  assign(x = i,value = temp)
}

exprs(GSE188799eset) <- log2(exprs(GSE188799eset)+1)


# Keller_et_al <- readRDS("./data/Keller_et_al.rds")
# Krieg_et_al <- readRDS("./data/Krieg_et_al.rds")
# EMIFeset <- readRDS("./data/EMIFeset.rds")
# Schleinitz_et_al <- readRDS("./data/Schleinitz_et_al.rds")
# Mazaki_et_al <- readRDS("./data/Mazaki_et_al.rds")
# Hoggard_et_al <- readRDS("./data/Hoggard_et_al.rds")

# 
# Hruska_2021 <- readRDS("./data/Hruska_2021.rds")
# Hruska_2023 <- readRDS("./data/Hruska_2023.rds")
# Krieg_et_al_proteomics <- readRDS("./data/Krieg_et_al_proteomics.rds")


PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
viridis <- viridis(51)
portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")
portal_col_ex <- colorRampPalette(portalcol2)