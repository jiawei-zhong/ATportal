# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))


trait_vector <- c("BMI" = "BMI", "HOMA" = "HOMA", "Age" = "Age", "WHR" = "WHR", "Waist" = "Waist", "Hip" = "Hip", "Glucose" = "Glucose", "Insulin" = "Insulin", "LEP Protein" = "LEP_protein", "TG" = "TG", "Chol" = "Chol", 
                  "HDL" = "HDL", "LDL" = "LDL", "CRP" = "CRP", "Hba1c" = "Hba1c", "TNF Protein" = "TNF_protein", "MCP1 Protein" = "MCP1_protein", "cell_vol" = "cell_vol", "basal_TG" = "basal_TG", "iso_TG" = "iso_TG", "iso_basal" = "iso_basal")


# cohort_forest <- c("DEOSH"="DEOSHeset_Baseline", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF sc"="EMIFeset_sc", "Krieg et al. sc"="Krieg_et_al_sc", "Keller et al. sc"="Keller_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset","METSIM 770"="METSIM_770eset","METSIM 434"="METSIM_434eset","METSIM 200"="METSIM_200eset","EMIF om"="EMIFeset_om", "Krieg et al. om"="Krieg_et_al_om", "Keller et al. om"="Keller_et_al_om")
cohort_name <- c("Keller, M. (2017)"="Keller_et_al",
                 "Krieg, L. (2021)"="Krieg_et_al",
                 "Arner, P. (2016)"="EMIFeset",
                 "Schleinitz, D. (2020)"="Schleinitz_et_al",
                 "Mazaki-Tovi, S. (2016)"="Mazaki_et_al",
                 "Hoggard, N. (2012)"="Hoggard_et_al"
                 )

citation <- readRDS("./data/citation.RDS")

## Depot Comparisons ----

Keller_et_al <- readRDS("./data/Keller_et_al.rds")
Krieg_et_al <- readRDS("./data/Krieg_et_al.rds")
EMIFeset <- readRDS("./data/EMIFeset.rds")
Schleinitz_et_al <- readRDS("./data/Schleinitz.rds")
Mazaki_et_al <- readRDS("./data/Mazaki.rds")
Hoggard_et_al <- readRDS("./data/Hoggard.rds")



