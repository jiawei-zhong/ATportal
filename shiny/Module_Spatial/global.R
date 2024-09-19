
# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))

discrete_color <- c("T, NK & NKT" = "#8ca5a5", "B" = "#7dbfb3", "mono. & macro." = "#4f736a", "mast" = "#1a4659", "FAPs" = "#f59b7c", "adipocytes" = "#ba5c28", "vascular" = "#e2c744")

discrete_color_gradient <- colorRampPalette(discrete_color)


trait_vector <- c("BMI" = "BMI", "HOMA" = "HOMA", "Age" = "Age", "WHR" = "WHR", "Waist" = "Waist", "Hip" = "Hip", "Glucose" = "Glucose", "Insulin" = "Insulin", "LEP Protein" = "LEP_protein", "TG" = "TG", "Chol" = "Chol", 
                  "HDL" = "HDL", "LDL" = "LDL", "CRP" = "CRP", "Hba1c" = "Hba1c", "TNF Protein" = "TNF_protein", "MCP1 Protein" = "MCP1_protein", "cell_vol" = "cell_vol", "basal_TG" = "basal_TG", "iso_TG" = "iso_TG", "iso_basal" = "iso_basal")

## Clinical association ----

# cohort_forest <- c("DEOSH"="DEOSHeset_Baseline", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF sc"="EMIFeset_sc", "Krieg et al. sc"="Krieg_et_al_sc", "Keller et al. sc"="Keller_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset","METSIM 770"="METSIM_770eset","METSIM 434"="METSIM_434eset","METSIM 200"="METSIM_200eset","EMIF om"="EMIFeset_om", "Krieg et al. om"="Krieg_et_al_om", "Keller et al. om"="Keller_et_al_om")
cohort_name <- c("Kerr, A. (2020) Baseline"="DEOSHeset_Baseline",            #  sc  phenotype
                 "Petrus, P. (2018) Baseline"="POeset_Baseline",             #  sc  phenotype 
                 "Arner, P. (2018)"="SOWOTeset",                             #  sc  phenotype 
                 "Keller, M. (2017) om"="Keller_et_al_om",                   #  om  phenotype  sex
                 "Keller, M. (2017) sc"="Keller_et_al_sc",                   #  sc  phenotype  sex
                 "Arner, E. (2012)"="RIKENeset",                             #  sc  phenotype
                 "Stančáková, A. (2012)"="METSIM_204eset",                   #  sc  phenotype
                 "Raulerson, C. (2019)"="METSIM_434eset",                    #  sc  phenotype
                 "Civelek, M. (2017)"="METSIM_770eset",                      #  sc  phenotype
                 "Krieg, L. (2021) sc"="Krieg_et_al_sc",                     #  sc  phenotype  sex
                 "Krieg, L. (2021) om"="Krieg_et_al_om",                     #  om  phenotype  sex
                 "Arner, P. (2016) om"="EMIFeset_om",                        #  om  phenotype
                 "Arner, P. (2016) sc"="EMIFeset_sc",                        #  sc  phenotype
                 "Imbert, A. (2022)"="GSE141221_diogenes1",                  #  sc  phenotype  sex
                 "Armenise, C. (2017)"="GSE95640_diogenes2",                 #  sc  phenotype  sex
                 "Winnier, DA. (2015)"="GSE64567eset",                       #  sc  phenotype  sex
                 "Keller, P. (2011)"="GSE27949eset",                         #  sc  phenotype
                 "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",        #  sc  phenotype  sex
                 "Barberio, MD. (2019)"="GSE88837eset",                      #  om  phenotype       ethnicity
                 "Sharma, NK. (2016)"="GSE95674eset",                        #  sc             sex
                 "GTEx microarray"="GTEx_microarray_sc",                     #  sc             sex
                 "Lonsdale, J. (2013) om"="GTEx_v4_om",                      #  om             sex
                 "Lonsdale, J. (2013) sc"="GTEx_v4_sc",                      #  sc             sex
                 "Aguet, F. (2017) om"="GTEx_v6_om",                         #  om             sex
                 "Aguet, F. (2017) sc"="GTEx_v6_sc",                         #  sc             sex
                 "GTEx v7 om"="GTEx_v7_om",                                  #  om             sex
                 "GTEx v7 sc"="GTEx_v7_sc",                                  #  sc             sex
                 "Aguet, F. (2020) om"="GTEx_v8_om",                         #  om             sex
                 "Aguet, F. (2020) sc"="GTEx_v8_sc",                         #  sc             sex
                 "Drong, A. (2013)"="E_MTAB_54",                             #  sc             sex
                 "Das, SK. (2015)"="GSE65221eset",                           #  sc             sex  ethnicity
                 "Naukkarinen, J. (2013)"="E_MTAB_1895",                     #  sc             sex
                 "Nookaew, I. (2013)"="GSE27916eset"                         #  sc             sex
)



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