# Load packages ----
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))




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



# snRNAseq ----

META_all <- readRDS("./data/META_all.sub.RDS")
META_all_sc <- readRDS("./data/META_all.s.sub.RDS")
META_all_om <- readRDS("./data/META_all.o.sub.RDS")
META_all_pvat <- readRDS("./data/META_all.p.sub.RDS")


temp <- as.character(Idents(META_all))
temp[temp=="Lymphoid"] <- "T, NK & NKT"
temp[temp=="B_cell"] <- "B"
temp[temp=="Myeloid"] <- "mono. & macro."
temp[temp=="Mast_cell"] <- "mast"
temp[temp=="Adipocyte"] <- "adipocytes"
temp[temp=="Vascular"] <- "vascular"
temp <- factor(temp,levels = c("T, NK & NKT", "B", "mono. & macro.", "mast", "FAPs", "adipocytes", "vascular"))
Idents(META_all) <- temp
rm(temp)


temp <- as.character(Idents(META_all_sc))
temp[temp=="Lymphoid"] <- "T, NK & NKT"
temp[temp=="B_cell"] <- "B"
temp[temp=="Myeloid"] <- "mono. & macro."
temp[temp=="Mast_cell"] <- "mast"
temp[temp=="Adipocyte"] <- "adipocytes"
temp[temp=="Vascular"] <- "vascular"
temp <- factor(temp,levels = c("T, NK & NKT", "B", "mono. & macro.", "mast", "FAPs", "adipocytes", "vascular"))
Idents(META_all_sc) <- temp
rm(temp)


temp <- as.character(Idents(META_all_om))
temp[temp=="Lymphoid"] <- "T, NK & NKT"
temp[temp=="B_cell"] <- "B"
temp[temp=="Myeloid"] <- "mono. & macro."
temp[temp=="Mast_cell"] <- "mast"
temp[temp=="Adipocyte"] <- "adipocytes"
temp[temp=="Vascular"] <- "vascular"
temp <- factor(temp,levels = c("T, NK & NKT", "B", "mono. & macro.", "mast", "FAPs", "adipocytes", "vascular"))
Idents(META_all_om) <- temp
rm(temp)


temp <- as.character(Idents(META_all_pvat))
temp[temp=="Lymphoid"] <- "T, NK & NKT"
# temp[temp=="B_cell"] <- "B"
temp[temp=="Myeloid"] <- "mono. & macro."
# temp[temp=="Mast_cell"] <- "mast"
temp[temp=="Adipocyte"] <- "adipocytes"
temp[temp=="Vascular"] <- "vascular"
temp <- factor(temp,levels = c("T, NK & NKT", "mono. & macro.", "FAPs", "adipocytes", "vascular"))
Idents(META_all_pvat) <- temp
rm(temp)





WAT_all_Lymphoid <- readRDS("./data/WAT_all_Lymphoid.RDS")
WAT_all_Myeloid <- readRDS( "./data/WAT_all_Myeloid.RDS")
WAT_all_Vascular <- readRDS("./data/WAT_all_Vascular.RDS")
WAT_all_B <- readRDS("./data/WAT_all_B.RDS")

WAT_sc_Lymphoid <- readRDS("./data/WAT_sc_Lymphoid.RDS")
WAT_sc_Myeloid <- readRDS( "./data/WAT_sc_Myeloid.RDS")
WAT_sc_FAPs <- readRDS(    "./data/WAT_sc_FAPs.RDS")
WAT_sc_Vascular <- readRDS("./data/WAT_sc_Vascular.RDS")
WAT_sc_B <- readRDS("./data/WAT_sc_B.RDS")

WAT_om_Lymphoid <- readRDS("./data/WAT_om_Lymphoid.RDS")
WAT_om_Myeloid <- readRDS( "./data/WAT_om_Myeloid.RDS")
WAT_om_FAPs <- readRDS(    "./data/WAT_om_FAPs.RDS")
WAT_om_Vascular <- readRDS("./data/WAT_om_Vascular.RDS")
WAT_om_B <- readRDS("./data/WAT_om_B.RDS")

WAT_pvat_Lymphoid <- readRDS("./data/WAT_pvat_Lymphoid.RDS")
WAT_pvat_Myeloid <- readRDS( "./data/WAT_pvat_Myeloid.RDS")
WAT_pvat_FAPs <- readRDS(    "./data/WAT_pvat_FAPs.RDS")
WAT_pvat_Vascular <- readRDS("./data/WAT_pvat_Vascular.RDS")

META_all_marker <- readRDS("./data/META_all_marker.RDS")
META_all_marker <- META_all_marker[META_all_marker$avg_log2FC>0.25 & (META_all_marker$pct.1>0.1 | META_all_marker$pct.2>0.1) & META_all_marker$p_val_adj<0.05,]
META_all_marker$cluster <- as.character(META_all_marker$cluster)
META_all_marker$cluster[META_all_marker$cluster=="Lymphoid"] <- "T, NK & NKT"
META_all_marker$cluster[META_all_marker$cluster=="B_cell"] <- "B"
META_all_marker$cluster[META_all_marker$cluster=="Myeloid"] <- "mono. & macro."
META_all_marker$cluster[META_all_marker$cluster=="Mast_cell"] <- "mast"
META_all_marker$cluster[META_all_marker$cluster=="Adipocyte"] <- "adipocytes"
META_all_marker$cluster[META_all_marker$cluster=="Vascular"] <- "vascular"


META_all_sc_marker <- readRDS("./data/META_all.s_marker.RDS")
META_all_sc_marker <- META_all_sc_marker[META_all_sc_marker$avg_log2FC>0.25 & (META_all_sc_marker$pct.1>0.1 | META_all_sc_marker$pct.2>0.1) & META_all_sc_marker$p_val_adj<0.05,]
META_all_sc_marker$cluster <- as.character(META_all_sc_marker$cluster)
META_all_sc_marker$cluster[META_all_sc_marker$cluster=="Lymphoid"] <- "T, NK & NKT"
META_all_sc_marker$cluster[META_all_sc_marker$cluster=="B_cell"] <- "B"
META_all_sc_marker$cluster[META_all_sc_marker$cluster=="Myeloid"] <- "mono. & macro."
META_all_sc_marker$cluster[META_all_sc_marker$cluster=="Mast_cell"] <- "mast"
META_all_sc_marker$cluster[META_all_sc_marker$cluster=="Adipocyte"] <- "adipocytes"
META_all_sc_marker$cluster[META_all_sc_marker$cluster=="Vascular"] <- "vascular"


META_all_om_marker <- readRDS("./data/META_all.o_marker.RDS")
META_all_om_marker <- META_all_om_marker[META_all_om_marker$avg_log2FC>0.25 & (META_all_om_marker$pct.1>0.1 | META_all_om_marker$pct.2>0.1) & META_all_om_marker$p_val_adj<0.05,]
META_all_om_marker$cluster <- as.character(META_all_om_marker$cluster)
META_all_om_marker$cluster[META_all_om_marker$cluster=="Lymphoid"] <- "T, NK & NKT"
META_all_om_marker$cluster[META_all_om_marker$cluster=="B_cell"] <- "B"
META_all_om_marker$cluster[META_all_om_marker$cluster=="Myeloid"] <- "mono. & macro."
META_all_om_marker$cluster[META_all_om_marker$cluster=="Mast_cell"] <- "mast"
META_all_om_marker$cluster[META_all_om_marker$cluster=="Adipocyte"] <- "adipocytes"
META_all_om_marker$cluster[META_all_om_marker$cluster=="Vascular"] <- "vascular"


META_all_pvat_marker <- readRDS("./data/META_all.p_marker.RDS")
META_all_pvat_marker <- META_all_pvat_marker[META_all_pvat_marker$avg_log2FC>0.25 & (META_all_pvat_marker$pct.1>0.1 | META_all_pvat_marker$pct.2>0.1) & META_all_pvat_marker$p_val_adj<0.05,]
META_all_pvat_marker$cluster <- as.character(META_all_pvat_marker$cluster)
META_all_pvat_marker$cluster[META_all_pvat_marker$cluster=="Lymphoid"] <- "T, NK & NKT"
META_all_pvat_marker$cluster[META_all_pvat_marker$cluster=="Myeloid"] <- "mono. & macro."
META_all_pvat_marker$cluster[META_all_pvat_marker$cluster=="Adipocyte"] <- "adipocytes"
META_all_pvat_marker$cluster[META_all_pvat_marker$cluster=="Vascular"] <- "vascular"


WAT_all_B_marker <- readRDS("./data/WAT_all_B_marker.RDS")
WAT_all_B_marker <- WAT_all_B_marker[WAT_all_B_marker$avg_log2FC>0.25 & (WAT_all_B_marker$pct.1>0.1 | WAT_all_B_marker$pct.2>0.1) & WAT_all_B_marker$p_val_adj<0.05,]
WAT_all_Lymphoid_marker <- readRDS("./data/WAT_all_Lymphoid_marker.RDS")
WAT_all_Lymphoid_marker <- WAT_all_Lymphoid_marker[WAT_all_Lymphoid_marker$avg_log2FC>0.25 & (WAT_all_Lymphoid_marker$pct.1>0.1 | WAT_all_Lymphoid_marker$pct.2>0.1) & WAT_all_Lymphoid_marker$p_val_adj<0.05,]
WAT_all_Myeloid_marker <- readRDS("./data/WAT_all_Myeloid_marker.RDS")
WAT_all_Myeloid_marker <- WAT_all_Myeloid_marker[WAT_all_Myeloid_marker$avg_log2FC>0.25 & (WAT_all_Myeloid_marker$pct.1>0.1 | WAT_all_Myeloid_marker$pct.2>0.1) & WAT_all_Myeloid_marker$p_val_adj<0.05,]
WAT_all_Vascular_marker <- readRDS("./data/WAT_all_Vascular_marker.RDS")
WAT_all_Vascular_marker <- WAT_all_Vascular_marker[WAT_all_Vascular_marker$avg_log2FC>0.25 & (WAT_all_Vascular_marker$pct.1>0.1 | WAT_all_Vascular_marker$pct.2>0.1) & WAT_all_Vascular_marker$p_val_adj<0.05,]

WAT_sc_B_marker <- readRDS("./data/WAT_sc_B_marker.RDS")
WAT_sc_B_marker <- WAT_sc_B_marker[WAT_sc_B_marker$avg_log2FC>0.25 & (WAT_sc_B_marker$pct.1>0.1 | WAT_sc_B_marker$pct.2>0.1) & WAT_sc_B_marker$p_val_adj<0.05,]
WAT_sc_FAPs_marker <- readRDS("./data/WAT_sc_FAPs_marker.RDS")
WAT_sc_FAPs_marker <- WAT_sc_FAPs_marker[WAT_sc_FAPs_marker$avg_log2FC>0.25 & (WAT_sc_FAPs_marker$pct.1>0.1 | WAT_sc_FAPs_marker$pct.2>0.1) & WAT_sc_FAPs_marker$p_val_adj<0.05,]
WAT_sc_Lymphoid_marker <- readRDS("./data/WAT_sc_Lymphoid_marker.RDS")
WAT_sc_Lymphoid_marker <- WAT_sc_Lymphoid_marker[WAT_sc_Lymphoid_marker$avg_log2FC>0.25 & (WAT_sc_Lymphoid_marker$pct.1>0.1 | WAT_sc_Lymphoid_marker$pct.2>0.1) & WAT_sc_Lymphoid_marker$p_val_adj<0.05,]
WAT_sc_Myeloid_marker <- readRDS("./data/WAT_sc_Myeloid_marker.RDS")
WAT_sc_Myeloid_marker <- WAT_sc_Myeloid_marker[WAT_sc_Myeloid_marker$avg_log2FC>0.25 & (WAT_sc_Myeloid_marker$pct.1>0.1 | WAT_sc_Myeloid_marker$pct.2>0.1) & WAT_sc_Myeloid_marker$p_val_adj<0.05,]
WAT_sc_Vascular_marker <- readRDS("./data/WAT_sc_Vascular_marker.RDS")
WAT_sc_Vascular_marker <- WAT_sc_Vascular_marker[WAT_sc_Vascular_marker$avg_log2FC>0.25 & (WAT_sc_Vascular_marker$pct.1>0.1 | WAT_sc_Vascular_marker$pct.2>0.1) & WAT_sc_Vascular_marker$p_val_adj<0.05,]

WAT_om_B_marker <- readRDS("./data/WAT_om_B_marker.RDS")
WAT_om_B_marker <- WAT_om_B_marker[WAT_om_B_marker$avg_log2FC>0.25 & (WAT_om_B_marker$pct.1>0.1 | WAT_om_B_marker$pct.2>0.1) & WAT_om_B_marker$p_val_adj<0.05,]
WAT_om_FAPs_marker <- readRDS("./data/WAT_om_FAPs_marker.RDS")
WAT_om_FAPs_marker <- WAT_om_FAPs_marker[WAT_om_FAPs_marker$avg_log2FC>0.25 & (WAT_om_FAPs_marker$pct.1>0.1 | WAT_om_FAPs_marker$pct.2>0.1) & WAT_om_FAPs_marker$p_val_adj<0.05,]
WAT_om_Lymphoid_marker <- readRDS("./data/WAT_om_Lymphoid_marker.RDS")
WAT_om_Lymphoid_marker <- WAT_om_Lymphoid_marker[WAT_om_Lymphoid_marker$avg_log2FC>0.25 & (WAT_om_Lymphoid_marker$pct.1>0.1 | WAT_om_Lymphoid_marker$pct.2>0.1) & WAT_om_Lymphoid_marker$p_val_adj<0.05,]
WAT_om_Myeloid_marker <- readRDS("./data/WAT_om_Myeloid_marker.RDS")
WAT_om_Myeloid_marker <- WAT_om_Myeloid_marker[WAT_om_Myeloid_marker$avg_log2FC>0.25 & (WAT_om_Myeloid_marker$pct.1>0.1 | WAT_om_Myeloid_marker$pct.2>0.1) & WAT_om_Myeloid_marker$p_val_adj<0.05,]
WAT_om_Vascular_marker <- readRDS("./data/WAT_om_Vascular_marker.RDS")
WAT_om_Vascular_marker <- WAT_om_Vascular_marker[WAT_om_Vascular_marker$avg_log2FC>0.25 & (WAT_om_Vascular_marker$pct.1>0.1 | WAT_om_Vascular_marker$pct.2>0.1) & WAT_om_Vascular_marker$p_val_adj<0.05,]

WAT_pvat_FAPs_marker <- readRDS("./data/WAT_pvat_FAPs_marker.RDS")
WAT_pvat_FAPs_marker <- WAT_pvat_FAPs_marker[WAT_pvat_FAPs_marker$avg_log2FC>0.25 & (WAT_pvat_FAPs_marker$pct.1>0.1 | WAT_pvat_FAPs_marker$pct.2>0.1) & WAT_pvat_FAPs_marker$p_val_adj<0.05,]
WAT_pvat_Lymphoid_marker <- readRDS("./data/WAT_pvat_Lymphoid_marker.RDS")
WAT_pvat_Lymphoid_marker <- WAT_pvat_Lymphoid_marker[WAT_pvat_Lymphoid_marker$avg_log2FC>0.25 & (WAT_pvat_Lymphoid_marker$pct.1>0.1 | WAT_pvat_Lymphoid_marker$pct.2>0.1) & WAT_pvat_Lymphoid_marker$p_val_adj<0.05,]
WAT_pvat_Myeloid_marker <- readRDS("./data/WAT_pvat_Myeloid_marker.RDS")
WAT_pvat_Myeloid_marker <- WAT_pvat_Myeloid_marker[WAT_pvat_Myeloid_marker$avg_log2FC>0.25 & (WAT_pvat_Myeloid_marker$pct.1>0.1 | WAT_pvat_Myeloid_marker$pct.2>0.1) & WAT_pvat_Myeloid_marker$p_val_adj<0.05,]
WAT_pvat_Vascular_marker <- readRDS("./data/WAT_pvat_Vascular_marker.RDS")
WAT_pvat_Vascular_marker <- WAT_pvat_Vascular_marker[WAT_pvat_Vascular_marker$avg_log2FC>0.25 & (WAT_pvat_Vascular_marker$pct.1>0.1 | WAT_pvat_Vascular_marker$pct.2>0.1) & WAT_pvat_Vascular_marker$p_val_adj<0.05,]


for (i in grep("marker$",ls(),value = T)) {
  eval(parse(text=paste0(i,'$p_val <- signif(',i,'$p_val,4)')))
  eval(parse(text=paste0(i,'$p_val_adj <- signif(',i,'$p_val_adj,4)')))
  eval(parse(text=paste0(i,'$avg_log2FC <- signif(',i,'$avg_log2FC,4)')))
}


gene_all <- sort(rownames(META_all))
gene_all <- c(gene_all[grep("^AC\\d{6}|^AD\\d{6}|^AE\\d{6}|^AF\\d{6}|^AJ\\d{6}|^AL\\d{6}|^AP\\d{6}|^CATG\\d{6}",gene_all,invert = T)],
              gene_all[grep("^AC\\d{6}|^AD\\d{6}|^AE\\d{6}|^AF\\d{6}|^AJ\\d{6}|^AL\\d{6}|^AP\\d{6}|^CATG\\d{6}",gene_all,invert = F)])
