# Load packages ----
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(extrafont))


## data ----

# trait_vector <- c("Body mass index" = "BMI", 
#                   "Homeostatic model of insulin resistance" = "HOMA", 
#                   "Age" = "Age", 
#                   "Waist-to-hip ratio" = "WHR", 
#                   "Waist circumference" = "Waist", 
#                   "Hip circumference" = "Hip", 
#                   "Glucose" = "Glucose", 
#                   "Insulin" = "Insulin", 
#                   "Triglyceride" = "TG", 
#                   "Cholesterol" = "Chol", 
#                   "High-density lipoprotein" = "HDL", 
#                   "Low-density lipoprotein" = "LDL", 
#                   "C-reactive protein" = "CRP", 
#                   "Hemoglobin A1C" = "Hba1c", 
#                   "Adipose LEP protein secretion" = "LEP_protein", 
#                   "Adipose TNFα protein secretion" = "TNF_protein", 
#                   "Adipose MCP1 protein secretion" = "MCP1_protein", 
#                   "Fat cell volume" = "cell_vol", 
#                   "basal_TG" = "basal_TG", 
#                   "iso_TG" = "iso_TG", 
#                   "iso_basal" = "iso_basal")


trait_vector <- c("BMI" = "BMI", 
                  "HOMA-IR" = "HOMA", 
                  "age" = "Age", 
                  "WHR" = "WHR", 
                  "waist" = "Waist", 
                  "hip" = "Hip", 
                  "circ-glucose" = "Glucose", 
                  "circ-insulin" = "Insulin", 
                  "circ-TG" = "TG", 
                  "circ-cholesterol" = "Chol", 
                  "circ-HDL" = "HDL", 
                  "circ-LDL" = "LDL", 
                  "circ-CRP" = "CRP", 
                  "HbA1c" = "Hba1c", 
                  "WAT LEP secretion" = "LEP_protein", 
                  "WAT TNF secretion" = "TNF_protein", 
                  "WAT MCP1 secretion" = "MCP1_protein", 
                  "fat cell volume" = "cell_vol", 
                  "basal lipolysis" = "basal_TG", 
                  "iso lipolysis" = "iso_TG", 
                  "iso/basal" = "iso_basal")



cohort_name <- c("Kerr, A. (2020) Baseline"="DEOSHeset_Baseline",                     #  sc  phenotype
                 "Petrus, P. (2018) Baseline"="POeset_Baseline",                      #  sc  phenotype
                 "Arner, P. (2018)"="SOWOTeset",                                      #  sc  phenotype
                 "Keller, M. (2017) om"="Keller_et_al_om",                            #  om  phenotype  sex
                 "Keller, M. (2017) sc"="Keller_et_al_sc",                            #  sc  phenotype  sex
                 "Arner, E. (2012)"="RIKENeset",                                      #  sc  phenotype
                 "Stančáková, A. (2012)"="METSIM_204eset",                            #  sc  phenotype
                 "Raulerson, C. (2019)"="METSIM_434eset",                             #  sc  phenotype
                 "Civelek, M. (2017)"="METSIM_770eset",                               #  sc  phenotype
                 "Krieg, L. (2021) sc"="Krieg_et_al_sc",                              #  sc  phenotype  sex
                 "Krieg, L. (2021) om"="Krieg_et_al_om",                              #  om  phenotype  sex
                 "Arner, P. (2016) om"="EMIFeset_om",                                 #  om  phenotype
                 "Arner, P. (2016) sc"="EMIFeset_sc",                                 #  sc  phenotype
                 "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",         #  sc  phenotype  sex
                 "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",        #  sc  phenotype  sex
                 "Winnier, DA. (2015)"="GSE64567eset",                                #  sc  phenotype  sex
                 "Nono Nankam, PA. (2020)"="Nankam_et_al",                            #  sc  phenotype
                 "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",                 #  sc  phenotype  sex
                 "Barberio, MD. (2019)"="GSE88837eset",                               #  om  phenotype       ethnicity
                 "Sharma, NK. (2016)"="GSE95674eset",                                 #  sc             sex
                 "GTEx microarray"="GTEx_microarray_sc",                              #  sc             sex
                 "Lonsdale, J. (2013) om"="GTEx_v4_om",                               #  om             sex
                 "Lonsdale, J. (2013) sc"="GTEx_v4_sc",                               #  sc             sex
                 "Aguet, F. (2017) om"="GTEx_v6_om",                                  #  om             sex
                 "Aguet, F. (2017) sc"="GTEx_v6_sc",                                  #  sc             sex
                 "GTEx v7 om"="GTEx_v7_om",                                           #  om             sex
                 "GTEx v7 sc"="GTEx_v7_sc",                                           #  sc             sex
                 "Aguet, F. (2020) om"="GTEx_v8_om",                                  #  om             sex
                 "Aguet, F. (2020) sc"="GTEx_v8_sc",                                  #  sc             sex
                 "Drong, A. (2013)"="E_MTAB_54",                                      #  sc             sex
                 "Das, SK. (2015)"="GSE65221eset",                                    #  sc             sex  ethnicity
                 "Naukkarinen, J. (2013)"="E_MTAB_1895",                              #  sc             sex
                 "Nookaew, I. (2013)"="GSE27916eset",                                 #  sc             sex
                 "Kerr, A. (2020) weight loss"="DEOSHeset_WL",                        #  wl
                 "Petrus, P. (2018) weight loss"="POeset_WL",                         #  wl
                 "Imbert, A. (2022) weight loss"="GSE141221_diogenes1_WL",            #  wl
                 "Armenise, C. (2017) weight loss"="GSE95640_diogenes2_WL",           #  wl
                 "Vink, RG. (2017) weight loss"="GSE77962eset_WL"                     #  wl
)




citation <- readRDS("./data/citation.RDS")

## Clinical association ----


cor_list_DEOSHeset_Baseline_pearson <- readRDS("./data/cor_list_DEOSHeset_Baseline_pearson.RDS")
cor_list_EMIFeset_om_pearson <- readRDS("./data/cor_list_EMIFeset_om_pearson.RDS")
cor_list_EMIFeset_sc_pearson <- readRDS("./data/cor_list_EMIFeset_sc_pearson.RDS")
cor_list_GSE141221_diogenes1_Baseline_pearson <- readRDS("./data/cor_list_GSE141221_diogenes1_Baseline_pearson.RDS")
cor_list_GSE64567eset_pearson <- readRDS("./data/cor_list_GSE64567eset_pearson.RDS")
cor_list_GSE77962eset_Baseline_pearson <- readRDS("./data/cor_list_GSE77962eset_Baseline_pearson.RDS")
cor_list_GSE88837eset_pearson <- readRDS("./data/cor_list_GSE88837eset_pearson.RDS")
cor_list_GSE95640_diogenes2_Baseline_pearson <- readRDS("./data/cor_list_GSE95640_diogenes2_Baseline_pearson.RDS")
cor_list_Keller_et_al_om_pearson <- readRDS("./data/cor_list_Keller_et_al_om_pearson.RDS")
cor_list_Keller_et_al_sc_pearson <- readRDS("./data/cor_list_Keller_et_al_sc_pearson.RDS")
cor_list_Krieg_et_al_om_pearson <- readRDS("./data/cor_list_Krieg_et_al_om_pearson.RDS")
cor_list_Krieg_et_al_sc_pearson <- readRDS("./data/cor_list_Krieg_et_al_sc_pearson.RDS")
cor_list_METSIM_204eset_pearson <- readRDS("./data/cor_list_METSIM_204eset_pearson.RDS")
cor_list_METSIM_434eset_pearson <- readRDS("./data/cor_list_METSIM_434eset_pearson.RDS")
cor_list_METSIM_770eset_pearson <- readRDS("./data/cor_list_METSIM_770eset_pearson.RDS")
cor_list_Nankam_et_al_pearson <- readRDS("./data/cor_list_Nankam_et_al_pearson.RDS")
cor_list_POeset_Baseline_pearson <- readRDS("./data/cor_list_POeset_Baseline_pearson.RDS")
cor_list_RIKENeset_pearson <- readRDS("./data/cor_list_RIKENeset_pearson.RDS")
cor_list_SOWOTeset_pearson <- readRDS("./data/cor_list_SOWOTeset_pearson.RDS")
cor_list_META_sc_pearson <- readRDS("./data/cor_list_META_sc_pearson.RDS")
cor_list_META_om_pearson <- readRDS("./data/cor_list_META_om_pearson.RDS")
cor_list_Proteomicseset_pearson <- readRDS("./data/cor_list_Proteomicseset_pearson.RDS")


cor_list_DEOSHeset_Baseline_spearman <- readRDS("./data/cor_list_DEOSHeset_Baseline_spearman.RDS")
cor_list_EMIFeset_om_spearman <- readRDS("./data/cor_list_EMIFeset_om_spearman.RDS")
cor_list_EMIFeset_sc_spearman <- readRDS("./data/cor_list_EMIFeset_sc_spearman.RDS")
cor_list_GSE141221_diogenes1_Baseline_spearman <- readRDS("./data/cor_list_GSE141221_diogenes1_Baseline_spearman.RDS")
cor_list_GSE64567eset_spearman <- readRDS("./data/cor_list_GSE64567eset_spearman.RDS")
cor_list_GSE77962eset_Baseline_spearman <- readRDS("./data/cor_list_GSE77962eset_Baseline_spearman.RDS")
cor_list_GSE88837eset_spearman <- readRDS("./data/cor_list_GSE88837eset_spearman.RDS")
cor_list_GSE95640_diogenes2_Baseline_spearman <- readRDS("./data/cor_list_GSE95640_diogenes2_Baseline_spearman.RDS")
cor_list_Keller_et_al_om_spearman <- readRDS("./data/cor_list_Keller_et_al_om_spearman.RDS")
cor_list_Keller_et_al_sc_spearman <- readRDS("./data/cor_list_Keller_et_al_sc_spearman.RDS")
cor_list_Krieg_et_al_om_spearman <- readRDS("./data/cor_list_Krieg_et_al_om_spearman.RDS")
cor_list_Krieg_et_al_sc_spearman <- readRDS("./data/cor_list_Krieg_et_al_sc_spearman.RDS")
cor_list_METSIM_204eset_spearman <- readRDS("./data/cor_list_METSIM_204eset_spearman.RDS")
cor_list_METSIM_434eset_spearman <- readRDS("./data/cor_list_METSIM_434eset_spearman.RDS")
cor_list_METSIM_770eset_spearman <- readRDS("./data/cor_list_METSIM_770eset_spearman.RDS")
cor_list_Nankam_et_al_spearman <- readRDS("./data/cor_list_Nankam_et_al_spearman.RDS")
cor_list_POeset_Baseline_spearman <- readRDS("./data/cor_list_POeset_Baseline_spearman.RDS")
cor_list_RIKENeset_spearman <- readRDS("./data/cor_list_RIKENeset_spearman.RDS")
cor_list_SOWOTeset_spearman <- readRDS("./data/cor_list_SOWOTeset_spearman.RDS")
cor_list_META_sc_spearman <- readRDS("./data/cor_list_META_sc_spearman.RDS")
cor_list_META_om_spearman <- readRDS("./data/cor_list_META_om_spearman.RDS")
cor_list_Proteomicseset_spearman <- readRDS("./data/cor_list_Proteomicseset_spearman.RDS")





cor_list_DEOSHeset_Baseline_pearson$r <- cor_list_DEOSHeset_Baseline_pearson$r[,colSums(is.na(cor_list_DEOSHeset_Baseline_pearson$r))<nrow(cor_list_DEOSHeset_Baseline_pearson$r)]
cor_list_DEOSHeset_Baseline_pearson$p <- cor_list_DEOSHeset_Baseline_pearson$p[,colSums(is.na(cor_list_DEOSHeset_Baseline_pearson$p))<nrow(cor_list_DEOSHeset_Baseline_pearson$p)]
cor_list_DEOSHeset_Baseline_pearson$n <- cor_list_DEOSHeset_Baseline_pearson$n[,colSums(is.na(cor_list_DEOSHeset_Baseline_pearson$n))<nrow(cor_list_DEOSHeset_Baseline_pearson$n)]
cor_list_DEOSHeset_Baseline_spearman$r <- cor_list_DEOSHeset_Baseline_spearman$r[,colSums(is.na(cor_list_DEOSHeset_Baseline_spearman$r))<nrow(cor_list_DEOSHeset_Baseline_spearman$r)]
cor_list_DEOSHeset_Baseline_spearman$p <- cor_list_DEOSHeset_Baseline_spearman$p[,colSums(is.na(cor_list_DEOSHeset_Baseline_spearman$p))<nrow(cor_list_DEOSHeset_Baseline_spearman$p)]
cor_list_DEOSHeset_Baseline_spearman$n <- cor_list_DEOSHeset_Baseline_spearman$n[,colSums(is.na(cor_list_DEOSHeset_Baseline_spearman$n))<nrow(cor_list_DEOSHeset_Baseline_spearman$n)]
cor_list_EMIFeset_om_pearson$r <- cor_list_EMIFeset_om_pearson$r[,colSums(is.na(cor_list_EMIFeset_om_pearson$r))<nrow(cor_list_EMIFeset_om_pearson$r)]
cor_list_EMIFeset_om_pearson$p <- cor_list_EMIFeset_om_pearson$p[,colSums(is.na(cor_list_EMIFeset_om_pearson$p))<nrow(cor_list_EMIFeset_om_pearson$p)]
cor_list_EMIFeset_om_pearson$n <- cor_list_EMIFeset_om_pearson$n[,colSums(is.na(cor_list_EMIFeset_om_pearson$n))<nrow(cor_list_EMIFeset_om_pearson$n)]
cor_list_EMIFeset_om_spearman$r <- cor_list_EMIFeset_om_spearman$r[,colSums(is.na(cor_list_EMIFeset_om_spearman$r))<nrow(cor_list_EMIFeset_om_spearman$r)]
cor_list_EMIFeset_om_spearman$p <- cor_list_EMIFeset_om_spearman$p[,colSums(is.na(cor_list_EMIFeset_om_spearman$p))<nrow(cor_list_EMIFeset_om_spearman$p)]
cor_list_EMIFeset_om_spearman$n <- cor_list_EMIFeset_om_spearman$n[,colSums(is.na(cor_list_EMIFeset_om_spearman$n))<nrow(cor_list_EMIFeset_om_spearman$n)]
cor_list_EMIFeset_sc_pearson$r <- cor_list_EMIFeset_sc_pearson$r[,colSums(is.na(cor_list_EMIFeset_sc_pearson$r))<nrow(cor_list_EMIFeset_sc_pearson$r)]
cor_list_EMIFeset_sc_pearson$p <- cor_list_EMIFeset_sc_pearson$p[,colSums(is.na(cor_list_EMIFeset_sc_pearson$p))<nrow(cor_list_EMIFeset_sc_pearson$p)]
cor_list_EMIFeset_sc_pearson$n <- cor_list_EMIFeset_sc_pearson$n[,colSums(is.na(cor_list_EMIFeset_sc_pearson$n))<nrow(cor_list_EMIFeset_sc_pearson$n)]
cor_list_EMIFeset_sc_spearman$r <- cor_list_EMIFeset_sc_spearman$r[,colSums(is.na(cor_list_EMIFeset_sc_spearman$r))<nrow(cor_list_EMIFeset_sc_spearman$r)]
cor_list_EMIFeset_sc_spearman$p <- cor_list_EMIFeset_sc_spearman$p[,colSums(is.na(cor_list_EMIFeset_sc_spearman$p))<nrow(cor_list_EMIFeset_sc_spearman$p)]
cor_list_EMIFeset_sc_spearman$n <- cor_list_EMIFeset_sc_spearman$n[,colSums(is.na(cor_list_EMIFeset_sc_spearman$n))<nrow(cor_list_EMIFeset_sc_spearman$n)]
cor_list_GSE141221_diogenes1_Baseline_pearson$r <- cor_list_GSE141221_diogenes1_Baseline_pearson$r[,colSums(is.na(cor_list_GSE141221_diogenes1_Baseline_pearson$r))<nrow(cor_list_GSE141221_diogenes1_Baseline_pearson$r)]
cor_list_GSE141221_diogenes1_Baseline_pearson$p <- cor_list_GSE141221_diogenes1_Baseline_pearson$p[,colSums(is.na(cor_list_GSE141221_diogenes1_Baseline_pearson$p))<nrow(cor_list_GSE141221_diogenes1_Baseline_pearson$p)]
cor_list_GSE141221_diogenes1_Baseline_pearson$n <- cor_list_GSE141221_diogenes1_Baseline_pearson$n[,colSums(is.na(cor_list_GSE141221_diogenes1_Baseline_pearson$n))<nrow(cor_list_GSE141221_diogenes1_Baseline_pearson$n)]
cor_list_GSE141221_diogenes1_Baseline_spearman$r <- cor_list_GSE141221_diogenes1_Baseline_spearman$r[,colSums(is.na(cor_list_GSE141221_diogenes1_Baseline_spearman$r))<nrow(cor_list_GSE141221_diogenes1_Baseline_spearman$r)]
cor_list_GSE141221_diogenes1_Baseline_spearman$p <- cor_list_GSE141221_diogenes1_Baseline_spearman$p[,colSums(is.na(cor_list_GSE141221_diogenes1_Baseline_spearman$p))<nrow(cor_list_GSE141221_diogenes1_Baseline_spearman$p)]
cor_list_GSE141221_diogenes1_Baseline_spearman$n <- cor_list_GSE141221_diogenes1_Baseline_spearman$n[,colSums(is.na(cor_list_GSE141221_diogenes1_Baseline_spearman$n))<nrow(cor_list_GSE141221_diogenes1_Baseline_spearman$n)]
cor_list_GSE64567eset_pearson$r <- cor_list_GSE64567eset_pearson$r[,colSums(is.na(cor_list_GSE64567eset_pearson$r))<nrow(cor_list_GSE64567eset_pearson$r)]
cor_list_GSE64567eset_pearson$p <- cor_list_GSE64567eset_pearson$p[,colSums(is.na(cor_list_GSE64567eset_pearson$p))<nrow(cor_list_GSE64567eset_pearson$p)]
cor_list_GSE64567eset_pearson$n <- cor_list_GSE64567eset_pearson$n[,colSums(is.na(cor_list_GSE64567eset_pearson$n))<nrow(cor_list_GSE64567eset_pearson$n)]
cor_list_GSE64567eset_spearman$r <- cor_list_GSE64567eset_spearman$r[,colSums(is.na(cor_list_GSE64567eset_spearman$r))<nrow(cor_list_GSE64567eset_spearman$r)]
cor_list_GSE64567eset_spearman$p <- cor_list_GSE64567eset_spearman$p[,colSums(is.na(cor_list_GSE64567eset_spearman$p))<nrow(cor_list_GSE64567eset_spearman$p)]
cor_list_GSE64567eset_spearman$n <- cor_list_GSE64567eset_spearman$n[,colSums(is.na(cor_list_GSE64567eset_spearman$n))<nrow(cor_list_GSE64567eset_spearman$n)]
cor_list_GSE77962eset_Baseline_pearson$r <- cor_list_GSE77962eset_Baseline_pearson$r[,colSums(is.na(cor_list_GSE77962eset_Baseline_pearson$r))<nrow(cor_list_GSE77962eset_Baseline_pearson$r)]
cor_list_GSE77962eset_Baseline_pearson$p <- cor_list_GSE77962eset_Baseline_pearson$p[,colSums(is.na(cor_list_GSE77962eset_Baseline_pearson$p))<nrow(cor_list_GSE77962eset_Baseline_pearson$p)]
cor_list_GSE77962eset_Baseline_pearson$n <- cor_list_GSE77962eset_Baseline_pearson$n[,colSums(is.na(cor_list_GSE77962eset_Baseline_pearson$n))<nrow(cor_list_GSE77962eset_Baseline_pearson$n)]
cor_list_GSE77962eset_Baseline_spearman$r <- cor_list_GSE77962eset_Baseline_spearman$r[,colSums(is.na(cor_list_GSE77962eset_Baseline_spearman$r))<nrow(cor_list_GSE77962eset_Baseline_spearman$r)]
cor_list_GSE77962eset_Baseline_spearman$p <- cor_list_GSE77962eset_Baseline_spearman$p[,colSums(is.na(cor_list_GSE77962eset_Baseline_spearman$p))<nrow(cor_list_GSE77962eset_Baseline_spearman$p)]
cor_list_GSE77962eset_Baseline_spearman$n <- cor_list_GSE77962eset_Baseline_spearman$n[,colSums(is.na(cor_list_GSE77962eset_Baseline_spearman$n))<nrow(cor_list_GSE77962eset_Baseline_spearman$n)]
cor_list_GSE88837eset_pearson$r <- cor_list_GSE88837eset_pearson$r[,colSums(is.na(cor_list_GSE88837eset_pearson$r))<nrow(cor_list_GSE88837eset_pearson$r)]
cor_list_GSE88837eset_pearson$p <- cor_list_GSE88837eset_pearson$p[,colSums(is.na(cor_list_GSE88837eset_pearson$p))<nrow(cor_list_GSE88837eset_pearson$p)]
cor_list_GSE88837eset_pearson$n <- cor_list_GSE88837eset_pearson$n[,colSums(is.na(cor_list_GSE88837eset_pearson$n))<nrow(cor_list_GSE88837eset_pearson$n)]
cor_list_GSE88837eset_spearman$r <- cor_list_GSE88837eset_spearman$r[,colSums(is.na(cor_list_GSE88837eset_spearman$r))<nrow(cor_list_GSE88837eset_spearman$r)]
cor_list_GSE88837eset_spearman$p <- cor_list_GSE88837eset_spearman$p[,colSums(is.na(cor_list_GSE88837eset_spearman$p))<nrow(cor_list_GSE88837eset_spearman$p)]
cor_list_GSE88837eset_spearman$n <- cor_list_GSE88837eset_spearman$n[,colSums(is.na(cor_list_GSE88837eset_spearman$n))<nrow(cor_list_GSE88837eset_spearman$n)]
cor_list_GSE95640_diogenes2_Baseline_pearson$r <- cor_list_GSE95640_diogenes2_Baseline_pearson$r[,colSums(is.na(cor_list_GSE95640_diogenes2_Baseline_pearson$r))<nrow(cor_list_GSE95640_diogenes2_Baseline_pearson$r)]
cor_list_GSE95640_diogenes2_Baseline_pearson$p <- cor_list_GSE95640_diogenes2_Baseline_pearson$p[,colSums(is.na(cor_list_GSE95640_diogenes2_Baseline_pearson$p))<nrow(cor_list_GSE95640_diogenes2_Baseline_pearson$p)]
cor_list_GSE95640_diogenes2_Baseline_pearson$n <- cor_list_GSE95640_diogenes2_Baseline_pearson$n[,colSums(is.na(cor_list_GSE95640_diogenes2_Baseline_pearson$n))<nrow(cor_list_GSE95640_diogenes2_Baseline_pearson$n)]
cor_list_GSE95640_diogenes2_Baseline_spearman$r <- cor_list_GSE95640_diogenes2_Baseline_spearman$r[,colSums(is.na(cor_list_GSE95640_diogenes2_Baseline_spearman$r))<nrow(cor_list_GSE95640_diogenes2_Baseline_spearman$r)]
cor_list_GSE95640_diogenes2_Baseline_spearman$p <- cor_list_GSE95640_diogenes2_Baseline_spearman$p[,colSums(is.na(cor_list_GSE95640_diogenes2_Baseline_spearman$p))<nrow(cor_list_GSE95640_diogenes2_Baseline_spearman$p)]
cor_list_GSE95640_diogenes2_Baseline_spearman$n <- cor_list_GSE95640_diogenes2_Baseline_spearman$n[,colSums(is.na(cor_list_GSE95640_diogenes2_Baseline_spearman$n))<nrow(cor_list_GSE95640_diogenes2_Baseline_spearman$n)]
cor_list_Keller_et_al_om_pearson$r <- cor_list_Keller_et_al_om_pearson$r[,colSums(is.na(cor_list_Keller_et_al_om_pearson$r))<nrow(cor_list_Keller_et_al_om_pearson$r)]
cor_list_Keller_et_al_om_pearson$p <- cor_list_Keller_et_al_om_pearson$p[,colSums(is.na(cor_list_Keller_et_al_om_pearson$p))<nrow(cor_list_Keller_et_al_om_pearson$p)]
cor_list_Keller_et_al_om_pearson$n <- cor_list_Keller_et_al_om_pearson$n[,colSums(is.na(cor_list_Keller_et_al_om_pearson$n))<nrow(cor_list_Keller_et_al_om_pearson$n)]
cor_list_Keller_et_al_om_spearman$r <- cor_list_Keller_et_al_om_spearman$r[,colSums(is.na(cor_list_Keller_et_al_om_spearman$r))<nrow(cor_list_Keller_et_al_om_spearman$r)]
cor_list_Keller_et_al_om_spearman$p <- cor_list_Keller_et_al_om_spearman$p[,colSums(is.na(cor_list_Keller_et_al_om_spearman$p))<nrow(cor_list_Keller_et_al_om_spearman$p)]
cor_list_Keller_et_al_om_spearman$n <- cor_list_Keller_et_al_om_spearman$n[,colSums(is.na(cor_list_Keller_et_al_om_spearman$n))<nrow(cor_list_Keller_et_al_om_spearman$n)]
cor_list_Keller_et_al_sc_pearson$r <- cor_list_Keller_et_al_sc_pearson$r[,colSums(is.na(cor_list_Keller_et_al_sc_pearson$r))<nrow(cor_list_Keller_et_al_sc_pearson$r)]
cor_list_Keller_et_al_sc_pearson$p <- cor_list_Keller_et_al_sc_pearson$p[,colSums(is.na(cor_list_Keller_et_al_sc_pearson$p))<nrow(cor_list_Keller_et_al_sc_pearson$p)]
cor_list_Keller_et_al_sc_pearson$n <- cor_list_Keller_et_al_sc_pearson$n[,colSums(is.na(cor_list_Keller_et_al_sc_pearson$n))<nrow(cor_list_Keller_et_al_sc_pearson$n)]
cor_list_Keller_et_al_sc_spearman$r <- cor_list_Keller_et_al_sc_spearman$r[,colSums(is.na(cor_list_Keller_et_al_sc_spearman$r))<nrow(cor_list_Keller_et_al_sc_spearman$r)]
cor_list_Keller_et_al_sc_spearman$p <- cor_list_Keller_et_al_sc_spearman$p[,colSums(is.na(cor_list_Keller_et_al_sc_spearman$p))<nrow(cor_list_Keller_et_al_sc_spearman$p)]
cor_list_Keller_et_al_sc_spearman$n <- cor_list_Keller_et_al_sc_spearman$n[,colSums(is.na(cor_list_Keller_et_al_sc_spearman$n))<nrow(cor_list_Keller_et_al_sc_spearman$n)]
cor_list_Krieg_et_al_om_pearson$r <- cor_list_Krieg_et_al_om_pearson$r[,colSums(is.na(cor_list_Krieg_et_al_om_pearson$r))<nrow(cor_list_Krieg_et_al_om_pearson$r)]
cor_list_Krieg_et_al_om_pearson$p <- cor_list_Krieg_et_al_om_pearson$p[,colSums(is.na(cor_list_Krieg_et_al_om_pearson$p))<nrow(cor_list_Krieg_et_al_om_pearson$p)]
cor_list_Krieg_et_al_om_pearson$n <- cor_list_Krieg_et_al_om_pearson$n[,colSums(is.na(cor_list_Krieg_et_al_om_pearson$n))<nrow(cor_list_Krieg_et_al_om_pearson$n)]
cor_list_Krieg_et_al_om_spearman$r <- cor_list_Krieg_et_al_om_spearman$r[,colSums(is.na(cor_list_Krieg_et_al_om_spearman$r))<nrow(cor_list_Krieg_et_al_om_spearman$r)]
cor_list_Krieg_et_al_om_spearman$p <- cor_list_Krieg_et_al_om_spearman$p[,colSums(is.na(cor_list_Krieg_et_al_om_spearman$p))<nrow(cor_list_Krieg_et_al_om_spearman$p)]
cor_list_Krieg_et_al_om_spearman$n <- cor_list_Krieg_et_al_om_spearman$n[,colSums(is.na(cor_list_Krieg_et_al_om_spearman$n))<nrow(cor_list_Krieg_et_al_om_spearman$n)]
cor_list_Krieg_et_al_sc_pearson$r <- cor_list_Krieg_et_al_sc_pearson$r[,colSums(is.na(cor_list_Krieg_et_al_sc_pearson$r))<nrow(cor_list_Krieg_et_al_sc_pearson$r)]
cor_list_Krieg_et_al_sc_pearson$p <- cor_list_Krieg_et_al_sc_pearson$p[,colSums(is.na(cor_list_Krieg_et_al_sc_pearson$p))<nrow(cor_list_Krieg_et_al_sc_pearson$p)]
cor_list_Krieg_et_al_sc_pearson$n <- cor_list_Krieg_et_al_sc_pearson$n[,colSums(is.na(cor_list_Krieg_et_al_sc_pearson$n))<nrow(cor_list_Krieg_et_al_sc_pearson$n)]
cor_list_Krieg_et_al_sc_spearman$r <- cor_list_Krieg_et_al_sc_spearman$r[,colSums(is.na(cor_list_Krieg_et_al_sc_spearman$r))<nrow(cor_list_Krieg_et_al_sc_spearman$r)]
cor_list_Krieg_et_al_sc_spearman$p <- cor_list_Krieg_et_al_sc_spearman$p[,colSums(is.na(cor_list_Krieg_et_al_sc_spearman$p))<nrow(cor_list_Krieg_et_al_sc_spearman$p)]
cor_list_Krieg_et_al_sc_spearman$n <- cor_list_Krieg_et_al_sc_spearman$n[,colSums(is.na(cor_list_Krieg_et_al_sc_spearman$n))<nrow(cor_list_Krieg_et_al_sc_spearman$n)]
cor_list_METSIM_204eset_pearson$r <- cor_list_METSIM_204eset_pearson$r[,colSums(is.na(cor_list_METSIM_204eset_pearson$r))<nrow(cor_list_METSIM_204eset_pearson$r)]
cor_list_METSIM_204eset_pearson$p <- cor_list_METSIM_204eset_pearson$p[,colSums(is.na(cor_list_METSIM_204eset_pearson$p))<nrow(cor_list_METSIM_204eset_pearson$p)]
cor_list_METSIM_204eset_pearson$n <- cor_list_METSIM_204eset_pearson$n[,colSums(is.na(cor_list_METSIM_204eset_pearson$n))<nrow(cor_list_METSIM_204eset_pearson$n)]
cor_list_METSIM_204eset_spearman$r <- cor_list_METSIM_204eset_spearman$r[,colSums(is.na(cor_list_METSIM_204eset_spearman$r))<nrow(cor_list_METSIM_204eset_spearman$r)]
cor_list_METSIM_204eset_spearman$p <- cor_list_METSIM_204eset_spearman$p[,colSums(is.na(cor_list_METSIM_204eset_spearman$p))<nrow(cor_list_METSIM_204eset_spearman$p)]
cor_list_METSIM_204eset_spearman$n <- cor_list_METSIM_204eset_spearman$n[,colSums(is.na(cor_list_METSIM_204eset_spearman$n))<nrow(cor_list_METSIM_204eset_spearman$n)]
cor_list_METSIM_434eset_pearson$r <- cor_list_METSIM_434eset_pearson$r[,colSums(is.na(cor_list_METSIM_434eset_pearson$r))<nrow(cor_list_METSIM_434eset_pearson$r)]
cor_list_METSIM_434eset_pearson$p <- cor_list_METSIM_434eset_pearson$p[,colSums(is.na(cor_list_METSIM_434eset_pearson$p))<nrow(cor_list_METSIM_434eset_pearson$p)]
cor_list_METSIM_434eset_pearson$n <- cor_list_METSIM_434eset_pearson$n[,colSums(is.na(cor_list_METSIM_434eset_pearson$n))<nrow(cor_list_METSIM_434eset_pearson$n)]
cor_list_METSIM_434eset_spearman$r <- cor_list_METSIM_434eset_spearman$r[,colSums(is.na(cor_list_METSIM_434eset_spearman$r))<nrow(cor_list_METSIM_434eset_spearman$r)]
cor_list_METSIM_434eset_spearman$p <- cor_list_METSIM_434eset_spearman$p[,colSums(is.na(cor_list_METSIM_434eset_spearman$p))<nrow(cor_list_METSIM_434eset_spearman$p)]
cor_list_METSIM_434eset_spearman$n <- cor_list_METSIM_434eset_spearman$n[,colSums(is.na(cor_list_METSIM_434eset_spearman$n))<nrow(cor_list_METSIM_434eset_spearman$n)]
cor_list_METSIM_770eset_pearson$r <- cor_list_METSIM_770eset_pearson$r[,colSums(is.na(cor_list_METSIM_770eset_pearson$r))<nrow(cor_list_METSIM_770eset_pearson$r)]
cor_list_METSIM_770eset_pearson$p <- cor_list_METSIM_770eset_pearson$p[,colSums(is.na(cor_list_METSIM_770eset_pearson$p))<nrow(cor_list_METSIM_770eset_pearson$p)]
cor_list_METSIM_770eset_pearson$n <- cor_list_METSIM_770eset_pearson$n[,colSums(is.na(cor_list_METSIM_770eset_pearson$n))<nrow(cor_list_METSIM_770eset_pearson$n)]
cor_list_METSIM_770eset_spearman$r <- cor_list_METSIM_770eset_spearman$r[,colSums(is.na(cor_list_METSIM_770eset_spearman$r))<nrow(cor_list_METSIM_770eset_spearman$r)]
cor_list_METSIM_770eset_spearman$p <- cor_list_METSIM_770eset_spearman$p[,colSums(is.na(cor_list_METSIM_770eset_spearman$p))<nrow(cor_list_METSIM_770eset_spearman$p)]
cor_list_METSIM_770eset_spearman$n <- cor_list_METSIM_770eset_spearman$n[,colSums(is.na(cor_list_METSIM_770eset_spearman$n))<nrow(cor_list_METSIM_770eset_spearman$n)]
cor_list_Nankam_et_al_pearson$r <- cor_list_Nankam_et_al_pearson$r[,colSums(is.na(cor_list_Nankam_et_al_pearson$r))<nrow(cor_list_Nankam_et_al_pearson$r)]
cor_list_Nankam_et_al_pearson$p <- cor_list_Nankam_et_al_pearson$p[,colSums(is.na(cor_list_Nankam_et_al_pearson$p))<nrow(cor_list_Nankam_et_al_pearson$p)]
cor_list_Nankam_et_al_pearson$n <- cor_list_Nankam_et_al_pearson$n[,colSums(is.na(cor_list_Nankam_et_al_pearson$n))<nrow(cor_list_Nankam_et_al_pearson$n)]
cor_list_Nankam_et_al_spearman$r <- cor_list_Nankam_et_al_spearman$r[,colSums(is.na(cor_list_Nankam_et_al_spearman$r))<nrow(cor_list_Nankam_et_al_spearman$r)]
cor_list_Nankam_et_al_spearman$p <- cor_list_Nankam_et_al_spearman$p[,colSums(is.na(cor_list_Nankam_et_al_spearman$p))<nrow(cor_list_Nankam_et_al_spearman$p)]
cor_list_Nankam_et_al_spearman$n <- cor_list_Nankam_et_al_spearman$n[,colSums(is.na(cor_list_Nankam_et_al_spearman$n))<nrow(cor_list_Nankam_et_al_spearman$n)]
cor_list_POeset_Baseline_pearson$r <- cor_list_POeset_Baseline_pearson$r[,colSums(is.na(cor_list_POeset_Baseline_pearson$r))<nrow(cor_list_POeset_Baseline_pearson$r)]
cor_list_POeset_Baseline_pearson$p <- cor_list_POeset_Baseline_pearson$p[,colSums(is.na(cor_list_POeset_Baseline_pearson$p))<nrow(cor_list_POeset_Baseline_pearson$p)]
cor_list_POeset_Baseline_pearson$n <- cor_list_POeset_Baseline_pearson$n[,colSums(is.na(cor_list_POeset_Baseline_pearson$n))<nrow(cor_list_POeset_Baseline_pearson$n)]
cor_list_POeset_Baseline_spearman$r <- cor_list_POeset_Baseline_spearman$r[,colSums(is.na(cor_list_POeset_Baseline_spearman$r))<nrow(cor_list_POeset_Baseline_spearman$r)]
cor_list_POeset_Baseline_spearman$p <- cor_list_POeset_Baseline_spearman$p[,colSums(is.na(cor_list_POeset_Baseline_spearman$p))<nrow(cor_list_POeset_Baseline_spearman$p)]
cor_list_POeset_Baseline_spearman$n <- cor_list_POeset_Baseline_spearman$n[,colSums(is.na(cor_list_POeset_Baseline_spearman$n))<nrow(cor_list_POeset_Baseline_spearman$n)]
cor_list_RIKENeset_pearson$r <- cor_list_RIKENeset_pearson$r[,colSums(is.na(cor_list_RIKENeset_pearson$r))<nrow(cor_list_RIKENeset_pearson$r)]
cor_list_RIKENeset_pearson$p <- cor_list_RIKENeset_pearson$p[,colSums(is.na(cor_list_RIKENeset_pearson$p))<nrow(cor_list_RIKENeset_pearson$p)]
cor_list_RIKENeset_pearson$n <- cor_list_RIKENeset_pearson$n[,colSums(is.na(cor_list_RIKENeset_pearson$n))<nrow(cor_list_RIKENeset_pearson$n)]
cor_list_RIKENeset_spearman$r <- cor_list_RIKENeset_spearman$r[,colSums(is.na(cor_list_RIKENeset_spearman$r))<nrow(cor_list_RIKENeset_spearman$r)]
cor_list_RIKENeset_spearman$p <- cor_list_RIKENeset_spearman$p[,colSums(is.na(cor_list_RIKENeset_spearman$p))<nrow(cor_list_RIKENeset_spearman$p)]
cor_list_RIKENeset_spearman$n <- cor_list_RIKENeset_spearman$n[,colSums(is.na(cor_list_RIKENeset_spearman$n))<nrow(cor_list_RIKENeset_spearman$n)]
cor_list_SOWOTeset_pearson$r <- cor_list_SOWOTeset_pearson$r[,colSums(is.na(cor_list_SOWOTeset_pearson$r))<nrow(cor_list_SOWOTeset_pearson$r)]
cor_list_SOWOTeset_pearson$p <- cor_list_SOWOTeset_pearson$p[,colSums(is.na(cor_list_SOWOTeset_pearson$p))<nrow(cor_list_SOWOTeset_pearson$p)]
cor_list_SOWOTeset_pearson$n <- cor_list_SOWOTeset_pearson$n[,colSums(is.na(cor_list_SOWOTeset_pearson$n))<nrow(cor_list_SOWOTeset_pearson$n)]
cor_list_SOWOTeset_spearman$r <- cor_list_SOWOTeset_spearman$r[,colSums(is.na(cor_list_SOWOTeset_spearman$r))<nrow(cor_list_SOWOTeset_spearman$r)]
cor_list_SOWOTeset_spearman$p <- cor_list_SOWOTeset_spearman$p[,colSums(is.na(cor_list_SOWOTeset_spearman$p))<nrow(cor_list_SOWOTeset_spearman$p)]
cor_list_SOWOTeset_spearman$n <- cor_list_SOWOTeset_spearman$n[,colSums(is.na(cor_list_SOWOTeset_spearman$n))<nrow(cor_list_SOWOTeset_spearman$n)]

cor_list_META_sc_pearson$r <- cor_list_META_sc_pearson$r[,colSums(is.na(cor_list_META_sc_pearson$r))<nrow(cor_list_META_sc_pearson$r)]
cor_list_META_sc_pearson$p <- cor_list_META_sc_pearson$p[,colSums(is.na(cor_list_META_sc_pearson$p))<nrow(cor_list_META_sc_pearson$p)]
cor_list_META_sc_pearson$n <- cor_list_META_sc_pearson$n[,colSums(is.na(cor_list_META_sc_pearson$n))<nrow(cor_list_META_sc_pearson$n)]
cor_list_META_sc_spearman$r <- cor_list_META_sc_spearman$r[,colSums(is.na(cor_list_META_sc_spearman$r))<nrow(cor_list_META_sc_spearman$r)]
cor_list_META_sc_spearman$p <- cor_list_META_sc_spearman$p[,colSums(is.na(cor_list_META_sc_spearman$p))<nrow(cor_list_META_sc_spearman$p)]
cor_list_META_sc_spearman$n <- cor_list_META_sc_spearman$n[,colSums(is.na(cor_list_META_sc_spearman$n))<nrow(cor_list_META_sc_spearman$n)]

cor_list_META_om_pearson$r <- cor_list_META_om_pearson$r[,colSums(is.na(cor_list_META_om_pearson$r))<nrow(cor_list_META_om_pearson$r)]
cor_list_META_om_pearson$p <- cor_list_META_om_pearson$p[,colSums(is.na(cor_list_META_om_pearson$p))<nrow(cor_list_META_om_pearson$p)]
cor_list_META_om_pearson$n <- cor_list_META_om_pearson$n[,colSums(is.na(cor_list_META_om_pearson$n))<nrow(cor_list_META_om_pearson$n)]
cor_list_META_om_spearman$r <- cor_list_META_om_spearman$r[,colSums(is.na(cor_list_META_om_spearman$r))<nrow(cor_list_META_om_spearman$r)]
cor_list_META_om_spearman$p <- cor_list_META_om_spearman$p[,colSums(is.na(cor_list_META_om_spearman$p))<nrow(cor_list_META_om_spearman$p)]
cor_list_META_om_spearman$n <- cor_list_META_om_spearman$n[,colSums(is.na(cor_list_META_om_spearman$n))<nrow(cor_list_META_om_spearman$n)]

cor_list_Proteomicseset_pearson$r <- cor_list_Proteomicseset_pearson$r[,colSums(is.na(cor_list_Proteomicseset_pearson$r))<nrow(cor_list_Proteomicseset_pearson$r)]
cor_list_Proteomicseset_pearson$p <- cor_list_Proteomicseset_pearson$p[,colSums(is.na(cor_list_Proteomicseset_pearson$p))<nrow(cor_list_Proteomicseset_pearson$p)]
cor_list_Proteomicseset_pearson$n <- cor_list_Proteomicseset_pearson$n[,colSums(is.na(cor_list_Proteomicseset_pearson$n))<nrow(cor_list_Proteomicseset_pearson$n)]
cor_list_Proteomicseset_spearman$r <- cor_list_Proteomicseset_spearman$r[,colSums(is.na(cor_list_Proteomicseset_spearman$r))<nrow(cor_list_Proteomicseset_spearman$r)]
cor_list_Proteomicseset_spearman$p <- cor_list_Proteomicseset_spearman$p[,colSums(is.na(cor_list_Proteomicseset_spearman$p))<nrow(cor_list_Proteomicseset_spearman$p)]
cor_list_Proteomicseset_spearman$n <- cor_list_Proteomicseset_spearman$n[,colSums(is.na(cor_list_Proteomicseset_spearman$n))<nrow(cor_list_DEOSHeset_Baseline_spearman$n)]


## Gender ----



Keller_et_al_sc <- readRDS("./data/Keller_et_al_sc.rds")
Krieg_et_al_sc <- readRDS("./data/Krieg_et_al_sc.rds")
GSE141221_diogenes1_Baseline <- readRDS("./data/GSE141221_diogenes1_Baseline.rds")
GSE95640_diogenes2_Baseline <- readRDS("./data/GSE95640_diogenes2_Baseline.rds")
GSE64567eset <- readRDS("./data/GSE64567eset.rds")
GSE77962eset_Baseline <- readRDS("./data/GSE77962eset_Baseline.rds")
GSE95674eset <- readRDS("./data/GSE95674eset.rds")
GTEx_microarray_sc <- readRDS("./data/GTEx_microarray_sc.rds")
GTEx_v4_sc <- readRDS("./data/GTEx_v4_sc.rds")
GTEx_v6_sc <- readRDS("./data/GTEx_v6_sc.rds")
GTEx_v7_sc <- readRDS("./data/GTEx_v7_sc.rds")
GTEx_v8_sc <- readRDS("./data/GTEx_v8_sc.rds")
E_MTAB_54 <- readRDS("./data/E_MTAB_54.rds")
GSE65221eset <- readRDS("./data/GSE65221eset.rds")
E_MTAB_1895 <- readRDS("./data/E_MTAB_1895.rds")
GSE27916eset <- readRDS("./data/GSE27916eset.rds")
Keller_et_al_om <- readRDS("./data/Keller_et_al_om.rds")
Krieg_et_al_om <- readRDS("./data/Krieg_et_al_om.rds")
GTEx_v4_om <- readRDS("./data/GTEx_v4_om.rds")
GTEx_v6_om <- readRDS("./data/GTEx_v6_om.rds")
GTEx_v7_om <- readRDS("./data/GTEx_v7_om.rds")
GTEx_v8_om <- readRDS("./data/GTEx_v8_om.rds")



## BMI cross-sectional ----

cohort_BMI_list <- readRDS("./data/cohort_BMI_list.rds")

DEOSHeset_Baseline <- readRDS("./data/DEOSHeset_Baseline.rds")
POeset_Baseline <- readRDS("./data/POeset_Baseline.rds")
SOWOTeset <- readRDS("./data/SOWOTeset.rds")
# Keller_et_al_om <- readRDS("./data/Keller_et_al_om.rds")
# Keller_et_al_sc <- readRDS("./data/Keller_et_al_sc.rds")
RIKENeset <- readRDS("./data/RIKENeset.rds")
METSIM_204eset <- readRDS("./data/METSIM_204eset.rds")
METSIM_434eset <- readRDS("./data/METSIM_434eset.rds")
METSIM_770eset <- readRDS("./data/METSIM_770eset.rds")
# Krieg_et_al_sc <- readRDS("./data/Krieg_et_al_sc.rds")
# Krieg_et_al_om <- readRDS("./data/Krieg_et_al_om.rds")
EMIFeset_om <- readRDS("./data/EMIFeset_om.rds")
EMIFeset_sc <- readRDS("./data/EMIFeset_sc.rds")
# GSE141221_diogenes1 <- readRDS("./data/GSE141221_diogenes1.rds")
# GSE95640_diogenes2 <- readRDS("./data/GSE95640_diogenes2.rds")
# GSE64567eset <- readRDS("./data/GSE64567eset.rds")
Nankam_et_al <- readRDS('./data/Nankam_et_al.rds')
GSE77962eset_Baseline <- readRDS("./data/GSE77962eset_Baseline.rds")
GSE88837eset <- readRDS("./data/GSE88837eset.rds")



## BMI longitudinal ----

DEOSHeset_WL <- readRDS("./data/DEOSHeset_WL.rds")
POeset_WL <- readRDS("./data/POeset_WL.rds")
GSE77962eset_WL <- readRDS("./data/GSE77962eset_WL.rds")
GSE141221_diogenes1_WL <- readRDS("./data/GSE141221_diogenes1_WL.rds")
GSE95640_diogenes2_WL <- readRDS("./data/GSE95640_diogenes2_WL.rds")


GSE77962eset_WL$time_point <- GSE77962eset_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,5,9))
GSE77962eset_WL <- GSE77962eset_WL[,colnames(GSE77962eset_WL)[GSE77962eset_WL$subject %in% names(table(GSE77962eset_WL$subject)[table(GSE77962eset_WL$subject)>2])]]


GSE141221_diogenes1_WL$time_point <- GSE141221_diogenes1_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,8))
# GSE141221_diogenes1_WL <- GSE141221_diogenes1_WL[,colnames(GSE141221_diogenes1_WL)[GSE141221_diogenes1_WL$subject %in% names(table(GSE141221_diogenes1_WL$subject)[table(GSE141221_diogenes1_WL$subject)>1])]]

GSE95640_diogenes2_WL$time_point <- GSE95640_diogenes2_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,8))
# GSE95640_diogenes2_WL <- GSE95640_diogenes2_WL[,colnames(GSE95640_diogenes2_WL)[GSE95640_diogenes2_WL$subject %in% names(table(GSE95640_diogenes2_WL$subject)[table(GSE95640_diogenes2_WL$subject)>1])]]


pData(DEOSHeset_WL) <- pData(DEOSHeset_WL)[,colSums(is.na(pData(DEOSHeset_WL))) < nrow(pData(DEOSHeset_WL))]
pData(POeset_WL) <- pData(POeset_WL)[,colSums(is.na(pData(POeset_WL))) < nrow(pData(POeset_WL))]
pData(GSE77962eset_WL) <- pData(GSE77962eset_WL)[,colSums(is.na(pData(GSE77962eset_WL))) < nrow(pData(GSE77962eset_WL))]
pData(GSE141221_diogenes1_WL) <- pData(GSE141221_diogenes1_WL)[,colSums(is.na(pData(GSE141221_diogenes1_WL))) < nrow(pData(GSE141221_diogenes1_WL))]
pData(GSE95640_diogenes2_WL) <- pData(GSE95640_diogenes2_WL)[,colSums(is.na(pData(GSE95640_diogenes2_WL))) < nrow(pData(GSE95640_diogenes2_WL))]


trait_vector_bmilong <- trait_vector[trait_vector %in% unique(c(colnames(pData(DEOSHeset_WL)),
                                                                colnames(pData(POeset_WL)),
                                                                colnames(pData(GSE77962eset_WL)),
                                                                colnames(pData(GSE141221_diogenes1_WL)),
                                                                colnames(pData(GSE95640_diogenes2_WL))))]





## RNA-seq log2 + 1

exprs(GTEx_v4_sc) <- log2(exprs(GTEx_v4_sc)+1)
exprs(GTEx_v6_sc) <- log2(exprs(GTEx_v6_sc)+1)
exprs(GTEx_v7_sc) <- log2(exprs(GTEx_v7_sc)+1)
exprs(GTEx_v8_sc) <- log2(exprs(GTEx_v8_sc)+1)

exprs(GTEx_v4_om) <- log2(exprs(GTEx_v4_om)+1)
exprs(GTEx_v6_om) <- log2(exprs(GTEx_v6_om)+1)
exprs(GTEx_v7_om) <- log2(exprs(GTEx_v7_om)+1)
exprs(GTEx_v8_om) <- log2(exprs(GTEx_v8_om)+1)


exprs(GSE141221_diogenes1_Baseline) <- log2(exprs(GSE141221_diogenes1_Baseline)+1)
exprs(GSE95640_diogenes2_Baseline) <- log2(exprs(GSE95640_diogenes2_Baseline)+1)
exprs(METSIM_434eset) <- log2(exprs(METSIM_434eset)+1)

exprs(GSE141221_diogenes1_WL) <- log2(exprs(GSE141221_diogenes1_WL)+1)
exprs(GSE95640_diogenes2_WL) <- log2(exprs(GSE95640_diogenes2_WL)+1)




Proteomicseset <- readRDS("./data/Proteomicseset.rds")


## color ----


PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")
portal_col_ex <- colorRampPalette(portalcol2)