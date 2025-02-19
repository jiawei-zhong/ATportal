# Load packages ----
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(extrafont))
suppressMessages(library(dplyr))



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
                 "Raulerson, C. (2019)"="METSIM_434eset",                             #  sc  phenotype       RNAseq
                 "Civelek, M. (2017)"="METSIM_770eset",                               #  sc  phenotype
                 "Krieg, L. (2021) sc"="Krieg_et_al_sc",                              #  sc  phenotype  sex
                 "Krieg, L. (2021) om"="Krieg_et_al_om",                              #  om  phenotype  sex
                 "Arner, P. (2016) om"="EMIFeset_om",                                 #  om  phenotype
                 "Arner, P. (2016) sc"="EMIFeset_sc",                                 #  sc  phenotype
                 "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",         #  sc  phenotype  sex  RNAseq
                 "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",        #  sc  phenotype  sex  RNAseq
                 "Winnier, DA. (2015)"="GSE64567eset",                                #  sc  phenotype  sex
                 "Nono Nankam, PA. (2020)"="Nankam_et_al",                            #  sc  phenotype
                 "Salcedo-Tacuma, D. (2022) sc"="GSE188799eset_sc",                   #  sc  phenotype       RNAseq
                 "Salcedo-Tacuma, D. (2022) om"="GSE188799eset_om",                   #  om  phenotype       RNAseq
                 "MacLaren, RE. (2010) sc"="GSE15524eset_sc",                         #  sc  phenotype  sex
                 "MacLaren, RE. (2010) om"="GSE15524eset_om",                         #  om  phenotype  sex
                 "Du Plessis, J. (2015) sc"="GSE58979eset_sc",                        #  sc  phenotype  sex
                 "Du Plessis, J. (2015) om"="GSE58979eset_om",                        #  om  phenotype  sex
                 "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",                 #  sc  phenotype  sex
                 "Johansson, LE. (2012) Baseline"="GSE35411eset_Baseline",            #  sc  phenotype  sex
                 "Matualatupauw, JC. (2017)"="GSE87382eset",                          #  sc  phenotype  sex
                 "Barberio, MD. (2019)"="GSE88837eset",                               #  om  phenotype       ethnicity
                 "Aguilera, CM. (2015)"="GSE9624eset",                                #  om  phenotype
                 "Defour, M. (2020)"="GSE154610eset",                                 #  sc  phenotype  sex
                 "Van Bussel, IPG. (2017)"="GSE84046eset_Baseline",                   #  sc  phenotype  sex
                 "Grundberg, E. (2012)"="E_TABM_1140",                                #  sc  phenotype
                 "Heinonen, S. (2017)"="GSE92405eset",                                #  sc  phenotype  sex
                 "Sharma, NK. (2016)"="GSE95674eset",                                 #  sc  phenotype  sex
                 "GTEx microarray"="GTEx_microarray_sc",                              #  sc             sex
                 "Lonsdale, J. (2013) om"="GTEx_v4_om",                               #  om             sex  RNAseq
                 "Lonsdale, J. (2013) sc"="GTEx_v4_sc",                               #  sc             sex  RNAseq
                 "Aguet, F. (2017) om"="GTEx_v6_om",                                  #  om             sex  RNAseq
                 "Aguet, F. (2017) sc"="GTEx_v6_sc",                                  #  sc             sex  RNAseq
                 "GTEx v7 om"="GTEx_v7_om",                                           #  om             sex  RNAseq
                 "GTEx v7 sc"="GTEx_v7_sc",                                           #  sc             sex  RNAseq
                 "Aguet, F. (2020) om"="GTEx_v8_om",                                  #  om             sex  RNAseq
                 "Aguet, F. (2020) sc"="GTEx_v8_sc",                                  #  sc             sex  RNAseq
                 "Drong, A. (2013)"="E_MTAB_54",                                      #  sc             sex
                 "Das, SK. (2015)"="GSE65221eset",                                    #  sc             sex  ethnicity
                 "Naukkarinen, J. (2013)"="E_MTAB_1895",                              #  sc             sex
                 "Nookaew, I. (2013)"="GSE27916eset",                                 #  sc             sex
                 "Hardy, OT. (2011) sc"="GSE20950eset_sc",                            #  sc             sex
                 "Hardy, OT. (2011) om"="GSE20950eset_om",                            #  om             sex
                 "Bollepalli, S. (2018)"="GSE103766eset_Baseline",                    #  sc             sex
                 "Rey, F. (2021)"="GSE166047eset",                                    #  sc             sex
                 "Kerr, A. (2020) weight loss"="DEOSHeset_WL",                        #  wl
                 "Petrus, P. (2018) weight loss"="POeset_WL",                         #  wl                  
                 "Imbert, A. (2022) weight loss"="GSE141221_diogenes1_WL",            #  wl                  RNAseq
                 "Armenise, C. (2017) weight loss"="GSE95640_diogenes2_WL",           #  wl                  RNAseq
                 "Vink, RG. (2017) weight loss"="GSE77962eset_WL",                    #  wl
                 "Johansson, LE. (2012) weight loss"="GSE35411eset_WL"                #  wl                 
)


## load all ExpressionSet ----

for (i in list.files(path = "./data") %>% grep(pattern = ".rds$",value = T)) {
  a <- paste0(gsub(pattern = ".rds",replacement = "",x = i),' <- readRDS("./data/',i,'")')
  eval(parse(text=a))
}



citation <- readRDS("./data/citation.RDS")
meta_sc_combination <- readRDS("./data/meta_sc_combination.RDS")
meta_om_combination <- readRDS("./data/meta_om_combination.RDS")

## Clinical association ----




for (i in list.files(path = "./data") %>% grep(pattern = "^cor_list",value = T)) {
  a <- paste0(gsub(pattern = ".RDS",replacement = "",x = i),' <- readRDS("./data/',i,'")')
  eval(parse(text=a))
  
  temp <- get(gsub(pattern = ".RDS",replacement = "",x = i))
  temp$r <- temp$r[,colSums(is.na(temp$r))<nrow(temp$r),drop=F]
  temp$p <- temp$p[,colSums(is.na(temp$p))<nrow(temp$p),drop=F]
  temp$n <- temp$n[,colSums(is.na(temp$n))<nrow(temp$n),drop=F]
  assign(x = gsub(pattern = ".RDS",replacement = "",x = i),value = temp)
}






## Gender ----





## BMI cross-sectional ----


cohort_BMI_list <- list()
for (i in list.files(path = "./data") %>% grep(pattern = ".rds$",value = T)) {
  # plot_list <- append(x = plot_list,values = list(p))  # 这个也可以
  cohort_BMI_list[[gsub(pattern = ".rds",replacement = "",x = i)]] <- get(gsub(pattern = ".rds",replacement = "",x = i))$BMI_catelogy %>% unique() %>% sort() %>% as.character()
  # cohort_BMI_list[[gsub(pattern = ".rds",replacement = "",x = i)]] <- get(gsub(pattern = ".rds",replacement = "",x = i))$BMI_catelogy %>% table() 
}




## BMI longitudinal ----




GSE141221_diogenes1_WL$time_point <- GSE141221_diogenes1_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,8))
# GSE141221_diogenes1_WL <- GSE141221_diogenes1_WL[,colnames(GSE141221_diogenes1_WL)[GSE141221_diogenes1_WL$subject %in% names(table(GSE141221_diogenes1_WL$subject)[table(GSE141221_diogenes1_WL$subject)>1])]]

GSE95640_diogenes2_WL$time_point <- GSE95640_diogenes2_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,8))
# GSE95640_diogenes2_WL <- GSE95640_diogenes2_WL[,colnames(GSE95640_diogenes2_WL)[GSE95640_diogenes2_WL$subject %in% names(table(GSE95640_diogenes2_WL$subject)[table(GSE95640_diogenes2_WL$subject)>1])]]


GSE77962eset_WL$time_point <- GSE77962eset_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,5,12))
GSE77962eset_WL <- GSE77962eset_WL[,colnames(GSE77962eset_WL)[GSE77962eset_WL$subject %in% names(table(GSE77962eset_WL$subject)[table(GSE77962eset_WL$subject)>1])]]


GSE35411eset_WL$time_point <- GSE35411eset_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,14))
GSE35411eset_WL <- GSE35411eset_WL[,colnames(GSE35411eset_WL)[GSE35411eset_WL$subject %in% names(table(GSE35411eset_WL$subject)[table(GSE35411eset_WL$subject)>1])]]


pData(DEOSHeset_WL) <- pData(DEOSHeset_WL)[,colSums(is.na(pData(DEOSHeset_WL))) < nrow(pData(DEOSHeset_WL))]
pData(POeset_WL) <- pData(POeset_WL)[,colSums(is.na(pData(POeset_WL))) < nrow(pData(POeset_WL))]
pData(GSE141221_diogenes1_WL) <- pData(GSE141221_diogenes1_WL)[,colSums(is.na(pData(GSE141221_diogenes1_WL))) < nrow(pData(GSE141221_diogenes1_WL))]
pData(GSE95640_diogenes2_WL) <- pData(GSE95640_diogenes2_WL)[,colSums(is.na(pData(GSE95640_diogenes2_WL))) < nrow(pData(GSE95640_diogenes2_WL))]
pData(GSE77962eset_WL) <- pData(GSE77962eset_WL)[,colSums(is.na(pData(GSE77962eset_WL))) < nrow(pData(GSE77962eset_WL))]
pData(GSE35411eset_WL) <- pData(GSE35411eset_WL)[,colSums(is.na(pData(GSE35411eset_WL))) < nrow(pData(GSE35411eset_WL))]


trait_vector_bmilong <- trait_vector[trait_vector %in% unique(c(colnames(pData(DEOSHeset_WL)),
                                                                colnames(pData(POeset_WL)),
                                                                colnames(pData(GSE141221_diogenes1_WL)),
                                                                colnames(pData(GSE95640_diogenes2_WL)),
                                                                colnames(pData(GSE77962eset_WL)),
                                                                colnames(pData(GSE35411eset_WL))))]





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

exprs(GSE188799eset_sc) <- log2(exprs(GSE188799eset_sc)+1)
exprs(GSE188799eset_om) <- log2(exprs(GSE188799eset_om)+1)


exprs(GSE166047eset) <- log2(exprs(GSE166047eset)+1)


## color ----


PRGn <- colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
RdBu <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))(51)
virid <- viridis(51)
portalcol <- colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)
portalcol2 <- c("#1a4659", "#e2c744", "#ba5c28", "#4f736a", "#f0cbbf", "#b2dfda", "#f59b7c", "#7dbfb3", "#939a64", "#8ca5a5", "#adbec4", "#c3b6b6")
portal_col_ex <- colorRampPalette(portalcol2)