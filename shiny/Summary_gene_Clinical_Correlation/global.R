# Load packages ----

suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(dplyr))




## Clinical ----

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
                  "iso/basal lipolysis" = "iso_basal")

meta_r <- readRDS("./data/meta_sc_r_spearman.RDS")
meta_p <- readRDS("./data/meta_sc_p_spearman.RDS")
meta_n <- readRDS("./data/meta_sc_n_spearman.RDS")

# protein_r <- readRDS("./data/cor_list_Proteomicseset_spearman.RDS")$r
# protein_p <- readRDS("./data/cor_list_Proteomicseset_spearman.RDS")$p
# protein_n <- readRDS("./data/cor_list_Proteomicseset_spearman.RDS")$n





protein_r <- readRDS("./data/meta_sc_r_spearman_proteomics.RDS") %>% as.matrix()
protein_p <- readRDS("./data/meta_sc_p_spearman_proteomics.RDS") %>% as.matrix()
protein_n <- readRDS("./data/meta_sc_n_spearman_proteomics.RDS") %>% as.matrix()



protein_r <- protein_r[,colSums(is.na(protein_r))<nrow(protein_r)]
protein_p <- protein_p[,colSums(is.na(protein_p))<nrow(protein_p)]
protein_n <- protein_n[,colSums(is.na(protein_n))<nrow(protein_n)]
# protein_n <- protein_n[,colSums(protein_n)>0]

# protein_n[!is.na(protein_n)] <- 1

