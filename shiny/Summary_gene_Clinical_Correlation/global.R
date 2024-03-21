# Load packages ----

suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))



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
                  "iso/basal" = "iso_basal")

meta_r <- readRDS("./data/meta_sc_r_spearman.RDS")
meta_p <- readRDS("./data/meta_sc_p_spearman.RDS")


