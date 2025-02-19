# Load packages ----

suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(dplyr))


## Clinical ----


DEOSHeset_WL <- readRDS("./data/DEOSHeset_WL.rds")
POeset_WL <- readRDS("./data/POeset_WL.rds")
# GSE77962eset_WL <- readRDS("./data/GSE77962eset_WL.rds")
GSE141221_diogenes1_WL <- readRDS("./data/GSE141221_diogenes1_WL.rds")
GSE95640_diogenes2_WL <- readRDS("./data/GSE95640_diogenes2_WL.rds")



DEOSHeset_WL <- readRDS("./data/DEOSHeset_WL.rds")
POeset_WL <- readRDS("./data/POeset_WL.rds")
# GSE77962eset_WL <- readRDS("./data/GSE77962eset_WL.rds")
GSE141221_diogenes1_WL <- readRDS("./data/GSE141221_diogenes1_WL.rds")
GSE95640_diogenes2_WL <- readRDS("./data/GSE95640_diogenes2_WL.rds")


# GSE77962eset_WL$time_point <- GSE77962eset_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,5,9))
# GSE77962eset_WL <- GSE77962eset_WL[,colnames(GSE77962eset_WL)[GSE77962eset_WL$subject %in% names(table(GSE77962eset_WL$subject)[table(GSE77962eset_WL$subject)>2])]]


GSE141221_diogenes1_WL$time_point <- GSE141221_diogenes1_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,8))
# GSE141221_diogenes1_WL <- GSE141221_diogenes1_WL[,colnames(GSE141221_diogenes1_WL)[GSE141221_diogenes1_WL$subject %in% names(table(GSE141221_diogenes1_WL$subject)[table(GSE141221_diogenes1_WL$subject)>1])]]

GSE95640_diogenes2_WL$time_point <- GSE95640_diogenes2_WL$time_point %>% as.character() %>% gsub(pattern = " week",replacement = "") %>% as.numeric() %>% factor(levels = c(0,8))
# GSE95640_diogenes2_WL <- GSE95640_diogenes2_WL[,colnames(GSE95640_diogenes2_WL)[GSE95640_diogenes2_WL$subject %in% names(table(GSE95640_diogenes2_WL$subject)[table(GSE95640_diogenes2_WL$subject)>1])]]


#define function to scale values between 0 and 1
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}