# Load packages ----

suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(shiny))
suppressMessages(library(viridis))
suppressMessages(library(dplyr))


## Clinical ----



# adipogenesis <- readRDS("./data/adipogenesis.RDS")
# 
# adipogenesis$log_timepoint <- log1p(adipogenesis$timepoint)
# 
# exprs(adipogenesis) <- scale(exprs(adipogenesis), center = TRUE, scale = TRUE)


# adipogenesis <- read.table("./data/TimeCourse.txt",header = 1)

# adipogenesis$log_timepoint <- log1p(adipogenesis$timepoint)

adipogenesis <- readRDS("./data/adipogenesis.RDS")

