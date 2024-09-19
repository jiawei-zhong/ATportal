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

perturbation <- readRDS("./data/All_Pertubation_ForSummaryVolcano.rds")

perturbation$group <- "Pharmacologic"
perturbation$group[perturbation$perturbation %in% c("Hypoxia", "TNF", "TNFa", "IL6", "TGFB1", "Inflammation")] <- "Cell stress"
perturbation$group[perturbation$perturbation %in% c("Insulin in vivo", "Insulin", "IGF1", "Adiponectin", "Leptin")] <- "Hormonal"
perturbation$group[perturbation$perturbation %in% c("Glucose", "RetinoicAcid", "Lauroylcarnitine", "Decanoyllcarnithine")] <- "Metabolic"
# perturbation$group[perturbation$perturbation %in% c()] <- "Metabolic"

perturbation$log10.adj.P.Val <- (-log10(perturbation$adj.P.Val))




  
# 
# p <- EnhancedVolcano(perturbation,
#                      lab = perturbation$perturbation,
#                      x = 'logFC',
#                      y = 'adj.P.Val',
#                      selectLab = NULL,
#                      colCustom = color_mapping,
#                      xlim = c(min(perturbation$logFC), max(perturbation$logFC)),
#                      ylim = c(0, max(perturbation$log10.adj.P.Val)),
#                      pointSize = 3.0,
#                      colAlpha = 1,
#                      labSize = 6.0,
#                      title = NULL,
#                      subtitle = NULL,
#                      caption = NULL,
#                      pCutoff = NA,
#                      FCcutoff = NA,
#                      legendPosition = 'none',
#                      gridlines.major = FALSE,
#                      gridlines.minor = FALSE) +
#   facet_grid(cols = vars(group)) +
#   geom_text_repel(aes(label = labels), 
#                   box.padding = 0.35, point.padding = 0.3, 
#                   max.overlaps = 10, size = 5, na.rm = TRUE) 
# 
# p <- as.ggplot(p)

