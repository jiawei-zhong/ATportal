# Load packages ----
suppressMessages(library(shiny))
suppressMessages(library(shinythemes))
suppressMessages(library(shinyhelper))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))  # for box()
suppressMessages(library(shinycustomloader))
suppressMessages(library(DT))
suppressMessages(library(meta))
suppressMessages(library(ppcor))
suppressMessages(library(Hmisc))
suppressMessages(library(psych))
suppressMessages(library(MKmisc))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(grid))
suppressMessages(library(plotly))
suppressMessages(library(pheatmap))
suppressMessages(library(pals))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(cowplot))
suppressMessages(library(ggplotify))
suppressMessages(library(ggpubr))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(scCustomize))

shinyServer(function(input, output, session) {
  
  #### shinyhelper control #### 
  shinyhelper::observe_helpers()
  
  # #### variables ####
  # variables <- reactiveValues(
  #     #### Profiling ####
  #     PRO.exp.data.demo = NULL,
  #     PRO.lipid.char.tab.demo = NULL,
  #     
  #     PRO.exp.data.user = NULL,
  #     PRO.lipid.char.tab.user = NULL,
  #     PRO.sample.count.user = NULL,
  #     PRO.exp.user.col1 = NULL,
  #     PRO.lipid.char.user.col1 = NULL,
  #     
  #     PRO.pca.result = NULL,
  #     PRO.plsda.result = NULL,
  #     PRO.tsne.result = NULL,
  #     PRO.umap.result = NULL,
  #     
  #     PRO.corr.heatmap.result = NULL
  #     
  # 
  # ) #variables <- reactiveValues
  
  
  
  # source(file = "server_Home.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Clinical_Heatmap.R" ,local = TRUE ,encoding = "UTF-8")
  
  source(file = "server_Clinical_Forest.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Clinical_Scatter.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Clinical_Weightloss.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Clinical_BMI.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Clinical_HOMA.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Clinical_Sex.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Depots_Forest.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Depots_Boxplot.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Depots_Heatmap.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Depots_Volcano.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Depots_GO.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Characterization.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Singlecell.R", local = TRUE, encoding = "UTF-8")
  
  source(file = "server_Spatial.R", local = TRUE, encoding = "UTF-8")
  
  
  
  session$allowReconnect(TRUE)
  
})
