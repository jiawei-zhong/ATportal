# Load packages ----

suppressMessages(library(shiny))
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
suppressMessages(library(fresh))
suppressMessages(library(shinycustomloader))



#### Set width and margin in CSS ####


shinyUI(fluidPage(
  tags$head(
    # tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/iframe-resizer/4.3.1/iframeResizer.contentWindow.min.js"),
    tags$style(HTML("
      .dataTables_wrapper .dataTable td {
        white-space: nowrap;
      }
    ")),
     tags$link(rel = "stylesheet", type = "text/css", href = "style_def.css")
  ),
  use_googlefont("Red Hat Display"),
  use_theme(create_theme(
    theme = "default",
    bs_vars_font(
      family_sans_serif = "'Red Hat Display', cursive"
    )
  )),
  fluidRow(
    column(width = 12,
           style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
           column(width = 2,
                  style='padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  selectInput(inputId = "Dataset_STx",
                              label = "Dataset: ",
                              width = "200px",
                              choices = c("Bäckdahl et al. baseline"="Jesper_et_al_baseline", "Bäckdahl et al. insulin"="Jesper_et_al_insulin"),
                              selected = "Jesper_et_al_baseline")),
           column(width = 2,
                  style='padding-left:15px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  selectizeInput(inputId = 'gene_for_featureplot_STx',
                                 label = 'Gene: ',
                                 choices = NULL,
                                 # selected = "LEP",
                                 width = "200px",
                                 # inline = TRUE,
                                 multiple = FALSE)),
           # column(width = 2,
           #        style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
           #        selectizeInput(inputId = 'deconvolution_for_featureplot_STx',
           #                       label = 'Deconvolution: ',
           #                       choices = c("None"="", "LEP"="LEP", "CCND1"="CCND1"),
           #                       selected = "",
           #                       width = "200px",
           #                       # inline = TRUE,
           #                       multiple = FALSE)),
           column(width = 2,
                  style='padding-left:15px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  selectizeInput(inputId = 'deconvolution_for_featureplot_STx',
                                 label = 'Deconvolution: ',
                                 choices = NULL,
                                 # selected = "LEP",
                                 width = "200px",
                                 # inline = TRUE,
                                 multiple = FALSE)),
           
           column(width = 6,
                  style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                  selectInput(inputId = "Slide_STx",
                              label = "Slide: ",
                              width = "200px",
                              choices = NA,
                              selected = NULL))
    ), # column
    
    column(width = 12,
           style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
           column(width = 6,
                  style='padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  tabsetPanel(id = "inTabset",
                              tabPanel(title = "FeaturePlot", 
                                       column(12,uiOutput("ui_FeaturePlot_STx"))),
                              tabPanel(title = "Violin Plot",
                                       column(12,uiOutput("ui_ViolinPlot_STx"))),
                              tabPanel(title = "Marker",
                                       column(12,uiOutput("ui_DT_STx")))
                  )),
           column(width = 6,
                  offset = 0, 
                  style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                  tabsetPanel(id = "inTabset",
                              tabPanel(title = "SpatialDimPlot",
                                       column(12,uiOutput("ui_SpatialDimPlot_STx"))),
                              tabPanel(title = "SpatialFeaturePlot", 
                                       column(12,uiOutput("ui_SpatialFeaturePlot_STx")))
                  ))
    ) # column
    
  ) # fluidRow
))




