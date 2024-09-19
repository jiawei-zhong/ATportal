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
           column(width = 3,
                  style='padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  selectInput(inputId = "Depots",
                              label = "Depots: ",
                              width = "200px",
                              choices = c("All"="META_all", "Subcutaneous WAT"="META_all_sc", "Omental WAT"="META_all_om", "Perivascular WAT"="META_all_pvat"),
                              selected = "META_all")),
           column(width = 3,
                  style='padding-left:15px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  selectizeInput(inputId = 'gene_for_featureplot',
                                 label = 'Gene: ',
                                 choices = NULL,
                                 selected = "LEP",
                                 width = "200px",
                                 # inline = TRUE,
                                 multiple = FALSE)),
           
           column(width = 6,
                  style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                  selectInput(inputId = "Subcluster",
                              label = "Subcluster: ",
                              width = "200px",
                              choices = c("T, NK & NKT"="Lymphoid", "mono. & macro."="Myeloid", "FAPs"="FAPs", "vascular"="Vascular", "B" = "B"),
                              selected = "Lymphoid"))
    ), # column
    
    column(width = 12,
           style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
           column(width = 6,
                  style='padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                  tabsetPanel(id = "inTabset",
                              tabPanel(title = "FeaturePlot", 
                                       column(12,uiOutput("ui_FeaturePlot_low"))),
                              tabPanel(title = "Violin Plot",
                                       column(12,uiOutput("ui_ViolinPlot_low"))),
                              tabPanel(title = "Marker",
                                       column(12,uiOutput("ui_DT_low"))))),
           column(width = 6,
                  style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                  tabsetPanel(id = "inTabset",
                              tabPanel(title = "FeaturePlot", 
                                       column(12,uiOutput("ui_FeaturePlot_high"))),
                              tabPanel(title = "Violin Plot",
                                       column(12,uiOutput("ui_ViolinPlot_high"))),
                              tabPanel(title = "Marker",
                                       column(12,uiOutput("ui_DT_high")))))
    ) # column
  
  ) # fluidRow
))



