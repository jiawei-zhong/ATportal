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
suppressMessages(library(shinyjs))

shinyUI(fluidPage(
  useShinyjs(),
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
  div(class="container" ,
      br(), 
      br(),
      sidebarLayout(
        div(class="sidebar",
            sidebarPanel(
              width=3,div(id="side-tab",
                          tabsetPanel(
                            tabPanel("Input",     
                                     #br(),
                                     div("Please select a single gene of interest or a cell type below", class="highlight"),
                                     #br(),
                                     #selectInput(inputId= "scseq_cohort", label="Cohort: ", choices = c("Massier et al.","Hinte et al.","Wang et al."), selected = "Massier et al."),
                                     br(),
                                     selectizeInput(inputId = 'gene_for_featureplot_STx',
                                                    label = 'Gene: ',
                                                    choices = NULL,
                                                    # selected = "LEP",
                                                    width = "200px",
                                                    # inline = TRUE,
                                                    multiple = FALSE),
                                     tags$hr(),
                                     selectizeInput(inputId = 'deconvolution_for_featureplot_STx',
                                                    label = 'Cell type: ',
                                                    choices = NULL,
                                                    # selected = "LEP",
                                                    width = "200px",
                                                    # inline = TRUE,
                                                    multiple = FALSE),
                                     br(),
                                     selectInput(inputId = "Dataset_STx",
                                                 label = "Dataset: ",
                                                 width = "200px",
                                                 choices = c("Bäckdahl et al. baseline"="Jesper_et_al_baseline", "Bäckdahl et al. insulin"="Jesper_et_al_insulin"),
                                                 selected = "Jesper_et_al_baseline"),
                                     
                                     br(),
                                     selectInput(inputId = "Slide_STx",
                                                 label = "Slide: ",
                                                 width = "200px",
                                                 choices = NA,
                                                 selected = NULL),
                                     br(),
                                     actionButton(inputId = "start", "Visualize!"),
                                     tags$br()),
                            tabPanel("Customization"
                                     # tags$br(),
                                     # div("You can use default settings or customize visualizations below:"),
                                     # tags$br(),
                                     # div("Gradient colors", class="highlight"),
                                     # div("applies to heatmaps"),
                                     # tags$br(),
                                     # radioButtons(inputId="gradient_col", choiceNames =  c("Default", "Purple to green", "Red to Blue" , "Viridis" ), choiceValues = c("default", "PRGn", "RdBu", "viridis"), selected = "default", label = NULL),
                                     # tags$br(),
                                     # div("Discrete colors",class="highlight"),
                                     # div("applies to line and bar charts"),
                                     # tags$br(),
                                     # radioButtons(inputId="discrete_col_sel", choiceNames =  c("Default", "Colorbrewer 2"), choiceValues = c("default_2", "CB2"), selected = "default_2", label = NULL),
                                     # div("if you selected ColorBrewer 2, please select a colorscheme below"),
                                     # selectInput(inputId= "discrete_col", label=NULL, choices = c("Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3","Spectral",
                                     #                                                              "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG"), selected = "Dark2"),
                                     # div("See", a(href="https://colorbrewer2.org/", "ColorBrewer2"), "for more"),
                                     # tags$br(),
                                     # div("Clustering method for heatmaps", class="highlight"),
                                     # tags$br(),
                                     # selectInput(inputId= "cluster_method", label=NULL, choices = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"), selected = "ward.D")
                            )
                          )))), # end sidebar
        div(class="main-content",
            mainPanel(
              navbarPage("Spatial", id ="spatial",
                         tabPanel("Baseline", value="scseq_base",
                                             includeHTML("htmls/description_baseline.html"), 
                                             tabsetPanel(
                                               
                                               tabPanel("Spatial Feature Plot",
                                                        #content Spatial FeaturePlot
                                                        uiOutput("ui_SpatialFeaturePlot_STx"),
                                                        tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                               tabPanel("Feature Plot",
                                                        uiOutput("ui_FeaturePlot_STx"),
                                                        #content Feature Plot
                                                        tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                               tabPanel("Violin Plot",
                                                        #content Violin Plot
                                                        uiOutput("ui_ViolinPlot_STx"),
                                                        tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                               tabPanel("Marker genes",
                                                        #content all unfiltered marker genes
                                                        uiOutput("ui_DT_STx")
                                                        ),
                                               tabPanel("Details",
                                                        includeHTML("htmls/details_baseline.html")
                                               ),
                                               tabPanel("Interpretation",
                                                        includeHTML("htmls/interpretation_baseline.html")
                                               )),
                                             tags$br(),
                                             tags$hr(),
                                             includeHTML("htmls/reference_baseline.html") 
                                    ), #end tabPanel
                         ), #end navbarPage
              )) #end div and mainPanel
      ) #end sidebar layout
        
      ) #end div class container
)
)#end shiny UI
