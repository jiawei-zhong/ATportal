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

shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style_def.css")
  ),
  # use_googlefont("Red Hat Display"),
  use_theme(create_theme(
    theme = "default",
    bs_vars_font(
      family_sans_serif = "'Red Hat Display', cursive"
    )
  )),
  tags$head(
    tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/iframe-resizer/4.3.1/iframeResizer.contentWindow.min.js"),
    # tags$style('
    #   /* .sidebar {position: sticky; top: 0;} /*
    #   .col-sm-3 { padding-left: 0px; padding-right: 15px; }
    #   .col-sm-9 { padding-left: 15px; padding-right: 0px; }
    # ')
    tags$style('
      .col-sm-3 { padding-left: 0px; padding-right: 15px; }
      .col-sm-9 { padding-left: 15px; padding-right: 0px; }
    ')
  ),
  tags$script(HTML("
    $(document).ready(function() {
        window.addEventListener('message', function(event) {
            if (event.origin !== 'https://www.adiposetissue.org/clinical') { // 替换为外部页面的实际来源
                return;
            }
            if (event.data.type === 'scroll') {
                $('.sticky-sidebar').css('top', event.data.scrollTop + 'px');
            } else if (event.data.type === 'resize') {
                $('.sticky-sidebar').css('height', event.data.height + 'px');
            }
        });
    });
  ")),
  
  sidebarLayout(
    div(class="sidebar",
        sidebarPanel(
          width=3,
          div(id="side-tab",
              tabsetPanel(
                tabPanel(title = "Input",     
                         #br(),
                         #div("Please select a single gene of interest", class="highlight"),
                         br(),
                         selectInput(inputId= "Cohort_bl", 
                                     label="Cohort", 
                                     choices = c("Massier, L. (2023)"="Massier_et_al",
                                                 "Hinte, LC. (2024)"="Hinte_et_al",
                                                 "Reinisch, I. (2024)"="Reinisch_et_al"), 
                                     selected = "Massier_et_al"),
                         # br(),
                         selectizeInput(inputId = 'Gene_bl',
                                        label = "Gene",
                                        choices = NULL,
                                        selected = "LEP",
                                        # width = "200px",
                                        # inline = TRUE,
                                        multiple = FALSE),
                         # tags$hr(),
                         selectInput(inputId= "Subcluster_bl", label="Subcluster", choices = NULL, selected = NULL),
                         # br(),
                         # conditional input panels, e.g., for depot selection 
                         #conditionalPanel(
                         #   condition = "input.perturbation === 'ko'",
                         #   pickerInput(inputId= "ko_gene", label="Select knock out/ knock down cell line", choices = ko_inputs, multiple=FALSE),
                         # ),
                         
                         actionButton(inputId = "SearchButton_bl", "Submit"),
                         tags$br()),
                tabPanel("Customization",
                         tags$br(),
                         div("You can use default settings or customize visualizations below:"),
                         tags$br(),
                         div("Gradient colors", class="highlight"),
                         div("applies to heatmaps"),
                         tags$br(),
                         radioButtons(inputId="gradient_col", choiceNames =  c("Default", "Purple to green", "Red to Blue" , "Viridis" ), choiceValues = c("default", "PRGn", "RdBu", "viridis"), selected = "default", label = NULL),
                         tags$br(),
                         div("Discrete colors",class="highlight"),
                         div("applies to line and bar charts"),
                         tags$br(),
                         radioButtons(inputId="discrete_col_sel", choiceNames =  c("Default", "Colorbrewer 2"), choiceValues = c("default_2", "CB2"), selected = "default_2", label = NULL),
                         div("if you selected ColorBrewer 2, please select a colorscheme below"),
                         selectInput(inputId= "discrete_col", label=NULL, choices = c("Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3","Spectral",
                                                                                      "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG"), selected = "Dark2"),
                         div("See", a(href="https://colorbrewer2.org/", "ColorBrewer2"), "for more"),
                         tags$br(),
                         div("Clustering method for heatmaps", class="highlight"),
                         tags$br(),
                         selectInput(inputId= "cluster_method", label=NULL, choices = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"), selected = "ward.D")
                )
              )))),
    div(class="main-content",
        mainPanel(
          width = 9,
          navbarPage("Single-cell", id ="singlecell",
                     tabPanel("Baseline", value="scseq_base",
                              includeHTML("htmls/description_baseline.html"), 
                              tabsetPanel(
                                tabPanel("Feature Plot",
                                         br(),
                                         uiOutput(outputId = "ui_FeaturePlot"),
                                         tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                tabPanel("Violin Plot",
                                         br(),
                                         uiOutput(outputId = "ui_ViolinPlot"),
                                         tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                # tabPanel("Marker",
                                #          br(),
                                #          uiOutput(outputId = "ui_DT"),
                                #          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                # tabPanel("Dot Plot",
                                #          br(),
                                #          uiOutput(outputId = "ui_DotPlot"),
                                #          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                # tabPanel("Density Plot",
                                #          br(),
                                #          uiOutput(outputId = "ui_DensityPlot"),
                                #          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                tabPanel("Details",
                                         includeHTML("htmls/details_baseline.html")
                                ),
                                tabPanel("Interpretation",
                                         includeHTML("htmls/interpretation_baseline.html")
                                )),
                              tags$br(),
                              tags$hr(),
                              uiOutput(outputId = "ui_citation")
                     ),
                     #Clinical
                     # navbarMenu("Clinical",
                     #            tabPanel("Weight loss", value = "scseq_weightloss",
                     #                     includeHTML("htmls/description_weightloss.html"),  
                     #                     
                     #                     tabsetPanel(
                     #                       
                     #                       tabPanel("Feature Plot",
                     #                                #content FeaturePlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #                       tabPanel("Violin Plot",
                     #                                #content ViolinPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #                       tabPanel("Dot Plot",
                     #                                #content DotPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #                       tabPanel("Density Plot",
                     #                                #content DensityPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #                       tabPanel("Details",
                     #                                includeHTML("htmls/details_weightloss.html")
                     #                       ),
                     #                       tabPanel("Interpretation",
                     #                                includeHTML("htmls/interpretation_weightloss.html")
                     #                       )),
                     #                     tags$br(),
                     #                     tags$hr(),
                     #                     includeHTML("htmls/reference_weightloss.html") 
                     #            ),
                     #            tabPanel("Healthy vs Unhealthy Obesity", value = "scseq_healthy",
                     #                     includeHTML("htmls/description_healthy.html"),  
                     #                     
                     #                     tabsetPanel(
                     #                       
                     #                       tabPanel("Feature Plot",
                     #                                #content FeaturePlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #                       tabPanel("Violin Plot",
                     #                                #content ViolinPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #                       tabPanel("Dot Plot",
                     #                                #content DotPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #                       tabPanel("Density Plot",
                     #                                #content DensityPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #                       tabPanel("Details",
                     #                                includeHTML("htmls/details_healthy.html")
                     #                       ),
                     #                       tabPanel("Interpretation",
                     #                                includeHTML("htmls/interpretation_healthy.html")
                     #                       )),
                     #                     tags$br(),
                     #                     tags$hr(),
                     #                     includeHTML("htmls/reference_healthy.html") 
                     #            ),
                     #            tabPanel("Sex", value = "scseq_sex",
                     #                     includeHTML("htmls/description_sex.html"),  
                     #                     
                     #                     tabsetPanel(
                     #                       
                     #                       tabPanel("Feature Plot",
                     #                                #content FeaturePlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #                       tabPanel("Violin Plot",
                     #                                #content ViolinPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #                       tabPanel("Dot Plot",
                     #                                #content DotPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #                       tabPanel("Density Plot",
                     #                                #content DensityPlot
                     #                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #                       tabPanel("Details",
                     #                                includeHTML("htmls/details_sex.html")
                     #                       ),
                     #                       tabPanel("Interpretation",
                     #                                includeHTML("htmls/interpretation_sex.html")
                     #                       )),
                     #                     tags$br(),
                     #                     tags$hr(),
                     #                     includeHTML("htmls/reference_sex.html") 
                     #            )
                     # ),
                     #Depots
                     # tabPanel("Depots", value = "scseq_depots",
                     #          includeHTML("htmls/description_depots.html"),  
                     #          
                     #          tabsetPanel(
                     #            
                     #            tabPanel("Feature Plot",
                     #                     #content FeaturePlot
                     #                     tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #            tabPanel("Violin Plot",
                     #                     #content ViolinPlot
                     #                     tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                     #            tabPanel("Dot Plot",
                     #                     #content DotPlot
                     #                     tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #            tabPanel("Density Plot",
                     #                     #content DensityPlot
                     #                     tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                     #            tabPanel("Details",
                     #                     includeHTML("htmls/details_depots.html")
                     #            ),
                     #            tabPanel("Interpretation",
                     #                     includeHTML("htmls/interpretation_depots.html")
                     #            )),
                     #          tags$br(),
                     #          tags$hr(),
                     #          includeHTML("htmls/reference_depots.html") 
                     # ),
          )
        )
    )
    
  )
)
)#end shiny UI
