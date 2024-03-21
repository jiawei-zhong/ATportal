suppressMessages(library(shiny))
suppressMessages(library(shinyhelper))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))  # for box()
suppressMessages(library(shinycustomloader))
suppressMessages(library(shinyjs))
suppressMessages(library(fresh))
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

shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style_def.css"),
  ),
  use_googlefont("Red Hat Display"),
  use_theme(create_theme(
    theme = "default",
    bs_vars_font(
      family_sans_serif = "'Red Hat Display', cursive"
    )
  )),
  div(class="container" , #important to use the same class ids
      br(), 
      br(),
      sidebarLayout(
        div(class="sidebar",
            sidebarPanel(
              width=3,
              div(id="side-tab",
                  tabsetPanel(#id = 'abcd',
                    tabPanel("Input",     
                             br(),
                             div("Cohorts", class="highlight"),
                             br(),
                             pickerInput(inputId = "cohort_dc",
                                         label = NULL,
                                         choices = list("Keller, M. (2017)"="Keller_et_al",
                                                        "Krieg, L. (2021)"="Krieg_et_al",
                                                        "Arner, P. (2016)"="EMIFeset",
                                                        "Schleinitz, D. (2020)"="Schleinitz_et_al",
                                                        "Mazaki-Tovi, S. (2016)"='Mazaki_et_al',
                                                        "Hoggard, N. (2012)"="Hoggard_et_al"),
                                         multiple = TRUE,
                                         options = pickerOptions(actionsBox = T, size = 10, liveSearch = T)),
                             #br(),
                             hr(),
                             #br(),
                             div("Gene Input", class="highlight"), #use class="highlight" for headings
                             br(),
                             textAreaInput(inputId = "gene_dc",
                                           label = NULL,
                                           value = "",
                                           # width = "300px",
                                           height = "200px",
                                           placeholder = "One per line, eg:\nPPARG\nADIPOQ"), # nolint
                             #tags$textarea(), #put your inputs here
                             hr(),        #use horizontal line to separate input options
                             div("Adipose Depot 1", class="highlight"), #use class="highlight" for headings
                             br(),             #additional empty lines for readability
                             pickerInput(inputId = "depot_ref_dc",
                                         label = NULL,
                                         choices = c('SAT Abdomen','VAT Omentum'),
                                         selected = 'SAT Abdomen',
                                         # width = "300px",
                                         options = pickerOptions(actionsBox = T, size = 10, maxOptions = 1,liveSearch = T,),
                                         multiple = FALSE),
                             hr(),        #use horizontal line to separate input options
                             div("Adipose Depot 2", class="highlight"), #use class="highlight" for headings
                             br(),             #additional empty lines for readability
                             pickerInput(inputId = "depot_qry_dc",
                                         label = NULL,
                                         choices = c('SAT Abdomen','VAT Omentum'),
                                         selected = 'VAT Omentum',
                                         # width = "300px",
                                         options = pickerOptions(actionsBox = T, size = 10, maxOptions = 1,liveSearch = T,),
                                         multiple = FALSE),
                             hr(),        #use horizontal line to separate input options
                             div("Correlation method", class="highlight"), #use class="highlight" for headings
                             br(),             #additional empty lines for readability
                             prettyRadioButtons(inputId = "correlation_method",
                                                label = NULL,
                                                selected = "spearman",
                                                choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                                                inline = TRUE),
                             hr(),        #use horizontal line to separate input options
                             div("Effect Size", class="highlight"), #use class="highlight" for headings
                             br(),             #additional empty lines for readability
                             prettyRadioButtons(inputId = "effect_size",
                                                label = NULL,
                                                selected = "correlation",
                                                choices = c("Correlation" = "correlation", "Standard Mean Difference" = "SMD"),
                                                inline = TRUE),
                             #br(),
                             tags$hr(),
                             #br(),
                             actionButton(inputId = "start", label = "Visualize!")),
                    tabPanel("Customization",                      #I keep this as is for the template
                             #####
                             tags$br(),
                             div("You can use default settings or customize visualizations below:"),
                             tags$br(),
                             div("Gradient colors", class="highlight"),
                             div("applies to heatmaps"),
                             tags$br(),
                             radioButtons(inputId="gradient_col", choiceNames =  c("Purple to green", "Red to Blue" , "Viridis" ), choiceValues = c("PRGn", "RdBu", "viridis"), selected = "PRGn", label = NULL),
                             tags$br(),
                             div("Discrete colors",class="highlight"),
                             div("applies to line and bar charts"),
                             tags$br(),
                             textInput(inputId="discrete_col", value="Dark2", label=NULL, placeholder = "colorbrewer palettes"),
                             div("See", a(href="https://colorbrewer2.org/", "ColorBrewer2"), "for all options"),
                             tags$br(),
                             div("Clustering method for heatmaps", class="highlight"),
                             tags$br(),
                             selectInput(inputId= "cluster_method", label=NULL, choices = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"), selected = "ward.D")
                             #####
                    )
                  )))),
        div(class="main-content",
            mainPanel(
              width=9,
              navbarPage("Depots", 
                         navbarMenu("Depot Association",
                                    tabPanel("Transcriptome",
                                             div("The depot association module provides an overview of expression differences between various WAT depots. It encompasses a total of 222 subjects spanning 6 studies and includes paired samples from abdominal subcutaneous and omental tissues.",class="detail-text"), #keep class, defines block writing
                                             br(),
                                             tabsetPanel(
                                               tabPanel("Plots", #you can add tabPanels for different plot types if needed, or put multiple plots in the first
                                                        br(),
                                                        plotOutput(outputId = "plot_dc")
                                               ),
                                               tabPanel("Tables",
                                                        br(),
                                                        DT::dataTableOutput(outputId = 'data_table_dc')
                                               ),
                                               tabPanel("Details",  #we will put Method details from the original manuscripts here
                                                        br(),
                                                        div("Method or cohort details",class="detail-text")
                                               ),
                                               tabPanel("Interpretation",
                                                        br(),
                                                        div("How to interpret results", class="detail-text")
                                               )
                                             ),
                                             tags$br(),
                                             tags$hr(),
                                             uiOutput(outputId = "ui_citation")
                                             )
                         )
              )
            )
        )
      )
  )#end div class container
))
#end shiny UI