library(shiny)
library(plotly)
library(shinyjs)
library(shinydashboard)
library(fresh)
library(shinyWidgets)
library(DT)

shinyUI(fluidPage(
  useShinyjs(),
  tags$head(
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
  #titlePanel("Module III: Characterization"),
  sidebarLayout(
    div(class="sidebar",
    sidebarPanel(
      width=3,div(id="side-tab",
      tabsetPanel(
          tabPanel("Input",     
      br(),
      div("Fill in your genes here:", class="highlight"),
      br(),
      tags$textarea(id = 'gene_id', label = "Input your genes of interest", placeholder = c("PPARG","\nADIPOQ"), rows = 8),
      tags$hr(),
      br(),
        conditionalPanel(
          condition = "input.perturbation === 'ko'",
          pickerInput(inputId= "ko_gene", label="Select knock out/ knock down cell line", choices = ko_inputs, multiple=FALSE),
        ),
      br(),
      div("Select data type", class="highlight"),

     br(),
    radioButtons(inputId = "mode", label = NULL, choiceNames = c("mean abundance", "normalized mean abundance"), 
                 choiceValues = c("raw","normalized"), selected = "normalized"),
      actionButton(inputId = "start", "Visualize!"),
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
      navbarPage("Perturbation", id ="perturbation",
                 navbarMenu("Cell stress",
                            tabPanel("Hypoxia", value="hypoxia",
                                     div("adipose tissue derived stem cells", class="highlight"),
                                     includeHTML("htmls/description_ADSCs.html"), 
                                     div(actionButton("toggle_table3", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                                     div(id= "gene_table3", style = "display:none;",includeHTML("htmls/outputs_hADSC.html")),
                                     tabsetPanel(
                                       
                                       tabPanel("Heatmap",
                                                br(),
                                                plotOutput(outputId = "heatmap_hyp"),
                                                br(),
                                                downloadButton(outputId = "download_pdf1", label = "Download PDF", class = "butt"),
                                                downloadButton(outputId = "download_gg1", label = "Download ggplot2", class = "butt"),
                                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                       tabPanel("Volcano",
                                                br(),
                                                shinycustomloader::withLoader(plotOutput(outputId = "volcano_hyp"), type="html", loader="dnaspin"),
                                                br(),
                                                downloadButton(outputId = "download_pdf2", label = "Download PDF", class = "butt"), 
                                                downloadButton(outputId = "download_gg2", label = "Download ggplot2", class = "butt"),
                                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                       tabPanel("Table",
                                                br(),
                                                DT::dataTableOutput(outputId = "dt_hyp", width = "100%"),
                                                br()),
                                       tabPanel("Details",
                                                includeHTML("htmls/details_ADSCs.html")
                                       ),
                                       tabPanel("Interpretation",
                                                includeHTML("htmls/interpretation_ADSCs.html")
                                       )),
                                     tags$br(),
                                     tags$hr(),
                                     "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/27803022/", "Ehrlund et al.; PMID: 27803022") ,  "(data) and Zhong et al. (portal).",
                                     br(),
                                     "Raw data for this data set can be downloaded from the DNA Data Bank of Japan under accession numbers", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA000991", "DRA000991"),",", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002711", "DRA002711"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002747", "DRA002747"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002748", "DRA002748"), 
                            ),
                            tabPanel("TNF", value="tnf",
                                     div("Ex vivo stroma vascular fraction", class="highlight"),
                                     includeHTML("htmls/description_SVF.html"),
                                     div(actionButton("toggle_table4", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                                     div(id = "gene_table4", style = "display:none;",includeHTML("htmls/outputs_SVF.html")),
                                     tabsetPanel(
                                       tabPanel("Heatmap",
                                                br(),
                                                plotOutput(outputId = "heatmap_tnf"), br(),
                                                downloadButton(outputId = "download_pdf3", label = "Download PDF", class = "butt"),
                                                downloadButton(outputId = "download_gg3", label = "Download ggplot2", class = "butt"),
                                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                       tabPanel("Volcano",
                                                br(),
                                                shinycustomloader::withLoader(plotOutput(outputId = "volcano_tnf"), type="html", loader="dnaspin"),
                                                br(),
                                                downloadButton(outputId = "download_pdf4", label = "Download PDF", class = "butt"),
                                                downloadButton(outputId = "download_gg4", label = "Download ggplot2", class = "butt"),
                                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                       tabPanel("Table",
                                                br(),
                                                DT::dataTableOutput(outputId = "dt_tnf", width = "100%"),
                                                br()),
                                       tabPanel("Details",
                                                includeHTML("htmls/details_SVF.html")
                                       ),
                                       tabPanel("Interpretation",
                                                includeHTML("htmls/interpretation_SVF.html")
                                       )),
                                     tags$br(),
                                     tags$hr(),
                                     "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/27803022/", "Ehrlund et al.; PMID: 27803022"),  "(data) and Zhong et al. (portal).",
                                     br(),
                                     "Raw data for this data set can be downloaded from the DNA Data Bank of Japan under accession numbers", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA000991", "DRA000991"),",", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002711", "DRA002711"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002747", "DRA002747"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002748", "DRA002748") 
                            )),
                 tabPanel(div("Insulin stimulation",class="tab_head"),value="insulin",
                          includeHTML("htmls/description_tissue.html"),
                          div(actionButton("toggle_table", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                          div(id= "gene_table", style = "display:none;",includeHTML("htmls/outputs_tissue.html")),   
                                     br(),
                                     tabsetPanel(
                                       tabPanel("Heatmap",
                                                br(),
                                                plotOutput(outputId = "heatmap_ins"),
                                                br(),
                                                downloadButton(outputId = "download_pdf8", label = "Download PDF", class = "butt"), 
                                                downloadButton(outputId = "download_gg8", label = "Download ggplot2", class = "butt"),
                                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                       tabPanel("Table",
                                                br(),
                                                DT::dataTableOutput("dt_ins")),
                                       tabPanel("Volcano",
                                                br(),
                                                shinycustomloader::withLoader(plotOutput(outputId = "volcano_ins"), type="html", loader="dnaspin"),
                                                br(),
                                                downloadButton(outputId = "download_pdf9", label = "Download PDF", class = "butt"),
                                                downloadButton(outputId = "download_gg9", label = "Download ggplot2", class = "butt"),
                                                tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                       tabPanel("Details",
                                                includeHTML("htmls/details_tissue.html")
                                       ),
                                       tabPanel("Interpretation",
                                                includeHTML("htmls/interpretation_tissue.html"))),
                                     tags$br(),
                                     tags$hr(),
                                     "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/27803022/", "Ehrlund et al.; PMID: 27803022"),  "(data) and Zhong et al. (portal).",
                                     br(),
                                     "Raw data for this data set can be downloaded fromthe DNA Data Bank of Japan under accession numbers", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA000991", "DRA000991"),",", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002711", "DRA002711"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002747", "DRA002747"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002748", "DRA002748"), 
                                     
                            ),  
             tabPanel(div("knock out/down", class="tab_head"),value="ko",
                      uiOutput("knockdowns")       
                        ),
)))
  
   
)
)#end div class container
)
)#end shiny UI
