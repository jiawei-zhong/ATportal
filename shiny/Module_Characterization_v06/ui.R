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
  tags$head(
    tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/iframe-resizer/4.3.1/iframeResizer.contentWindow.min.js"),
    tags$style('
      .col-sm-3 { padding-left: 0px; padding-right: 15px; }
      .col-sm-9 { padding-left: 15px; padding-right: 0px; }
    ')
  ),
  #titlePanel("Module III: Characterization"),
  sidebarLayout(
    div(class="sidebar",
        sidebarPanel(
          width=3,div(id="side-tab",
                      tabsetPanel(
                        tabPanel(title = "Input",     
                                 br(),
                                 # div("Upload a csv file:"),
                                 #  fileInput("file1", label=NULL,
                                 #           multiple = FALSE,
                                 #            accept = c("text/csv",
                                 #                     "text/comma-separated-values,text/plain",
                                 #                     ".csv")),
                                 # div("Fill in your genes here:", class="highlight"),
                                 # br(),
                                 # tags$textarea(id = 'gene_id', label = "Input your genes of interest", placeholder = c("PPARG","\nADIPOQ"), rows = 8),
                                 # tags$hr(),
                                 textAreaInput(inputId = "gene_id",
                                               label = "Input your genes of interest",
                                               value = "",
                                               # width = "300px",
                                               height = "200px",
                                               placeholder = "PPARG\nADIPOQ"), # nolint
                                 # br(),
                                 
                                 conditionalPanel(
                                   condition = "input.cell_type === 'adsc'",
                                   materialSwitch(inputId = "prot_switch_adsc", label = "Add proteomics plots", status="primary"),
                                   # br(),
                                   selectInput(inputId = "cell_line", label = "Select adipose tissue derived stem cell line", choices = c("hAT MSC", "SGBS")),
                                 ),

                                 conditionalPanel(
                                   condition = "input.cell_type === 'fraction'",
                                   materialSwitch(inputId = "prot_switch_fraction", label = "Add proteomics plots", status="primary"),
                                   # br(),
                                 ),

                                
                                 #div("Currently only available for Adipogensis: hADSC"),
                                 # br(),
                                 # div("Select data type", class="highlight"),
                                 
                                 # br(),
                                 radioButtons(inputId = "mode", label = "Select data type", choiceNames = c("mean abundance", "normalized mean abundance"), 
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
                                 div("See", a(href="https://colorbrewer2.org/", target = "_blank", "ColorBrewer2"), "for more"),
                                 tags$br(),
                                 div("Clustering method for heatmaps", class="highlight"),
                                 tags$br(),
                                 selectInput(inputId= "cluster_method", label=NULL, choices = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"), selected = "ward.D")
                        )
                      )))),
        mainPanel(
          width = 9,
          navbarPage("Cell type sub-modules", id ="cell_type",
                     tabPanel(div("Tissue comparison",class="tab_head"),value="tissue",
                              includeHTML("htmls/description_tissue.html"),
                              # div(actionButton("toggle_table", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                              # div(id= "gene_table", style = "display:none;",includeHTML("htmls/outputs_tissue.html")),   
                              br(),
                              tabsetPanel(
                                tabPanel("Rug plot",
                                         br(),
                                         plotOutput(outputId = "heatmap_fantom"),
                                         br(),
                                         downloadButton(outputId = "download_pdf8", label = "Download PDF", class = "butt"), 
                                         downloadButton(outputId = "download_gg8", label = "Download ggplot2", class = "butt"),
                                         tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                tabPanel("Table",
                                         br(),
                                         DT::dataTableOutput("table_fantom")),
                                tabPanel("Details",
                                         includeHTML("htmls/details_tissue.html")
                                ),
                                tabPanel("Interpretation",
                                         includeHTML("htmls/interpretation_tissue.html"))),
                              tags$br(),
                              tags$hr(),
                              "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/27803022/", target = "_blank", "Ehrlund et al.; PMID: 27803022"),  "(data) and Zhong et al. (portal).",
                              br(),
                              "Raw data for this data set can be downloaded fromthe DNA Data Bank of Japan under accession numbers", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA000991", target = "_blank", "DRA000991"),",", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002711", target = "_blank", "DRA002711"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002747", target = "_blank", "DRA002747"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002748", target = "_blank", "DRA002748"), 
                              
                     ),  
                     tabPanel(div("WAT fractionation", class="tab_head"),value="fraction",
                              includeHTML("htmls/description_novo.html"),
                              # div(actionButton("toggle_table2", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                              # div(id = "gene_table2", style = "display:none;",includeHTML("htmls/outputs_FACS.html")), 
                              #using class = "gene_table" for everything broke the code somehow
                              br(),
                              tabsetPanel(
                                tabPanel("Heatmap",
                                         br(),
                                         # plotOutput(outputId = "heatmap_novo"),
                                         # br(),
                                         # downloadButton(outputId = "download_pdf7", label = "Download PDF", class = "butt"), 
                                         # downloadButton(outputId = "download_gg7", label = "Download ggplot2", class = "butt"),
                                         # tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))
                                         uiOutput(outputId = "ui_novo_heat")
                                         ), 
                                tabPanel("Table",
                                         br(),
                                         DT::dataTableOutput("table_novo")),
                                tabPanel("Details",
                                         includeHTML("htmls/details_novo.html")
                                ),
                                tabPanel("Interpretation",
                                         includeHTML("htmls/interpretation_novo.html")
                                )),
                              tags$br(),
                              tags$hr(),
                              "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/29116032/", "Acosta et al.; PMID: 29116032"),  "(data) and Zhong et al. (portal).",
                              br(),
                              "Raw data for this data set can be downloaded from the GEO database:", tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", target = "_blank", "GSE100795")
                     ),
                     navbarMenu("Adipogenesis",
                                tabPanel("ADSCs", value="adsc",
                                         div("adipose tissue derived stem cells", class="highlight"),
                                         includeHTML("htmls/description_ADSCs.html"), 
                                         # div(actionButton("toggle_table3", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                                         # div(id= "gene_table3", style = "display:none;",includeHTML("htmls/outputs_hADSC.html")),
                                         br(),
                                         tabsetPanel(
                                           
                                           tabPanel("Line chart",
                                                    br(),
                                                    # fluidRow(
                                                    #   splitLayout(cellWidths = c("50%", "50%"), plotOutput(outputId = "linechart"), plotOutput("linechart_prot"))
                                                    # ),
                                                    # br(),
                                                    # downloadButton(outputId = "download_pdf1", label = "Download PDF", class = "butt"),
                                                    # downloadButton(outputId = "download_gg1", label = "Download ggplot2", class = "butt"),
                                                    # tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))
                                                    uiOutput(outputId = "ui_niga_line")
                                                    ), 
                                           tabPanel("Heatmap",
                                                    br(),
                                                    # fluidRow(
                                                    #   splitLayout(cellWidths = c("50%", "50%"), plotOutput("heatmap"), plotOutput("heatmap_prot"))
                                                    # ),
                                                    # br(),
                                                    # downloadButton(outputId = "download_pdf2", label = "Download PDF", class = "butt"), 
                                                    # downloadButton(outputId = "download_gg2", label = "Download ggplot2", class = "butt"),
                                                    # tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))
                                                    uiOutput(outputId = "ui_niga_heat")
                                                    ), 
                                           tabPanel("Table",
                                                    br(),
                                                    DT::dataTableOutput(outputId = "DT_NIGA", width = "100%"),
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
                                         "Raw data for this data set can be downloaded from the DNA Data Bank of Japan under accession numbers", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA000991", target = "_blank", "DRA000991"),",", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002711", target = "_blank", "DRA002711"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002747", target = "_blank", "DRA002747"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002748", target = "_blank", "DRA002748"), 
                                ),
                                tabPanel("SVF", value="svf",
                                         div("Ex vivo stroma vascular fraction", class="highlight"),
                                         includeHTML("htmls/description_SVF.html"),
                                         # div(actionButton("toggle_table4", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                                         # div(id = "gene_table4", style = "display:none;",includeHTML("htmls/outputs_SVF.html")),
                                         br(),
                                         tabsetPanel(
                                           tabPanel("Line chart",
                                                    br(),
                                                    plotOutput(outputId = "linechart2"), br(),
                                                    downloadButton(outputId = "download_pdf3", label = "Download PDF", class = "butt"),
                                                    downloadButton(outputId = "download_gg3", label = "Download ggplot2", class = "butt"),
                                                    tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                           tabPanel("Heatmap",
                                                    br(),
                                                    plotOutput("heatmap2", width="80%"),
                                                    br(),
                                                    downloadButton(outputId = "download_pdf4", label = "Download PDF", class = "butt"),
                                                    downloadButton(outputId = "download_gg4", label = "Download ggplot2", class = "butt"),
                                                    tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                           tabPanel("Table",
                                                    br(),
                                                    DT::dataTableOutput(outputId = "DT_SVF", width = "100%"),
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
                                         "Raw data for this data set can be downloaded from the DNA Data Bank of Japan under accession numbers", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA000991", target = "_blank", "DRA000991"),",", tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002711", target = "_blank", "DRA002711"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002747", target = "_blank", "DRA002747"),",",tags$a(href="https://ddbj.nig.ac.jp/resource/sra-submission/DRA002748", target = "_blank", "DRA002748") 
                                ))))
    
    
  )
)
)#end shiny UI
