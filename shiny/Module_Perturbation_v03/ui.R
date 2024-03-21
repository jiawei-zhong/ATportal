library(shiny)
library(plotly)
library(shinyjs)
library(shinydashboard)
library(fresh)

# Define UI for application that draws a histogram
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
  div(class="container" ,
      br(), 
      br(),
      #titlePanel("Perturbation"),
    sidebarLayout(
      div(class="sidebar",
      sidebarPanel(
        width=3, 
        div(id="side-tab",
            tabsetPanel(
              tabPanel("Input",
                       br(),
                #fileInput("file1", "Choose CSV File",
                #multiple = FALSE,
                #accept = c("text/csv",
                #           "text/comma-separated-values,text/plain",
                #           ".csv")),
                div("Fill in your genes here:", class="highlight"),
                br(),
                tags$textarea(id = 'gene_id', label = "Input your genes of interest", placeholder = c("PPARG","\nADIPOQ"), rows = 8),
                tags$hr(),
                br(),
                div("Select data type", class="highlight"),
                radioButtons(inputId = "mode", label = NULL, choiceNames = c("mean abundance", "normalized mean abundance"), 
                             choiceValues = c("raw","normalized"), selected = "normalized"),
                actionButton(inputId = "start", label = "Visualize!"),
                tags$br()),
              tabPanel("Customization",
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
              )
            )))),
    div(class="main-content", 
    mainPanel(
      width=9,
      navbarPage("Perturbation sub-modules",
                 navbarMenu("Cell stress",
                          tabPanel("Hypoxia",
                          div("Description",class="detail-text"),
                          br(),
                            tabsetPanel(
                              tabPanel("Heatmap",
                                       br(),
                                       plotOutput(outputId = "heatmap_hypoxia"),
                                       br(),
                                       downloadButton(outputId = "download_pdf5", label = "Download PDF"), 
                                       downloadButton(outputId = "download_gg5", label = "Download ggplot2")),
                              tabPanel("Table",
                                       br(),
                                       DT::dataTableOutput(outputId = "DT_hypoxia", width = "100%")),
                              tabPanel("Volcano",
                                       br(),
                                       shinycustomloader::withLoader(plotlyOutput(outputId = "volcano_hypoxia", width = "100%"), type="html", loader="dnaspin"),
                                       br(),
                                       downloadButton(outputId = "download_pdf6", label = "Download PDF"), 
                                       downloadButton(outputId = "download_gg6", label = "Download ggplot2")),
                              tabPanel("Details",
                                       br(),
                                       p(style="text-align: justify;","Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque arcu lacus, molestie id facilisis ut, porta et neque. Nullam id ultricies urna, quis pulvinar tellus. Duis dignissim purus et gravida aliquet. Curabitur vel sapien purus. Etiam nisi mi, cursus vel arcu id, lobortis faucibus lacus. Praesent mauris odio, vulputate eu sapien vitae, maximus mollis erat. Etiam at venenatis lorem, nec cursus tellus. Nulla rutrum consectetur iaculis. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras id lacus at dolor lacinia imperdiet eu quis nunc. Nunc nisi ex, tincidunt ut ultricies eget, porta sed libero.
                              Sed eget velit sed mauris eleifend euismod. Duis a risus ultrices, tristique eros sit amet, tempor sapien. Aliquam ut placerat erat, ut viverra tellus. Integer commodo placerat justo, dignissim viverra elit semper a. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc ornare commodo nisi, vitae facilisis nisi feugiat sed. Quisque imperdiet ultrices eleifend. Vivamus sed fringilla nibh, maximus ultricies erat. ")),
                              tabPanel("Interpretation",
                                       br(),
                                       p(style="text-align: justify;","Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque arcu lacus, molestie id facilisis ut, porta et neque. Nullam id ultricies urna, quis pulvinar tellus. Duis dignissim purus et gravida aliquet. Curabitur vel sapien purus. Etiam nisi mi, cursus vel arcu id, lobortis faucibus lacus. Praesent mauris odio, vulputate eu sapien vitae, maximus mollis erat. Etiam at venenatis lorem, nec cursus tellus. Nulla rutrum consectetur iaculis. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras id lacus at dolor lacinia imperdiet eu quis nunc. Nunc nisi ex, tincidunt ut ultricies eget, porta sed libero.
                              Sed eget velit sed mauris eleifend euismod. Duis a risus ultrices, tristique eros sit amet, tempor sapien. Aliquam ut placerat erat, ut viverra tellus. Integer commodo placerat justo, dignissim viverra elit semper a. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc ornare commodo nisi, vitae facilisis nisi feugiat sed. Quisque imperdiet ultrices eleifend. Vivamus sed fringilla nibh, maximus ultricies erat. "))
                            ),
                            tags$br(),
                            tags$hr(),
                            "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/27803022/", "Ehrlund et al.; PMID: 27803022"),  "(data) and Zhong et al. (portal).",
                            br(),
                            "Raw data for this data set can be downloaded from: "
                 ),
                 tabPanel("TNF stimulation INT",
                          tabsetPanel(
                            tabPanel("Heatmap",
                                     br(),
                                     plotOutput(outputId = "heatmap_tnf"),
                                     br(),
                                     downloadButton(outputId = "download_pdf7", label = "Download PDF"), 
                                     downloadButton(outputId = "download_gg7", label = "Download ggplot2")),
                            tabPanel("Table",
                                     br(),
                                     DT::dataTableOutput(outputId = "DT_tnf", width = "100%")),
                            tabPanel("Volcano",
                                     br(),
                                     shinycustomloader::withLoader(plotlyOutput(outputId = "volcano_tnf", width = "100%"), type="html", loader="dnaspin"),
                                     br(),
                                     downloadButton(outputId = "download_pdf8", label = "Download PDF"), 
                                     downloadButton(outputId = "download_gg8", label = "Download ggplot2")),
                            tabPanel("Details",
                                     br(),
                                     p(style="text-align: justify;","Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque arcu lacus, molestie id facilisis ut, porta et neque. Nullam id ultricies urna, quis pulvinar tellus. Duis dignissim purus et gravida aliquet. Curabitur vel sapien purus. Etiam nisi mi, cursus vel arcu id, lobortis faucibus lacus. Praesent mauris odio, vulputate eu sapien vitae, maximus mollis erat. Etiam at venenatis lorem, nec cursus tellus. Nulla rutrum consectetur iaculis. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras id lacus at dolor lacinia imperdiet eu quis nunc. Nunc nisi ex, tincidunt ut ultricies eget, porta sed libero.
                              Sed eget velit sed mauris eleifend euismod. Duis a risus ultrices, tristique eros sit amet, tempor sapien. Aliquam ut placerat erat, ut viverra tellus. Integer commodo placerat justo, dignissim viverra elit semper a. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc ornare commodo nisi, vitae facilisis nisi feugiat sed. Quisque imperdiet ultrices eleifend. Vivamus sed fringilla nibh, maximus ultricies erat. ")),
                            tabPanel("Interpretation",
                                     br(),
                                     p(style="text-align: justify;","Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque arcu lacus, molestie id facilisis ut, porta et neque. Nullam id ultricies urna, quis pulvinar tellus. Duis dignissim purus et gravida aliquet. Curabitur vel sapien purus. Etiam nisi mi, cursus vel arcu id, lobortis faucibus lacus. Praesent mauris odio, vulputate eu sapien vitae, maximus mollis erat. Etiam at venenatis lorem, nec cursus tellus. Nulla rutrum consectetur iaculis. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras id lacus at dolor lacinia imperdiet eu quis nunc. Nunc nisi ex, tincidunt ut ultricies eget, porta sed libero.
                              Sed eget velit sed mauris eleifend euismod. Duis a risus ultrices, tristique eros sit amet, tempor sapien. Aliquam ut placerat erat, ut viverra tellus. Integer commodo placerat justo, dignissim viverra elit semper a. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc ornare commodo nisi, vitae facilisis nisi feugiat sed. Quisque imperdiet ultrices eleifend. Vivamus sed fringilla nibh, maximus ultricies erat. "))
                               ),
                          tags$br(),
                          tags$hr(),
                          "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/27803022/", "Ehrlund et al.; PMID: 27803022"),  "(data) and Zhong et al. (portal).",
                          br(),
                          "Raw data for this data set can be downloaded from: ")),
                 navbarMenu("knock out cell lines",
                            tabPanel("siGLS",
                                     tabsetPanel(
                                       tabPanel("Heatmap", plotOutput(outputId = "heatmap_gls")), 
                                       tabPanel("Table", dataTableOutput("DT_gls")),
                                       tabPanel("Volcano", 
                                                shinycustomloader::withLoader(plotlyOutput("volcano_gls"), type="html", loader="dnaspin")),
                                       tabPanel("Details",
                                                br(),
                                                p(style="text-align: justify;","Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque arcu lacus, molestie id facilisis ut, porta et neque. Nullam id ultricies urna, quis pulvinar tellus. Duis dignissim purus et gravida aliquet. Curabitur vel sapien purus. Etiam nisi mi, cursus vel arcu id, lobortis faucibus lacus. Praesent mauris odio, vulputate eu sapien vitae, maximus mollis erat. Etiam at venenatis lorem, nec cursus tellus. Nulla rutrum consectetur iaculis. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras id lacus at dolor lacinia imperdiet eu quis nunc. Nunc nisi ex, tincidunt ut ultricies eget, porta sed libero.
                              Sed eget velit sed mauris eleifend euismod. Duis a risus ultrices, tristique eros sit amet, tempor sapien. Aliquam ut placerat erat, ut viverra tellus. Integer commodo placerat justo, dignissim viverra elit semper a. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc ornare commodo nisi, vitae facilisis nisi feugiat sed. Quisque imperdiet ultrices eleifend. Vivamus sed fringilla nibh, maximus ultricies erat. ")),
                                       tabPanel("Interpretation",
                                                br(),
                                                p(style="text-align: justify;","Lorem ipsum dolor sit amet, consectetur adipiscing elit. Quisque arcu lacus, molestie id facilisis ut, porta et neque. Nullam id ultricies urna, quis pulvinar tellus. Duis dignissim purus et gravida aliquet. Curabitur vel sapien purus. Etiam nisi mi, cursus vel arcu id, lobortis faucibus lacus. Praesent mauris odio, vulputate eu sapien vitae, maximus mollis erat. Etiam at venenatis lorem, nec cursus tellus. Nulla rutrum consectetur iaculis. Orci varius natoque penatibus et magnis dis parturient montes, nascetur ridiculus mus. Cras id lacus at dolor lacinia imperdiet eu quis nunc. Nunc nisi ex, tincidunt ut ultricies eget, porta sed libero.
                              Sed eget velit sed mauris eleifend euismod. Duis a risus ultrices, tristique eros sit amet, tempor sapien. Aliquam ut placerat erat, ut viverra tellus. Integer commodo placerat justo, dignissim viverra elit semper a. Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nunc ornare commodo nisi, vitae facilisis nisi feugiat sed. Quisque imperdiet ultrices eleifend. Vivamus sed fringilla nibh, maximus ultricies erat. ")),
                                       tags$br(),
                                       tags$hr(),
                                       footer =  "If you want to use this figure in your publication, please cite: XXX et al. (data) and XXX et al. (portal).")
                                     )))))

)
    )#end div class container
  )
)#end shiny UI
