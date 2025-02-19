suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyhelper))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))
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
suppressMessages(library(extrafont))


shinyUI(fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style_def.css"),
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
    tags$style('
      .col-sm-3 { padding-left: 0px; padding-right: 15px; }
      .col-sm-9 { padding-left: 15px; padding-right: 0px; }
    ')
  ),
  sidebarLayout(
    div(class="sidebar",
        sidebarPanel(
          width=3,
          div(id="side-tab",
              tabsetPanel(#id = 'abcd',
                tabPanel("Input",
                         br(),
                         pickerInput(inputId = "cohort_dc",
                                     label = "Cohort",
                                     choices = list("Keller, M. (2017)"="Keller_et_al",
                                                    "Krieg, L. (2021)"="Krieg_et_al",
                                                    "Arner, P. (2016)"="EMIFeset",
                                                    "Schleinitz, D. (2020)"="Schleinitz_et_al",
                                                    "Mazaki-Tovi, S. (2016)"='Mazaki_et_al',
                                                    "Hoggard, N. (2012)"="Hoggard_et_al",
                                                    "MacLaren, RE. (2010)"="GSE15524eset",
                                                    "Salcedo-Tacuma, D. (2022)"="GSE188799eset",
                                                    "Hardy, OT. (2011)"="GSE20950eset",
                                                    "Du Plessis, J. (2015)"="GSE58979eset"),
                                     multiple = TRUE,
                                     options = pickerOptions(actionsBox = T, size = 10, liveSearch = T)),
                         textAreaInput(inputId = "gene_dc",
                                       label = "Gene",
                                       value = "",
                                       # width = "300px",
                                       height = "200px",
                                       placeholder = "One per line, eg:\nPPARG\nADIPOQ"), # nolint
                         pickerInput(inputId = "depot_ref_dc",
                                     label = "Adipose Depot 1",
                                     choices = c('SAT Abdomen','VAT Omentum'),
                                     selected = 'SAT Abdomen',
                                     # width = "300px",
                                     options = pickerOptions(actionsBox = T, size = 10, maxOptions = 1,liveSearch = T,),
                                     multiple = FALSE),
                         pickerInput(inputId = "depot_qry_dc",
                                     label = "Adipose Depot 2",
                                     choices = c('SAT Abdomen','VAT Omentum'),
                                     selected = 'VAT Omentum',
                                     # width = "300px",
                                     options = pickerOptions(actionsBox = T, size = 10, maxOptions = 1,liveSearch = T,),
                                     multiple = FALSE),
                         # prettyRadioButtons(inputId = "effect_size",
                         #                    label = "Effect Size",
                         #                    selected = "SMD",
                         #                    choices = c("Correlation" = "correlation", "Standard Mean Difference" = "SMD", "log2 Fold Change" = "log2_FC"),
                         #                    inline = TRUE),
                         # prettyRadioButtons(inputId = "correlation_method",
                         #                    label = "Correlation method",
                         #                    selected = "spearman",
                         #                    choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                         #                    inline = TRUE),
                         uiOutput(outputId = "ui_effect"),
                         uiOutput(outputId = "ui_corr_method"),
                         #br(),
                         #tags$hr(),
                         #br(),
                         actionButton(inputId = "start", label = "Submit")),
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
          width = 9,
          navbarPage("Depots",id = 'navbar_id_depot',
                     tabPanel("Transcriptomics",
                              div("The depot association module provides an overview of expression differences between various WAT depots. All included subjects have samples form multiple depots allowing for an enxlusively paired data analysis encompassing direct intrasubject comparisons.",class="detail-text"), #keep class, defines block writing
                              br(),
                              tabsetPanel(
                                tabPanel("Plots", #you can add tabPanels for different plot types if needed, or put multiple plots in the first
                                         br(),
                                         plotOutput(outputId = "plot_dc_1"),
                                         downloadButton(outputId = "dc_pdf_dl1",label="Download")
                                ),
                                tabPanel("Tables",
                                         br(),
                                         DT::dataTableOutput(outputId = 'data_table_dc_1')
                                ),
                                tabPanel("Details",  #we will put Method details from the original manuscripts here
                                         div("Method or cohort details",class="detail-text")
                                ),
                                tabPanel("Interpretation",
                                         br(),
                                         "When only one cohort is selected, a single gene can be accepted in the Input field. Pressing the Visualize button will generate a boxplot comparing the expression levels of the given gene between two different depots with paired data points being indicated by a connecting dashed line. An exact p-value, derived from a two-sided Student’s t-test, is specified above the boxes. Mean expression values for each included depot can be acquired form the Table tab.",
                                         tags$br(),tags$br(),
                                         "Selecting multiple cohorts with a single gene results in the generation of a meta-analysis forest plot highlighting standard mean difference or correlation with a 95% confidence interval between the two depots being compared. Weighted average standard mean difference is calculated according to the common or random effects model, and the respective statistics should be chosen based on the heterogeneity testing below the plot.",
                                         tags$br(),tags$br(),
                                         "If multiple cohorts and multiple genes are selected, a heatmap will be generated displaying standard mean difference, correlation, or log2 fold change between each depot in each cohort. Significant values are indicated by asterisks, whereby: ***; p-value < 0.001, **; p-value < 0.01, *; p-value < 0.05."
                                )
                              )
                     ),
                     tabPanel("Proteomics",
                              div("The depot association module provides an overview of expression differences between various WAT depots. All included subjects have samples form multiple depots allowing for an enxlusively paired data analysis encompassing direct intrasubject comparisons.",class="detail-text"), #keep class, defines block writing
                              br(),
                              tabsetPanel(
                                tabPanel("Plots", #you can add tabPanels for different plot types if needed, or put multiple plots in the first
                                         br(),
                                         plotOutput(outputId = "plot_dc_2"),
                                         downloadButton(outputId = "dc_pdf_dl2",label="Download")
                                ),
                                tabPanel("Tables",
                                         br(),
                                         DT::dataTableOutput(outputId = 'data_table_dc_2')
                                ),
                                tabPanel("Details",  #we will put Method details from the original manuscripts here
                                         br(),
                                         div("Method or cohort details",class="detail-text")
                                ),
                                tabPanel("Interpretation",
                                         br(),
                                         "When only one cohort is selected, a single gene can be accepted in the Input field. Pressing the Visualize button will generate a boxplot comparing the expression levels of the given gene between two different depots with paired data points being indicated by a connecting dashed line. An exact p-value, derived from a one-sample, two-sided Student’s t-test comparing to a null hypothesis of a difference of 0, is specified above the boxes. Mean expression values for each included depot can be acquired form the Table tab.",
                                         tags$br(),tags$br(),
                                         "Selecting multiple cohorts with a single gene results in the generation of a meta-analysis forest plot highlighting standard mean difference or correlation with a 95% confidence interval between the two depots being compared. Weighted average standard mean difference is calculated according to the common or random effects model, and the respective statistics should be chosen based on the heterogeneity testing below the plot.",
                                         tags$br(),tags$br(),
                                         "If multiple cohorts and multiple genes are selected, a heatmap will be generated displaying standard mean difference, correlation, or log2 fold change between each depot in each cohort. Significant values are indicated by asterisks, whereby: ***; p-value < 0.001, **; p-value < 0.01, *; p-value < 0.05."
                                )
                              )
                     )
          ),
          tags$br(),
          tags$hr(),
          uiOutput(outputId = "ui_citation")
        )
    )
  )
))
#end shiny UI
