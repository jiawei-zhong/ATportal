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
          width = 3,
          div(id="side-tab",
              tabsetPanel(
                tabPanel(title = "Input",
                         conditionalPanel(
                           condition = "input.navbar_selected == 'Phenotype correlation'",
                           br(),
                           pickerInput(inputId = "cohort_pc_RNA",
                                       label = "Cohort",
                                       choices = list(META= c("META sc"="META_sc",
                                                              "META om"="META_om"),
                                                      subcutaneous=c("Kerr, A. (2020) Baseline"="DEOSHeset_Baseline",               #  sc  phenotype
                                                                     "Petrus, P. (2018) Baseline"="POeset_Baseline",                #  sc  phenotype
                                                                     "Arner, P. (2018)"="SOWOTeset",                                #  sc  phenotype
                                                                     "Keller, M. (2017) sc"="Keller_et_al_sc",                      #  sc  phenotype  sex
                                                                     "Arner, E. (2012)"="RIKENeset",                                #  sc  phenotype
                                                                     "Stančáková, A. (2012)"="METSIM_204eset",                      #  sc  phenotype
                                                                     "Raulerson, C. (2019)"="METSIM_434eset",                       #  sc  phenotype
                                                                     "Civelek, M. (2017)"="METSIM_770eset",                         #  sc  phenotype
                                                                     "Krieg, L. (2021) sc"="Krieg_et_al_sc",                        #  sc  phenotype  sex
                                                                     "Arner, P. (2016) sc"="EMIFeset_sc",                           #  sc  phenotype
                                                                     "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",   #  sc  phenotype  sex
                                                                     "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",  #  sc  phenotype  sex
                                                                     "Winnier, DA. (2015)"="GSE64567eset",                          #  sc  phenotype  sex
                                                                     "Nono Nankam, PA. (2020)"="Nankam_et_al",                      #  sc  phenotype
                                                                     "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",           #  sc  phenotype  sex
                                                                     "Salcedo-Tacuma, D. (2022) sc"="GSE188799eset_sc",             #  sc  phenotype
                                                                     "MacLaren, RE. (2010) sc"="GSE15524eset_sc",                   #  sc  phenotype  sex
                                                                     "Matualatupauw, JC. (2017)"="GSE87382eset",                    #  sc  phenotype  sex
                                                                     "Defour, M. (2020)"="GSE154610eset",                           #  sc  phenotype  sex
                                                                     "Johansson, LE. (2012) Baseline"="GSE35411eset_Baseline",      #  sc  phenotype  sex
                                                                     "Du Plessis, J. (2015) sc"="GSE58979eset_sc",                  #  sc  phenotype  sex
                                                                     "Van Bussel, IPG. (2017)"="GSE84046eset_Baseline",             #  sc  phenotype  sex
                                                                     "Grundberg, E. (2012)"="E_TABM_1140",                          #  sc  phenotype
                                                                     "Heinonen, S. (2017)"="GSE92405eset",                          #  sc  phenotype  sex
                                                                     "Sharma, NK. (2016)"="GSE95674eset"                            #  sc  phenotype  sex
                                                                     ),
                                                      omental=c("Keller, M. (2017) om"="Keller_et_al_om",                           #  om  phenotype  sex
                                                                "Krieg, L. (2021) om"="Krieg_et_al_om",                             #  om  phenotype  sex
                                                                "Arner, P. (2016) om"="EMIFeset_om",                                #  om  phenotype
                                                                "Barberio, MD. (2019)"="GSE88837eset",                              #  om  phenotype       ethnicity
                                                                "Salcedo-Tacuma, D. (2022) om"="GSE188799eset_om",                  #  om  phenotype
                                                                "MacLaren, RE. (2010) om"="GSE15524eset_om",                        #  om  phenotype  sex
                                                                "Du Plessis, J. (2015) om"="GSE58979eset_om",                       #  om  phenotype  sex
                                                                "Aguilera, CM. (2015)"="GSE9624eset"                                #  om  phenotype
                                                                )
                                                      ),
                                       # selected = "RIKENeset",
                                       multiple = TRUE,
                                       # width = "300px",
                                       options = pickerOptions(actionsBox = FALSE, size = 10, liveSearch = T),
                                       choicesOpt = list(subtext = c("sc","om",rep("sc",25),rep("om",8)))),
                           uiOutput("ui_warning_pc"), # 动态UI输出显示警告信息
                           pickerInput(inputId = "trait_pc_RNA",
                                       label = "Trait",
                                       choices = trait_vector,
                                       # width = "300px",
                                       options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
                                       multiple = TRUE),
                           uiOutput("ui_correction_pc_RNA"),
                           textAreaInput(inputId = "gene_pc",
                                         label = "Gene",
                                         value = "",
                                         # width = "300px",
                                         height = "200px",
                                         placeholder = "One per line, eg:\nPPARG\nADIPOQ"), # nolint
                           prettyRadioButtons(inputId = "correlation_method_pc",
                                              label = "Correlation method",
                                              selected = "spearman",
                                              choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                                              inline = TRUE),
                           tags$label(style = "font-weight: bold; padding-right: 10px; margin: 0;", "Add proteomics plots", `for` = "proteomics_pc"),
                           materialSwitch(inputId = "proteomics_pc", label = NULL, status="primary"),
                           actionButton(inputId = "SearchButton_pc", label = "Submit")
                         ),
                         conditionalPanel(
                           condition = "input.navbar_selected == 'Sex'",
                           br(),
                           pickerInput(inputId = "cohort_sex",
                                       label = "Cohort",
                                       choices = list(subcutaneous=c("Keller, M. (2017) sc"="Keller_et_al_sc",                            #  sc  phenotype  sex
                                                                     "Krieg, L. (2021) sc"="Krieg_et_al_sc",                              #  sc  phenotype  sex
                                                                     "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",         #  sc  phenotype  sex
                                                                     "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",        #  sc  phenotype  sex
                                                                     "Winnier, DA. (2015)"="GSE64567eset",                                #  sc  phenotype  sex
                                                                     "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",                 #  sc  phenotype  sex
                                                                     "Heinonen, S. (2017)"="GSE92405eset",                                #  sc  phenotype  sex
                                                                     "Sharma, NK. (2016)"="GSE95674eset",                                 #  sc  phenotype  sex
                                                                     "GTEx microarray"="GTEx_microarray_sc",                              #  sc             sex
                                                                     "Lonsdale, J. (2013) sc"="GTEx_v4_sc",                               #  sc             sex
                                                                     "Aguet, F. (2017) sc"="GTEx_v6_sc",                                  #  sc             sex
                                                                     "GTEx v7 sc"="GTEx_v7_sc",                                           #  sc             sex
                                                                     "Aguet, F. (2020) sc"="GTEx_v8_sc",                                  #  sc             sex
                                                                     "Drong, A. (2013)"="E_MTAB_54",                                      #  sc             sex
                                                                     "Das, SK. (2015)"="GSE65221eset",                                    #  sc             sex  ethnicity
                                                                     "Naukkarinen, J. (2013)"="E_MTAB_1895",                              #  sc             sex
                                                                     "Nookaew, I. (2013)"="GSE27916eset",                                 #  sc             sex
                                                                     "MacLaren, RE. (2010) sc"="GSE15524eset_sc",                         #  sc  phenotype  sex
                                                                     "Defour, M. (2020)"="GSE154610eset",                                 #  sc  phenotype  sex
                                                                     "Johansson, LE. (2012) Baseline"="GSE35411eset_Baseline",            #  sc  phenotype  sex
                                                                     "Matualatupauw, JC. (2017)"="GSE87382eset",                          #  sc  phenotype  sex
                                                                     "Hardy, OT. (2011) sc"="GSE20950eset_sc",                            #  sc             sex
                                                                     "Du Plessis, J. (2015) sc"="GSE58979eset_sc",                        #  sc  phenotype  sex
                                                                     "Van Bussel, IPG. (2017)"="GSE84046eset_Baseline",                   #  sc  phenotype  sex
                                                                     "Bollepalli, S. (2018)"="GSE103766eset_Baseline",                    #  sc             sex
                                                                     "Rey, F. (2021)"="GSE166047eset"                                     #  sc             sex
                                                                     ),
                                       omental=c("Keller, M. (2017) om"="Keller_et_al_om",                                 #  om  phenotype  sex
                                                 "Krieg, L. (2021) om"="Krieg_et_al_om",                                   #  om  phenotype  sex
                                                 "Lonsdale, J. (2013) om"="GTEx_v4_om",                                    #  om             sex
                                                 "Aguet, F. (2017) om"="GTEx_v6_om",                                       #  om             sex
                                                 "GTEx v7 om"="GTEx_v7_om",                                                #  om             sex
                                                 "Aguet, F. (2020) om"="GTEx_v8_om",                                       #  om             sex
                                                 "MacLaren, RE. (2010) om"="GSE15524eset_om",                              #  om  phenotype  sex
                                                 "Hardy, OT. (2011) om"="GSE20950eset_om",                                 #  om             sex
                                                 "Du Plessis, J. (2015) om"="GSE58979eset_om"                              #  om  phenotype  sex
                                                 )
                                       ),
                                       # selected = "RIKENeset",
                                       multiple = TRUE,
                                       # width = "300px",
                                       options = pickerOptions(actionsBox = T, size = 10, liveSearch = T),
                                       choicesOpt = list(subtext = c(rep("sc",26),rep("om",9)))),
                           textAreaInput(inputId = "gene_sex",
                                         label = "Gene",
                                         value = "",
                                         # width = "300px",
                                         height = "200px",
                                         placeholder = "One per line, eg:\nPPARG\nADIPOQ"), # nolint
                           uiOutput("ui_statistics_sex"),
                           uiOutput("ui_effect_sex"),
                           uiOutput("ui_statistics_sex_protein"),
                           tags$label(style = "font-weight: bold; padding-right: 10px; margin: 0;", "Add proteomics plots", `for` = "proteomics_sex"),
                           materialSwitch(inputId = "proteomics_sex", label = NULL, status="primary"),
                           actionButton(inputId = "SearchButton_sex", label = "Submit")
                         ),
                         conditionalPanel(
                           condition = "input.navbar_selected == 'BMI_cross'",
                           br(),
                           pickerInput(inputId = "cohort_bmicross",
                                       label = "Cohort",
                                       choices = list(subcutaneous=c("Kerr, A. (2020) Baseline"="DEOSHeset_Baseline",                     #  sc  phenotype
                                                                     "Petrus, P. (2018) Baseline"="POeset_Baseline",                      #  sc  phenotype
                                                                     "Arner, P. (2018)"="SOWOTeset",                                      #  sc  phenotype
                                                                     "Keller, M. (2017) sc"="Keller_et_al_sc",                            #  sc  phenotype  sex
                                                                     "Arner, E. (2012)"="RIKENeset",                                      #  sc  phenotype
                                                                     "Stančáková, A. (2012)"="METSIM_204eset",                            #  sc  phenotype
                                                                     "Raulerson, C. (2019)"="METSIM_434eset",                             #  sc  phenotype
                                                                     "Civelek, M. (2017)"="METSIM_770eset",                               #  sc  phenotype
                                                                     "Krieg, L. (2021) sc"="Krieg_et_al_sc",                              #  sc  phenotype  sex
                                                                     "Arner, P. (2016) sc"="EMIFeset_sc",                                 #  sc  phenotype
                                                                     "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",         #  sc  phenotype  sex
                                                                     "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",        #  sc  phenotype  sex
                                                                     "Winnier, DA. (2015)"="GSE64567eset",                                #  sc  phenotype  sex
                                                                     "Nono Nankam, PA. (2020)"="Nankam_et_al",                            #  sc  phenotype
                                                                     "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",                 #  sc  phenotype  sex
                                                                     "Salcedo-Tacuma, D. (2022) sc"="GSE188799eset_sc",                   #  sc  phenotype
                                                                     "MacLaren, RE. (2010) sc"="GSE15524eset_sc",                         #  sc  phenotype  sex
                                                                     "Johansson, LE. (2012) Baseline"="GSE35411eset_Baseline",            #  sc  phenotype  sex
                                                                     "Matualatupauw, JC. (2017)"="GSE87382eset",                          #  sc  phenotype  sex
                                                                     "Du Plessis, J. (2015) sc"="GSE58979eset_sc"                         #  sc  phenotype  sex
                                                                     ),
                                                      omental=c("Keller, M. (2017) om"="Keller_et_al_om",                                 #  om  phenotype  sex
                                                                "Krieg, L. (2021) om"="Krieg_et_al_om",                                   #  om  phenotype  sex
                                                                "Arner, P. (2016) om"="EMIFeset_om",                                      #  om  phenotype
                                                                "Barberio, MD. (2019)"="GSE88837eset",                                    #  om  phenotype       ethnicity
                                                                "Salcedo-Tacuma, D. (2022) om"="GSE188799eset_om",                        #  om  phenotype
                                                                "MacLaren, RE. (2010) om"="GSE15524eset_om",                              #  om  phenotype  sex
                                                                "Du Plessis, J. (2015) om"="GSE58979eset_om",                             #  om  phenotype  sex
                                                                "Aguilera, CM. (2015)"="GSE9624eset"                                      #  om  phenotype
                                                                )
                                                      ),
                                       
                                       # selected = "RIKENeset",
                                       multiple = TRUE,
                                       # width = "300px",
                                       options = pickerOptions(actionsBox = T, size = 10, liveSearch = T),
                                       choicesOpt = list(subtext = c(rep("sc",20),rep("om",8)))),
                           textAreaInput(inputId = "gene_bmicross",
                                         label = "Gene",
                                         value = "",
                                         # width = "300px",
                                         height = "200px",
                                         placeholder = "One per line, eg:\nPPARG\nADIPOQ"), # nolint
                           uiOutput("ui_group_bmicross_RNA"),
                           uiOutput("ui_statistics_bmicross"),
                           uiOutput("ui_effect_bmicross"),
                           uiOutput("ui_group_bmicross_protein"),
                           uiOutput("ui_statistics_bmicross_protein"),
                           tags$label(style = "font-weight: bold; padding-right: 10px; margin: 0;", "Add proteomics plots", `for` = "proteomics_bmicross"),
                           materialSwitch(inputId = "proteomics_bmicross", label = NULL, status="primary"),
                           actionButton(inputId = "SearchButton_bmicross", label = "Submit")
                         ),
                         conditionalPanel(
                           condition = "input.navbar_selected == 'BMI_long'",
                           br(),
                           pickerInput(inputId = "cohort_bmilong",
                                       label = "Cohort",
                                       choices = list(`Bariatric surgery`=c("Kerr, A. (2020) weight loss"="DEOSHeset_WL",                        #  wl
                                                                            "Petrus, P. (2018) weight loss"="POeset_WL"                          #  wl
                                                                            ),
                                                      Diet=c("Imbert, A. (2022) weight loss"="GSE141221_diogenes1_WL",            #  wl
                                                             "Armenise, C. (2017) weight loss"="GSE95640_diogenes2_WL",           #  wl
                                                             "Vink, RG. (2017) weight loss"="GSE77962eset_WL",                    #  wl
                                                             "Johansson, LE. (2012) weight loss"="GSE35411eset_WL"                #  wl
                                                             )
                                                      ),
                                       
                                       # selected = "RIKENeset",
                                       multiple = TRUE,
                                       # width = "300px",
                                       options = pickerOptions(actionsBox = T, size = 10, liveSearch = T, maxOptions = 100),
                                       choicesOpt = list(subtext = c(rep("sc",5)))),
                           textAreaInput(inputId = "gene_bmilong",
                                         label = "Gene",
                                         value = "",
                                         # width = "300px",
                                         height = "200px",
                                         placeholder = "One per line, eg:\nPPARG\nADIPOQ"), # nolint
                           br(),
                           "OR/AND",
                           br(),
                           br(),
                           pickerInput(inputId = "trait_bmilong",
                                       label = "Trait",
                                       choices = trait_vector_bmilong,
                                       multiple = T,
                                       # width = "300px",
                                       options = pickerOptions(actionsBox = T, size = 10, maxOptions = 100,liveSearch = T)),
                           prettyRadioButtons(inputId = "statistics_bmilong",
                                              label = "Statistical test",
                                              selected = "wilcox.test",
                                              choices = c("t-Test" = "t.test", "Wilcoxon Test" = "wilcox.test"),
                                              inline = TRUE),
                           # actionBttn(inputId = "SearchButton_bmilong",
                           #            label = "Submit",
                           #            # style = "simple",
                           #            # color = "primary",
                           #            size = "sm"),
                           actionButton(inputId = "SearchButton_bmilong", label = "Submit")
                         )
                ),  # end tabPanel Input
                tabPanel(title = "Customization",
                         br(),
                         div("You can use default settings or customize visualizations below:"),
                         br(),
                         radioButtons(inputId="gradient_col",
                                      label = "Gradient colors",
                                      choices = c("Default" = "portalcol", 
                                                  "Purple to green" = "PRGn", 
                                                  "Red to Blue" = "RdBu", 
                                                  "Viridis" = "virid"),
                                      selected = "portalcol"),
                         radioButtons(inputId="discrete_col",
                                      label = "Discrete colors", 
                                      choiceNames =  c("Default", "Colorbrewer 2"), 
                                      choiceValues = c("default", "CB2"), 
                                      selected = "default"),
                         uiOutput("ui_discrete_col_CB2"),
                         div("See", a(href="https://colorbrewer2.org/", target = "_blank", "ColorBrewer2"), "for more"),
                         br(),
                         selectInput(inputId= "cluster_method", label="Clustering method for heatmaps", choices = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"), selected = "ward.D")
                )  # end tabPanel Customization
              )  # end tabsetPanel 
          )  # end div id="side-tab"
        )  # end sidebarPanel
    )
    ,
    mainPanel(
      width = 9,
      navbarPage(title = "Clinical association",
                 id = "navbar_selected",
                 tabPanel(title = "Phenotype correlation", 
                          value = "Phenotype correlation",
                          # div("Brief description for Phenotype correlation",class="detail-text"),
                          includeHTML(path = "htmls/description_pheno_corr.html"),
                          br(),
                          tabsetPanel(
                            tabPanel("Plots",
                                     br(),
                                     uiOutput(outputId = "ui_plot_pc_RNA"),
                                     br(),
                                     uiOutput(outputId = "ui_plot_pc_protein")),
                            tabPanel("Tables",
                                     br(),
                                     uiOutput("ui_df_pc_RNA"),
                                     br(),
                                     uiOutput("ui_df_pc_protein")),
                            tabPanel("Details",
                                     # br(),
                                     # div("Method or cohort details", class="detail-text"),
                                     includeHTML(path = "htmls/details_pheno_corr.html")),
                            tabPanel("Interpretation",
                                     # br(),
                                     # div("How to interpret results", class="detail-text"),
                                     includeHTML(path = "htmls/interpretation_pheno_corr.html"))
                          ),
                          tags$br(),
                          tags$hr(),
                          uiOutput(outputId = "ui_citation_pc")
                 ),
                 tabPanel(title = "Sex", 
                          value = "Sex",
                          # div("Brief description for Sex",class="detail-text"),
                          includeHTML(path = "htmls/description_sex.html"),
                          br(),
                          tabsetPanel(
                            tabPanel("Plots",
                                     br(),
                                     uiOutput(outputId = "ui_plot_sex_RNA"),
                                     br(),
                                     uiOutput(outputId = "ui_plot_sex_protein")), 
                            tabPanel("Tables",
                                     br(),
                                     uiOutput(outputId = "ui_df_sex_RNA"),
                                     br(),
                                     uiOutput(outputId = "ui_df_sex_protein")), 
                            tabPanel("Details",
                                     # br(),
                                     # div("Method or cohort details", class="detail-text"),
                                     includeHTML(path = "htmls/details_sex.html")),
                            tabPanel("Interpretation",
                                     # br(),
                                     # div("How to interpret results", class="detail-text"),
                                     includeHTML(path = "htmls/interpretation_sex.html"))
                          ),
                          tags$br(),
                          tags$hr(),
                          uiOutput(outputId = "ui_citation_sex")
                 ),
                 tabPanel(title = "Obese", 
                          value = "BMI_cross",
                          # div("Brief description for Cross-sectional BMI",class="detail-text"),
                          includeHTML(path = "htmls/description_obesity.html"),
                          br(),
                          tabsetPanel(
                            tabPanel("Plots",
                                     br(),
                                     uiOutput(outputId = "ui_plot_bmicross_RNA"),
                                     br(),
                                     uiOutput(outputId = "ui_plot_bmicross_protein")),
                            tabPanel("Tables",
                                     br(),
                                     uiOutput(outputId = "ui_df_bmicross_RNA"),
                                     br(),
                                     uiOutput(outputId = "ui_df_bmicross_protein")),
                            tabPanel("Details",
                                     # br(),
                                     # div("Method or cohort details", class="detail-text"),
                                     includeHTML(path = "htmls/details_obesity.html")),
                            tabPanel("Interpretation",
                                     # br(),
                                     # div("How to interpret results", class="detail-text"),
                                     includeHTML(path = "htmls/interpretation_obesity.html"))
                          ),
                          tags$br(),
                          tags$hr(),
                          uiOutput(outputId = "ui_citation_bmicross")
                 ),  # end tabPanel title = "Cross-sectional"
                 tabPanel(title = "Weight loss", 
                          value = "BMI_long",
                          # div("Brief description for Longitudinal BMI",class="detail-text"),
                          includeHTML(path = "htmls/description_weightloss.html"),
                          br(),
                          tabsetPanel(
                            tabPanel("Plots",
                                     br(),
                                     uiOutput(outputId = "ui_plot_bmilong")),
                            tabPanel("Tables",
                                     br(),
                                     uiOutput(outputId = "ui_df_bmilong")),
                            tabPanel("Details",
                                     # br(),
                                     # div("Method or cohort details", class="detail-text"),
                                     includeHTML(path = "htmls/details_weightloss.html")),
                            tabPanel("Interpretation",
                                     # br(),
                                     # div("How to interpret results", class="detail-text"),
                                     includeHTML(path = "htmls/interpretation_weightloss.html"))
                          ),
                          tags$br(),
                          tags$hr(),
                          uiOutput(outputId = "ui_citation_bmilong")
                 )  # end tabPanel title = "Longitudinal"
                 # navbarMenu(title = "BMI",
                 #            tabPanel(title = "Cross-sectional", 
                 #                     value = "BMI_cross",
                 #                     div("Brief description for Cross-sectional BMI",class="detail-text"),
                 #                     br(),
                 #                     tabsetPanel(
                 #                       tabPanel("Plots",
                 #                                br(),
                 #                                uiOutput(outputId = "ui_plot_bmicross")),
                 #                       tabPanel("Tables",
                 #                                br(),
                 #                                uiOutput(outputId = "ui_df_bmicross")),
                 #                       tabPanel("Details",
                 #                                br(),
                 #                                div("Method or cohort details", class="detail-text")),
                 #                       tabPanel("Interpretation",
                 #                                br(),
                 #                                div("How to interpret results", class="detail-text"))
                 #                     ),
                 #                     tags$br(),
                 #                     tags$hr(),
                 #                     uiOutput(outputId = "ui_citation_bmicross")
                 #            ),  # end tabPanel title = "Cross-sectional"
                 #            tabPanel(title = "Longitudinal", 
                 #                     value = "BMI_long",
                 #                     div("Brief description for Longitudinal BMI",class="detail-text"),
                 #                     br(),
                 #                     tabsetPanel(
                 #                       tabPanel("Plots",
                 #                                br(),
                 #                                uiOutput(outputId = "ui_plot_bmilong")),
                 #                       tabPanel("Tables",
                 #                                br(),
                 #                                uiOutput(outputId = "ui_df_bmilong")),
                 #                       tabPanel("Details",
                 #                                br(),
                 #                                div("Method or cohort details", class="detail-text")),
                 #                       tabPanel("Interpretation",
                 #                                br(),
                 #                                div("How to interpret results", class="detail-text"))
                 #                     ),
                 #                     tags$br(),
                 #                     tags$hr(),
                 #                     uiOutput(outputId = "ui_citation_bmilong")
                 #            )  # end tabPanel title = "Longitudinal"
                 # )  # end navbarMenu title = "BMI"
      )  # end navbarPage
    )  # end mainPanel
  )  # end sidebarLayout
))
