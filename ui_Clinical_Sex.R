tabPanel(title = "Sex difference",
         titlePanel("Boxplot for different genders"),
         sidebarLayout(position = "left",
                       sidebarPanel = sidebarPanel(helpText("Create box plot about sex difference"), # nolint
                                                   textAreaInput(inputId = "gene_sex", 
                                                                 label = "Input gene name", 
                                                                 value = "LEP", 
                                                                 width = "300px", 
                                                                 height = "50px", 
                                                                 placeholder = "Only 1 gene"),
                                                   pickerInput(inputId = "cohort_sex", 
                                                               label = "Choose a cohort", 
                                                               choices = list(subcutaneous=c("DiOGenes1" = "GSE141221_diogenes1", "DiOGenes2" = "GSE95640_diogenes2", "GTEx sc" = "GTEx_sc","Keller et al. sc"="Keller_et_al_sc"),
                                                                              omental=c("GTEx om" = "GTEx_om","Keller et al. om"="Keller_et_al_om")
                                                                              # epiploic=c("Krieg et al."="Krieg_et_al_ep"),
                                                                              # `mesenteric duodenum`=c("Krieg et al."="Krieg_et_al_md")
                                                               ),
                                                               width = "300px",
                                                               options = list(size = 10, `live-search` = TRUE),
                                                               choicesOpt = list(subtext = c("sc","sc","sc","sc","om","om"))),
                                                   prettyRadioButtons(inputId = "statistical_test_sex", 
                                                                      label = "Statistical test", 
                                                                      selected = "t.test", 
                                                                      choices = c("t-Test" = "t.test", "Wilcoxon Test" = "wilcox.test"), 
                                                                      inline = TRUE),
                                                   actionBttn(inputId = "SearchButton_sex", 
                                                              label = "Search", 
                                                              style = "simple",
                                                              color = "primary",
                                                              size = "sm"),
                                                   width = 3),
                       mainPanel = mainPanel(fluidRow(column(12, uiOutput("ui_boxplot_sex"))))
         )
)