tabPanel(title = "Forest plot",
         titlePanel("Forest plot for gene expression and clinical measure"),
         sidebarLayout(position = "left", 
                       sidebarPanel = sidebarPanel(helpText("Create forest plot with correlation between gene expression and clinical measure."), # nolint
                                                   pickerInput(inputId = "cohort_fp",
                                                               label = "Choose cohorts",
                                                               choices = list(subcutaneous=c("DEOSH"="DEOSHeset_Baseline", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF"="EMIFeset_sc", "Krieg et al."="Krieg_et_al_sc", "Keller et al. sc"="Keller_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset", "METSIM 770"="METSIM_770eset", "METSIM 434"="METSIM_434eset", "METSIM 200"="METSIM_200eset"),
                                                                              omental=c("EMIF"="EMIFeset_om", "Krieg et al."="Krieg_et_al_om", "Keller et al."="Keller_et_al_om")
                                                                              # epiploic=c("Krieg et al."="Krieg_et_al_ep"),
                                                                              # `mesenteric duodenum`=c("Krieg et al."="Krieg_et_al_md")
                                                               ), 
                                                               selected = "RIKENeset",
                                                               width = "300px",
                                                               options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                                                               multiple = TRUE,
                                                               choicesOpt = list(subtext = c("sc","sc","sc","sc","sc","sc","sc","sc","sc","sc","sc","sc","om","om","om"))),
                                                   textAreaInput(inputId = "gene_fp", 
                                                                 label = "Input gene name", 
                                                                 value = "LEP", 
                                                                 width = "300px", 
                                                                 height = "50px", 
                                                                 placeholder = "Only 1 gene"), # nolint
                                                   pickerInput(inputId = "clinical_measure_fp", 
                                                               label = "Clinical measures", 
                                                               choices = clinical_measure_vector,
                                                               width = "300px"),
                                                   pickerInput(inputId = "clinical_measure_adjust_fp",
                                                               label = "Adjustment",
                                                               choices = NA,
                                                               width = "300px",
                                                               options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                                                               multiple = TRUE),
                                                   prettyRadioButtons(inputId = "correlation_method_fp", 
                                                                      label = "Correlation method", 
                                                                      selected = "pearson", 
                                                                      choices = c("Pearson" = "pearson", "Spearman" = "spearman"), 
                                                                      inline = TRUE),
                                                   actionBttn(inputId = "SearchButton_fp", 
                                                              label = "Search", 
                                                              style = "simple",
                                                              color = "primary",
                                                              size = "sm"),
                                                   
                                                   
                                                   width = 3),
                       mainPanel = mainPanel(fluidRow(column(12, uiOutput("ui_forestplot"))))
         )
)
