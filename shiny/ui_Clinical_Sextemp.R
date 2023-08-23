tabPanel(title = "Sex difference",
         titlePanel("Forest plot for different genders"),
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
                                                               choices = list(subcutaneous=c("AAGMEx"="GSE95674eset",                               #  sc             sex
                                                                                             "ADIPOSTRESS Baseline"="GSE77962eset_Baseline",        #  sc  phenotype  sex
                                                                                             "Das et al."="GSE65221eset",                           #  sc             sex  ethnicity
                                                                                             "DiOGenes1"="GSE141221_diogenes1",                     #  sc  phenotype  sex
                                                                                             "DiOGenes2"="GSE95640_diogenes2",                      #  sc  phenotype  sex
                                                                                             "GTEx microarray"="GTEx_microarray_sc",                #  sc             sex
                                                                                             "GTEx v4 sc"="GTEx_v4_sc",                             #  sc             sex
                                                                                             "GTEx v6 sc"="GTEx_v6_sc",                             #  sc             sex
                                                                                             "GTEx v7 sc"="GTEx_v7_sc",                             #  sc             sex
                                                                                             "GTEx v8 sc"="GTEx_v8_sc",                             #  sc             sex
                                                                                             "Keller et al. 2017 sc"="Keller_et_al_sc",             #  sc  phenotype  sex
                                                                                             "Krieg et al. sc"="Krieg_et_al_sc",                    #  sc  phenotype  sex
                                                                                             "MolOBB"="E_MTAB_54",                                  #  sc             sex
                                                                                             "Naukkarinen et al."="E_MTAB_1895",                    #  sc             sex
                                                                                             "SOS Sib-Pair"="GSE27916eset",                         #  sc             sex
                                                                                             "VAGES"="GSE64567eset"                                 #  sc  phenotype  sex
                                                                                             ),
                                                                              omental=c("GTEx v4 om"="GTEx_v4_om",                             #  om             sex
                                                                                        "GTEx v6 om"="GTEx_v6_om",                             #  om             sex
                                                                                        "GTEx v7 om"="GTEx_v7_om",                             #  om             sex
                                                                                        "GTEx v8 om"="GTEx_v8_om",                             #  om             sex
                                                                                        "Keller et al. 2017 om"="Keller_et_al_om",             #  om  phenotype  sex
                                                                                        "Krieg et al. om"="Krieg_et_al_om"                     #  om  phenotype  sex
                                                                                        )
                                                                              # epiploic=c("Krieg et al."="Krieg_et_al_ep"),
                                                                              # `mesenteric duodenum`=c("Krieg et al."="Krieg_et_al_md")
                                                               ),
                                                               selected = "RIKENeset",
                                                               multiple = TRUE,
                                                               width = "300px",
                                                               options = list(size = 10, `live-search` = TRUE,`actions-box` = TRUE),
                                                               choicesOpt = list(subtext = c(rep("sc",16),rep("om",6)))),
                                                   actionBttn(inputId = "SearchButton_sex", 
                                                              label = "Search", 
                                                              style = "simple",
                                                              color = "primary",
                                                              size = "sm"),
                                                   width = 3),
                       mainPanel = mainPanel(fluidRow(column(12, uiOutput("ui_forestplot_sex"))))
         )
)
