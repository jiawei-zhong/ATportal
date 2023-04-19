# Load packages ----
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Hmisc))  
suppressMessages(library(grid))
suppressMessages(library(cowplot))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))  # for box()
suppressMessages(library(ggplotify))



clinical_measure_vector <- c("BMI" = "BMI", "HOMA" = "HOMA", "Age" = "Age", "WHR" = "WHR", "Waist" = "Waist", "Hip" = "Hip", "Glucose" = "Glucose", "Insulin" = "Insulin", "LEP" = "LEP", "TG" = "TG", "Chol" = "Chol", 
                             "HDL" = "HDL", "LDL" = "LDL", "CRP" = "CRP", "Hba1c" = "Hba1c", "TNF" = "TNF", "MCP1" = "MCP1", "cell_vol" = "cell_vol", "basal_TG" = "basal_TG", "iso_TG" = "iso_TG", "iso_basal" = "iso_basal")

cohort_forest <- c("DEOSH"="DEOSHeset", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF sc"="EMIFeset_sc", "Krieg et al. sc"="Krieg_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset","EMIF om"="EMIFeset_om", "Krieg et al. om"="Krieg_et_al_om")

# Define UI ----
ui <- fluidPage(
  # We MUST load the ECharts javascript library in advance
  # loadEChartsLibrary(),
  theme = shinythemes::shinytheme("flatly"),
  
  sidebarLayout(position = "left", 
                sidebarPanel = sidebarPanel(helpText("Create forest plot with correlation between gene expression and clinical measure."), # nolint
                                            pickerInput(inputId = "cohort_fp",
                                                        label = "Choose cohorts",
                                                        choices = list(subcutaneous=c("DEOSH"="DEOSHeset", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF"="EMIFeset_sc", "Krieg et al."="Krieg_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset"),
                                                                       omental=c("EMIF"="EMIFeset_om", "Krieg et al."="Krieg_et_al_om")
                                                                       # epiploic=c("Krieg et al."="Krieg_et_al_ep"),
                                                                       # `mesenteric duodenum`=c("Krieg et al."="Krieg_et_al_md")
                                                        ), 
                                                        selected = "RIKENeset",
                                                        width = "300px",
                                                        options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                                                        multiple = TRUE,
                                                        choicesOpt = list(subtext = c("sc","sc","sc","sc","sc","sc","sc","sc","om","om"))),
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
                mainPanel = mainPanel(fluidRow(column(12, uiOutput("ui.forestplot")))
                                      
                )
  )
)





# Define server logic ----
server <- function(input, output, session) {
  
  
  ## forest plot
  
  plotdata_fp <- eventReactive(input$SearchButton_fp, {
    
    
    cor <- list()
    for (i in 1:length(input$cohort_fp)) {
      eset <- readRDS(paste0('../clinical_eset/',input$cohort_fp[i],'.rds'))
      if (input$gene_fp %in% rownames(exprs(eset))) {
        if (all(is.na(pData(eset)[,input$clinical_measure_fp]))) {
          
        } else {
          # exprs(eset)[input$gene_fp,]
          # pData(eset)[,input$clinical_measure_fp]
          temp <- Hmisc::rcorr(exprs(eset)[input$gene_fp,],pData(eset)[,input$clinical_measure_fp],type = input$correlation_method_fp)
          cor$r <- c(cor$r,temp$r[1,2])
          cor$n <- c(cor$n,temp$n[1,2])
          cor$ID <- c(cor$ID,names(cohort_forest[cohort_forest %in% input$cohort_fp[i]]))
        }
      }
      
    }
    
    test <- meta::metacor(cor = cor$r,n = cor$n, studlab = cor$ID)
    test
  })
  
  output$forestplot <- renderPlot({
    
    meta::forest(plotdata_fp())
    
  })
  
  output$downloadPlot_fp <- downloadHandler(filename = function() {"forestplot.pdf"},
                                            content = function(fname) {
                                              pdf(fname, width=10, height=6)
                                              meta::forest(plotdata_fp())
                                              dev.off()
                                              })
  
  
  
  output$ui.forestplot <- renderUI({ if (input$SearchButton_fp) {
    
    list(column(12,plotOutput(outputId = "forestplot")),
         column(12,downloadBttn(outputId = 'downloadPlot_fp', label = "Download forest plot", style= "simple", color = "primary", size = "sm")))
    

    }
  })
  
}

# Run the app ----
shinyApp(ui = ui, server = server)