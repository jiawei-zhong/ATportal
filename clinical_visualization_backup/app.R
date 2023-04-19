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
  
  navbarPage(title = "WAT portal",
             tabPanel(title = "Home", "This part is under progress"),
             navbarMenu(title = "Meta analysis",
                        tabPanel(title = "Heatmap",
                                 titlePanel("Heatmap for gene expression and clinical measures"),
                                 sidebarLayout(position = "left", 
                                               sidebarPanel = sidebarPanel(helpText("Create heatmap with correlation between gene expression and clinical measures."), # nolint
                                                                           pickerInput(inputId = "cohort_hm",
                                                                                       label = "Choose a cohort",
                                                                                       choices = list(subcutaneous=c("DEOSH"="DEOSHeset", "DiOGenes1"="GSE141221_diogenes1", "DiOGenes2"="GSE95640_diogenes2", "EMIF"="EMIFeset_sc", "Krieg et al."="Krieg_et_al_sc", "PO"="POeset_Baseline", "RIKEN"="RIKENeset", "SOWOT"="SOWOTeset"),
                                                                                                      omental=c("EMIF"="EMIFeset_om", "Krieg et al."="Krieg_et_al_om")
                                                                                                      # epiploic=c("Krieg et al."="Krieg_et_al_ep"),
                                                                                                      # `mesenteric duodenum`=c("Krieg et al."="Krieg_et_al_md")
                                                                                       ), 
                                                                                       selected = "RIKENeset",
                                                                                       width = "300px",
                                                                                       options = list(size = 10, `live-search` = TRUE),
                                                                                       choicesOpt = list(subtext = c("sc","sc","sc","sc","sc","sc","sc","sc","om","om"))),
                                                                           textAreaInput(inputId = "gene_hm", 
                                                                                         label = "Input list of genes", 
                                                                                         value = "", 
                                                                                         width = "300px", 
                                                                                         height = "300px", 
                                                                                         placeholder = "A least 2 genes. Input empty will use the first 30 genes of the cohort"), # nolint
                                                                           prettyCheckboxGroup(inputId = "clinical_measure_hm", 
                                                                                               label = "Clinical measures", 
                                                                                               choices = clinical_measure_vector,
                                                                                               # selected = c("BMI","HOMA"),
                                                                                               width = "300px",
                                                                                               inline = T),
                                                                           materialSwitch(inputId = "select_all_hm", label = "Select all", right = TRUE, value = F),
                                                                           prettyRadioButtons(inputId = "correlation_method_hm", 
                                                                                              label = "Correlation method", 
                                                                                              selected = "pearson", 
                                                                                              choices = c("Pearson" = "pearson", "Spearman" = "spearman"), 
                                                                                              inline = TRUE),
                                                                           actionBttn(inputId = "SearchButton_hm", 
                                                                                      label = "Search", 
                                                                                      style = "simple",
                                                                                      color = "primary",
                                                                                      size = "sm"),
                                                                           width = 3),
                                               mainPanel = mainPanel(fluidRow(column(12, uiOutput("ui.heatmap"))),
                                                                     br(),
                                                                     fluidRow(column(12, uiOutput("ui.cor_table")))
                                               )
                                 )
                        ),
                        tabPanel(title = "Forest plot",
                                 titlePanel("Forest plot for gene expression and clinical measure"),
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
                        
             ),
             tabPanel(title = "Module 2", "Module 2 is under progress"),
             tabPanel(title = "Module 3", "Module 3 is under progress"),
             tabPanel(title = "Module 4", "Module 4 is under progress"),
             tabPanel(title = "Module 5", "Module 5 is under progress"),
             tabPanel(title = "About", "This part is under progress")
  )
  
)





# Define server logic ----
server <- function(input, output, session) {
  
  
  ## heatmap
  observe({updatePrettyCheckboxGroup(session = session, 
                                     inputId = 'clinical_measure_hm',  
                                     choices = clinical_measure_vector,
                                     selected = if (input$select_all_hm) {clinical_measure_vector},
                                     inline = T )})
  
  
  
  plotdata_hm <- eventReactive(input$SearchButton_hm, {
    eset <- readRDS(paste0('../clinical_eset/',input$cohort_hm,'.rds'))
    exprs <- exprs(eset)
    if (nchar(input$gene_hm) > 0) {
      exprs <- exprs[rownames(exprs) %in% unlist(strsplit(input$gene_hm,'\n')),]
    } else {
      exprs <- exprs[1:30,]
    }
    pData <- pData(eset)
    pData <- pData[,colnames(pData) %in% input$clinical_measure_hm]
    pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
    
    
    
    cor_df <- t(rbind(exprs,t(pData)))
    cor_df <- rcorr(cor_df,type=input$correlation_method_hm)
    cor_df$r <- (cor_df$r)[rownames(exprs),colnames(pData)]
    cor_df$P <- (cor_df$P)[rownames(exprs),colnames(pData)]
    cor_df$P <- (-log10(cor_df$P))
    # cor_df$n <- (cor_df$n)[rownames(exprs),colnames(pData)]
    
    row_order <- hclust(d = dist(cor_df$r))$order
    column_order <- hclust(d = dist(t(cor_df$r)))$order
    # cor_df$order_row <- row_order
    # cor_df$order_col <- column_order
    cor_df$r <- cor_df$r[row_order,column_order]
    cor_df$P <- cor_df$P[row_order,column_order]
    
    p1 <- pheatmap::pheatmap(mat = cor_df$r,
                             cluster_cols = T,
                             cluster_rows = T,
                             cellwidth = 20,
                             cellheight = 20,
                             border_color = NA,
                             treeheight_col = 0,
                             treeheight_row = 30, 
                             fontsize = 15,
                             main = "Heatmap of correlation",
                             silent = T,
                             show_colnames = T,
                             color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
    )
    
    p2 <- pheatmap::pheatmap(mat = cor_df$P,
                             cluster_cols = F,
                             cluster_rows = F,
                             display_numbers = ifelse(test = cor_df$P > (-log10(0.001)),
                                                      yes =  "***",
                                                      no =  ifelse(test = cor_df$P > (-log10(0.01)),
                                                                   yes =  "**",
                                                                   no =  ifelse(test = cor_df$P > (-log10(0.05)),
                                                                                yes =  "*",
                                                                                no =  ""))),
                             cellwidth = 20,
                             cellheight = 20,
                             border_color = NA,
                             treeheight_col = 0,
                             treeheight_row = 0,
                             fontsize = 15,
                             main = "Heatmap of P-value (-log10)",
                             silent = T,
                             show_colnames = T
    )
    

    cor_df$heatmap <- plot_grid(as.ggplot(p1), as.ggplot(p2),ncol=2)
    
    cor_df
  })
  
  output$heat_map <- renderPlot({
    
    plotdata_hm()$heatmap
    
  })
  
  
  
  
  output$r <- DT::renderDataTable(server = FALSE,{return(plotdata_hm()$r)},
                                  extensions = c('Buttons'),
                                  options = list(scrollX = TRUE,
                                                 pageLength = 10,
                                                 lengthMenu = c(10, 25, 50, 100),
                                                 dom = 'Blfrtip',
                                                 buttons = list(
                                                   list(extend = "csv", text = "Download Current Page", filename = "page",
                                                        exportOptions = list(modifier = list(page = "current"))),
                                                   list(extend = "csv", text = "Download Full matrix", filename = "data",
                                                        exportOptions = list(modifier = list(page = "all"))))))
  output$P <- DT::renderDataTable(server = FALSE,{return(plotdata_hm()$P)},
                                  extensions = c('Buttons'),
                                  options = list(scrollX = TRUE,
                                                 pageLength = 10,
                                                 lengthMenu = c(10, 25, 50, 100),
                                                 dom = 'Blfrtip',
                                                 buttons = list(
                                                   list(extend = "csv", text = "Download Current Page", filename = "page",
                                                        exportOptions = list(modifier = list(page = "current"))),
                                                   list(extend = "csv", text = "Download Full matrix", filename = "data",
                                                        exportOptions = list(modifier = list(page = "all"))))))
  
  # output$P <- DT::renderDataTable({return(plotdata()$P)},extensions = c('Buttons'),options = list(scrollX = TRUE))
  # output$download_r <- downloadHandler(filename = function(){"matrix_r.rho.csv"}, 
  #                                      content = function(fname){write.csv(plotdata()$r, fname, quote = F, row.names = T)})
  # output$download_P <- downloadHandler(filename = function(){"matrix_P.csv"}, 
  #                                      content = function(fname){write.csv(plotdata()$P, fname, quote = F, row.names = T)})
  
  output$downloadPlot_hm <- downloadHandler(filename = function() {"heatmap.pdf"},
                                            content = function(fname) {ggsave(filename = fname, plot = plotdata_hm()$heatmap, height = 40, width = 20, units = "in")})
  
  
  
  output$ui.heatmap <- renderUI({ if (input$SearchButton_hm) {
    list(box(width=12,
             style=paste0('height: ',min(30*20+120,nrow(plotdata_hm()$r)*20+120),'px; overflow-x: scroll; overflow-y: scroll;'),
             plotOutput(outputId = "heat_map", height = paste0(nrow(plotdata_hm()$r)*20+100,'px'), width = paste0(ncol(plotdata_hm()$r)*20+900,'px'))),
         column(12,downloadBttn(outputId = 'downloadPlot_hm', label = "Download heatmap", style= "simple", color = "primary", size = "sm"))
    )}
  })
  
  
  
  
  output$ui.cor_table <- renderUI({if (input$SearchButton_hm) {
    tabsetPanel(id = "inTabset",
                tabPanel(title = "r/rho matrix",
                         column(12,DT::dataTableOutput("r"))),
                tabPanel(title = "P-value matrix (-log10)", 
                         column(12,DT::dataTableOutput("P"))))
  }})
  
  # output$ui.table <- renderUI({if (input$SearchButton) {
  #   tabsetPanel(id = "inTabset",
  #               tabPanel(title = "r/rho matrix",
  #                        column(12,DT::dataTableOutput("r")),
  #                        column(12,downloadButton('download_r',"Download r/rho matrix"))),
  #               tabPanel(title = "P-value matrix (-log10)", 
  #                        column(12,DT::dataTableOutput("P")),
  #                        column(12,downloadButton('download_P',"Download P-value matrix"))))
  # }})
  
  
  
  
  
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