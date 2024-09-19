# Load packages ----

suppressMessages(library(shiny))
suppressMessages(library(shinythemes))
suppressMessages(library(shinyhelper))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))  # for box()
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
suppressMessages(library(extrafont))


shinyServer(function(input, output, session) {
  
  gene_reactive <- reactiveVal("")
  
  observe({
    if (!is.null(input$data_from_html) && substr(input$data_from_html,1,5)=="gene:") {
      gene_reactive(gsub(pattern = "gene:",replacement = "",x = input$data_from_html))
    }
  })
  
  output$ui <- renderUI(
    if (gene_reactive()=="") {
      h1("undefined",align = "center")
    } else {
      if (gene_reactive() %in% perturbation$GENE) {
        
        output$plot <- renderPlot({
          
          perturbation_sub <- perturbation[perturbation$GENE==gene_reactive(),]
          
          
          
          
          color_mapping <- setNames(perturbation_sub$color, perturbation_sub$perturbation)
          
          labels <- ifelse(perturbation_sub$color != "darkgrey", as.character(perturbation_sub$perturbation), NA)
          
          
          ggscatter(data = perturbation_sub,x = "logFC",y = "log10.adj.P.Val",size = 3, color = "color") +
            # facet_grid(cols = vars(group)) +
            scale_color_manual(values = c(`#1C4759` = "#1C4759", "darkgrey" = "darkgrey", `#E2C744` = "#E2C744")) +  # 自定义颜色映射
            
            geom_text_repel(aes(label = labels), 
                            box.padding = 0.35, point.padding = 0.3, 
                            max.overlaps = 10, size = 5, na.rm = TRUE) +
            theme(legend.position = "none",aspect.ratio = 0.5) +
            scale_x_continuous(name = expression(Log[2]~fold~change)) +
            scale_y_continuous(name = expression(-Log[10]~adjusted~P~value))
          
        })
        
        fluidRow(
          column(width = 12,
                 style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                 plotOutput("plot", height = "400px")
          )
        )
        
      } else {
        h1(paste0(gene_reactive(),' is not in this module'),align = "center")
      }
    }
  )
})


