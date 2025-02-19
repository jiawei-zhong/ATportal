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
      if (gene_reactive() %in% gene_list) {
        
        output$plot1 <- renderPlot({


          temp <- exprs(adipogenesis)[gene_reactive(),,drop=F]
          df <- data.frame(pd,expression=t(temp))

          colnames(df)[4] <- "expression"

          df <- df[!df$timepoint==20160,]

          df$log_timepoint <- log1p(df$timepoint)

          # df$expression <- as.vector(scale(df$expression))

          ggline(df, x = "log_timepoint", y = "expression",numeric.x.axis = T, add = c("mean_se"), color = "#1a4659", xlab = "Log(timepoint)", ylab = paste0(gene_reactive(), " expression")) +
            theme(text = element_text(family = "Red Hat Display"))
          
          

        })

        gene_protein <- intersect(gene_reactive(), data_prot_raw$gene)

        if (length(gene_protein) > 0) {
          
          output$plot2 <- renderPlot({

            df <- data_prot_raw[data_prot_raw$gene %in% gene_protein,] %>% na.omit()

            ggline(data = df, x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "#1a4659", ylab = paste0(gene_reactive(), " expression")) +
              theme(text = element_text(family = "Red Hat Display"))
          })

          fluidRow(
            column(width = 6, 
                   style = 'padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                   "Transcriptomics",
                   plotOutput("plot1", height = "230px")
            ),
            column(width = 6, 
                   style='padding-left:15px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                   "Proteomics",
                   plotOutput("plot2", height = "230px")
            )
          )
        } else {
          fluidRow(
            column(width = 12,
                   style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                   plotOutput("plot1", height = "250px")
            )
          )
        }
        
      } else {
        h1(paste0(gene_reactive(),' is not in this module'),align = "center")
      }
    }
  )
})

