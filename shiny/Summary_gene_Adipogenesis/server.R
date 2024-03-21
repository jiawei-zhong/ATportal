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
  
  output$ui <- renderUI(
    if (is.null(input$data_from_html)) {
      h1("undefined",align = "center")
    } else {
      if (input$data_from_html %in% rownames(adipogenesis)) {
        
        output$plot <- renderPlot({
          # 
          # df <- data.frame(expression=exprs(adipogenesis)[input$data_from_html,],
          #                  timepoint=adipogenesis$timepoint,
          #                  log_timepoint=adipogenesis$log_timepoint)
          # df <- df[!df$timepoint==20160,]

          # ggline(adipogenesis[adipogenesis$gene==input$data_from_html,,drop=F], color = "#1b9e77", x = "log_timepoint", y = "normalized", numeric.x.axis = TRUE, xlab = "Log(timepoint)", ylab = paste0(input$data_from_html, " expression")) + 
          #   theme(text = element_text(family = "Arial"))

          df <- data.frame(pData(adipogenesis),expression=t(exprs(adipogenesis)[input$data_from_html,,drop=F]))

          colnames(df)[4] <- "expression"

          df <- df[!df$timepoint==20160,]

          df$log_timepoint <- log1p(df$timepoint)

          # df$expression <- as.vector(scale(df$expression))

          ggline(df, x = "log_timepoint", y = "expression",numeric.x.axis = T, add = c("mean_se"), color = "#1a4659", xlab = "Log(timepoint)", ylab = paste0(input$data_from_html, " expression")) +
            theme(text = element_text(family = "Red Hat Display"))
          

        })
        
        plotOutput("plot", height = "250px")
        
      } else {
        h1(paste0(input$data_from_html,' is not in this module'),align = "center")
      }
    }
  )
})

