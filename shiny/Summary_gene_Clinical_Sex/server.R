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
      if (file.exists(paste0('./data/',input$data_from_html,'.RDS'))) {
        
        output$plot <- renderPlot({
          temp <- readRDS(paste0('./data/',input$data_from_html,'.RDS'))

          if (is.na(temp$pval.Q)) {
            meta::forest(temp, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA", fontfamily = "Red Hat Display" ,col.square = "#1a4659", col.diamond.common = "#e2c744", col.diamond.random = "#e2c744", col.square.lines = "#000000")
          } else {
            if (temp$pval.Q<0.05) {
              if (temp$TE.random>0) {
                meta::forest(temp, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA", fontfamily = "Red Hat Display" ,col.square = "#1a4659", col.diamond.common = "#e2c744", col.diamond.random = "#e2c744", col.square.lines = "#000000")
              } else {
                meta::forest(temp, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA", fontfamily = "Red Hat Display" ,col.square = "#1a4659", col.diamond.common = "#e2c744", col.diamond.random = "#e2c744", col.square.lines = "#000000")
              }
            } else {
              if (temp$TE.common>0) {
                meta::forest(temp, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA", fontfamily = "Red Hat Display" ,col.square = "#1a4659", col.diamond.common = "#e2c744", col.diamond.random = "#e2c744", col.square.lines = "#000000")
              } else {
                meta::forest(temp, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA", fontfamily = "Red Hat Display" ,col.square = "#1a4659", col.diamond.common = "#e2c744", col.diamond.random = "#e2c744", col.square.lines = "#000000")
              }
            }
          }
          
        })
        
        plotOutput("plot", height = "400px")
        
      } else {
        h1(paste0(input$data_from_html,' is not in this module'),align = "center")
      }
    }
  )
})


