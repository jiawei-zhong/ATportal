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
      if (input$data_from_html %in% rownames(FANTOM)) {
        
        output$plot <- renderPlot({

          df <- data.frame(pData(FANTOM),expression=t(exprs(FANTOM)[input$data_from_html,,drop=F]))
          colnames(df)[4] <- "expression"

          df$classfication[df$ID %in% c("adipose, donor1", "adipose, donor2", "adipose, donor3", "adipose, donor4")] <- "other tissue"
          df$classfication[df$classfication %in% c("muture adipocyte", "NiGa adipocyte")] <- "adipocyte"


          df <- df[order(df$expression,decreasing = T),]

          df$rank <- 1:nrow(df)

          p1 <- ggbarplot(df, x = "rank", y = "expression", fill = "#adbec4", color = "#adbec4", ylab = paste0(input$data_from_html, " expression")) +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  axis.line.x = element_line(),  # 启用x轴线
                  axis.ticks.x = element_line(),  # 启用x轴刻度
                  axis.text.x = element_text())   # 启用x轴文本


          group_colors <- c("adipocyte" = "#1a4659", "adipose" = "#e2c744", "other tissue" = "grey")
          annotation_bar <- data.frame(
            xmin = 0:203 + 0.5,
            xmax = 1:204 + 0.5,
            ymin = -Inf,
            ymax = -max(df$expression) * 0.02,
            group = df$classfication
          )

          # 在条形图下方添加注释条
          p1 <- p1 + geom_rect(data = annotation_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), inherit.aes = FALSE) +
            scale_fill_manual(values = group_colors) +
            coord_cartesian(ylim = c(-max(df$expression) * 0.05, NA))

          # 添加x轴标题和标签
          p1 <- p1 + 
            scale_x_continuous(breaks = seq(0,200,50)+0.5, labels = seq(0,200,50), name = "Rank", expand = c(0.02, 0)) +
            theme(legend.position = c(0.9,0.8),
            legend.title = element_blank())

          p1

        })
        
        plotOutput("plot", height = "250px")
        
      } else {
        h1(paste0(input$data_from_html,' is not in this module'),align = "center")
      }
    }
  )
})


