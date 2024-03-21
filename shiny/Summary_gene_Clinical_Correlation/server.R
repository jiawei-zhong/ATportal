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
      if (input$data_from_html %in% rownames(meta_r)) {
        output$plot <- renderPlot({
          
          
          cor_df <- list()
          cor_df$r <- meta_r[input$data_from_html,,drop=F]
          cor_df$p <- meta_p[input$data_from_html,,drop=F]
          cor_df$p <- (-log10(cor_df$p))
          
          # row_order <- hclust(d = dist(cor_df$r))$order
          # column_order <- hclust(d = dist(t(cor_df$r)))$order
          # cor_df$r <- cor_df$r[row_order,column_order]
          # cor_df$p <- cor_df$p[row_order,column_order]
          # cor_df$r <- cor_df$r[,column_order,drop=F]
          # cor_df$p <- cor_df$p[,column_order,drop=F]
          
          temp <- data.frame(trait=colnames(cor_df$r),
                             trait_full=names(trait_vector)[match(colnames(cor_df$r), trait_vector)],
                             r=cor_df$r[1,],
                             p=cor_df$p[1,],
                             p_stars=ifelse(cor_df$p[1,] < -log10(0.05), "ns", ifelse(cor_df$p[1,] < -log10(0.01), "*", ifelse(cor_df$p[1,] < -log10(0.001), "**", ifelse(cor_df$p[1,] < -log10(0.0001), "***", "****")))))
          
          temp <- temp[order(temp$r,decreasing = T),] %>% na.omit()
          
          
          
          
          out_plot <- ggbarplot(temp, x = "trait_full", y = "r", fill = "p", color = NA) +
            geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
            geom_text(aes(label = p_stars, y = ifelse(r > 0, r + (range(temp$r)[2]-range(temp$r)[1])*0.01, r - (range(temp$r)[2]-range(temp$r)[1])*0.01)), 
                      vjust = ifelse(temp$r > 0, 0, 1), color = "black") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  #axis.text.x = element_blank(),  # 隐藏 X 轴文本
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                  axis.title.x = element_blank(),  # 隐藏 X 轴标题
                  # axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,colour = "black")
                  axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) +  
            labs(y = "Spearman correlation", fill = "- Log10 P-value") +
            scale_fill_gradientn(colours = colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)[26:51]) 
            # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51])
          
          
          return(out_plot)
          
          # p1 <- pheatmap::pheatmap(mat = cor_df$r,
          #                          cluster_cols = T,
          #                          cluster_rows = F,
          #                          cellwidth = 20,
          #                          cellheight = 20,
          #                          border_color = NA,
          #                          treeheight_col = 0,
          #                          treeheight_row = 30, 
          #                          fontsize = 15,
          #                          main = "Heatmap of correlation",
          #                          silent = T,
          #                          show_colnames = T,
          #                          color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
          # )
          
          # p2 <- pheatmap::pheatmap(mat = cor_df$p,
          #                          cluster_cols = F,
          #                          cluster_rows = F,
          #                          display_numbers = ifelse(test = cor_df$p > (-log10(0.001)),
          #                                                   yes =  "***",
          #                                                   no =  ifelse(test = cor_df$p > (-log10(0.01)),
          #                                                                yes =  "**",
          #                                                                no =  ifelse(test = cor_df$p > (-log10(0.05)),
          #                                                                             yes =  "*",
          #                                                                             no =  ""))),
          #                          cellwidth = 20,
          #                          cellheight = 20,
          #                          border_color = NA,
          #                          treeheight_col = 0,
          #                          treeheight_row = 0,
          #                          fontsize = 15,
          #                          main = "Heatmap of P-value (-log10)",
          #                          silent = T,
          #                          show_colnames = T
          # )
          
          
          # cor_df$heatmap <- plot_grid(plotlist = list(as.ggplot(p1), as.ggplot(p2)),ncol = 2,align = "hv")
          
          # cor_df$heatmap
          
        })
        
        plotOutput("plot",height="250px")
        
      } else {
        h1(paste0(input$data_from_html,' is not in this module'),align = "center")
      }
    }
  )
})



