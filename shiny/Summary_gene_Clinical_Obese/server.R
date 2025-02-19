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
      if (file.exists(paste0('./data/',gene_reactive(),'_RNA.RDS'))) {
        
        temp <- readRDS(paste0('./data/',gene_reactive(),'_RNA.RDS'))
        
        cohort_number <- reactiveVal()
        
        cohort_number(length(temp$n.c))
        
        
        output$plot1 <- renderPlot({
          
          temp <- readRDS(paste0('./data/',gene_reactive(),'_RNA.RDS'))

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
        
        
        output$plot2 <- renderPlot({
          gene <- gene_reactive()
          
          # 先处理SMD值小于-2和大于2的部分
          obese_df_RNA$SMD_binned <- ifelse(obese_df_RNA$SMD < -1, "<-1",
                                      ifelse(obese_df_RNA$SMD > 1, ">1", NA))
          
          # 处理其余部分，将其分成50个等距离的bin
          bins <- 50
          bin_breaks <- seq(-1, 1, length.out = bins + 1)
          obese_df_RNA$SMD_binned[is.na(obese_df_RNA$SMD_binned)] <- cut(obese_df_RNA$SMD[is.na(obese_df_RNA$SMD_binned)], 
                                                             breaks = bin_breaks, include.lowest = TRUE)
          
          # 计算每个bin中的计数
          bin_counts <- table(obese_df_RNA$SMD_binned)
          
          # 将数据转换为数据框以便使用ggplot2绘图
          bin_data <- as.data.frame(bin_counts)
          colnames(bin_data) <- c("bin", "count")
          bin_data$bin <- as.character(bin_data$bin)
          
          
          bin_data$bin <- factor(bin_data$bin,levels = c("<-1",1:bins,">1"))
          
          bin_data <- bin_data[order(bin_data$bin),]
          
          
          bin_data$start <- c(-Inf,seq(-1,1,(1--1)/bins),Inf)[1:52]
          bin_data$end <- c(-Inf,seq(-1,1,(1--1)/bins),Inf)[2:53]
          
          bin_data$group <- "out"
          bin_data$group[bin_data$start < obese_df_RNA$SMD[rownames(obese_df_RNA)==gene] & bin_data$end > obese_df_RNA$SMD[rownames(obese_df_RNA)==gene]] <- "in"
          bin_data$group <- factor(bin_data$group,levels = c("in","out"))
          
          
          bin_data$classfication <- "no"
          
          bin_data$classfication[bin_data$end<=-0.2|bin_data$start>=0.2] <- "small"
          bin_data$classfication[bin_data$end<=-0.5|bin_data$start>=0.5] <- "medium"
          bin_data$classfication[bin_data$end<=-0.8|bin_data$start>=0.8] <- "large"
          
          bin_data$classfication <- factor(bin_data$classfication,levels = c("no", "small", "medium", "large"))
          
          bin_data$rank <- 1:nrow(bin_data)
          
          
          p <- ggbarplot(bin_data, x = "rank", y = "count", fill = "group", color = NA, xlab = "SMD", ylab = "Count") +
            scale_fill_manual(values = c("in"="red","out"="grey"),name="Group") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  legend.text = element_text(size = 10, family = "Red Hat Display"),
                  axis.line.x = element_line(),  # 启用x轴线
                  axis.ticks.x = element_line(),  # 启用x轴刻度
                  axis.text.x = element_text()) +  # 启用x轴文本
            guides(fill = guide_legend(nrow = 4), override.aes = list(size=4))
          
          legend1 <- as_ggplot(get_legend(p))
          
          group_colors <- c("no" = "#7dbfb3", "small" = "#8ca5a5", "medium" = "#4f736a", "large" = "#1a4659", "in"="red", "out"="grey")
          annotation_bar <- data.frame(
            xmin = 0:51 + 0.5,
            xmax = 1:52 + 0.5,
            ymin = -Inf,
            ymax = -max(bin_data$count) * 0.02,
            group = bin_data$classfication
          )
          
          
          # 在条形图下方添加注释条
          p <- p + geom_rect(data = annotation_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), inherit.aes = FALSE) +
            scale_fill_manual(values = group_colors, name = "Effect Size") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"))+
            coord_cartesian(ylim = c(-max(bin_data$count) * 0.02, NA))
          
          # legend2 <- as_ggplot(get_legend(p + geom_rect(data = annotation_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), inherit.aes = FALSE) +
          #                                   scale_fill_manual(values = c("no" = "#7dbfb3", "small" = "#8ca5a5", "medium" = "#4f736a", "large" = "#1a4659"), name = "Effect Size") +
          #                                   coord_cartesian(ylim = c(-max(bin_data$count) * 0.02, NA))))
          
          legend2 <- as_ggplot(get_legend(ggbarplot(data = data.frame(value=c(1,2,3,4),group=factor(c("no", "small","medium","large"), levels = c("no", "small","medium","large"))),x = "group", y = "value",fill = "group", color = NA) +
                                            scale_fill_manual(values = c("no" = "#7dbfb3", "small" = "#8ca5a5", "medium" = "#4f736a", "large" = "#1a4659"), name = "Effect Size") +
                                            theme(legend.position = "right",
                                                  text = element_text(family = "Red Hat Display"),
                                                  legend.text = element_text(size = 10, family = "Red Hat Display")) +
                                            guides(fill = guide_legend(nrow = 4), override.aes = list(size=4)))
          )
          
          
          
          # 添加x轴标题和标签
          p <- p + 
            scale_x_continuous(breaks = c(1,6.5,13.5,20.5,26.5,32.5,39.5,46.5,52), labels = c("<-1",-0.8,-0.5,-0.2,0,0.2,0.5,0.8,">1"), name = "SMD", expand = c(0.02, 0)) +
            theme(legend.position = "none",
                  # legend.title = element_blank()
                  # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)
            )
          
          
          
          # 假设 p 是您的主图
          # 创建图例，这里是简化的例子，需要您根据实际情况调整
          legend1_grob <- ggplotGrob(legend1)
          legend2_grob <- ggplotGrob(legend2)
          
          # 调整图例尺寸，确保它们适合放置在图中
          legend1_grob$widths <- unit(rep(0.5, length(legend1_grob$widths)), "npc")
          legend2_grob$widths <- unit(rep(0.5, length(legend2_grob$widths)), "npc")
          
          
          p <- p + 
            # annotation_custom(grob = legend1_grob, 
            #                   xmin = layer_scales(p, 1)$x$get_limits()[1] + (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.0, 
            #                   xmax = layer_scales(p, 1)$x$get_limits()[1] + (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.2,
            #                   ymax = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.1,
            #                   ymin = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.2) +
            annotation_custom(grob = legend2_grob, 
                              xmin = layer_scales(p, 1)$x$get_limits()[2] - (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.3, 
                              xmax = layer_scales(p, 1)$x$get_limits()[2] - (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.0, 
                              ymax = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.1,
                              ymin = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.2)
          
          p
          
        })
        
        
        # gene_protein <- intersect(gene_reactive(), rownames(Proteomicseset))

        if (file.exists(paste0('./data/',gene_reactive(),'_protein.RDS'))) {
          
          output$plot3 <- renderPlot({

            # df <- rbind(
            #   data.frame(expression = Biobase::exprs(Proteomicseset)[gene_protein, colnames(Proteomicseset)[Proteomicseset$BMI_catelogy %in% c("Underweight", "Healthy", "Overweight")]], BMI_catelogy = "Non-obese"),
            #   data.frame(expression = Biobase::exprs(Proteomicseset)[gene_protein, colnames(Proteomicseset)[Proteomicseset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]], BMI_catelogy = "Obese")
            # ) %>% na.omit()
            # 
            # ggboxplot(data = df, x = "BMI_catelogy", y = "expression", add = "dotplot", fill = "BMI_catelogy", xlab = "", ylab = paste0(gene_protein, " expression")) +
            #   theme(axis.title.x = element_blank(),
            #         # aspect.ratio = 1, 
            #         text = element_text(family = "Red Hat Display"),
            #         legend.position = "none") + 
            #   stat_compare_means(comparisons = list(c("Non-obese", "Obese")), method = "wilcox.test") +
            #   scale_fill_manual(values = portalcol2)
            
            
            temp <- readRDS(paste0('./data/',gene_reactive(),'_protein.RDS'))
            
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
          
          
          output$plot4 <- renderPlot({
            gene <- gene_reactive()
            
            # 先处理SMD值小于-2和大于2的部分
            obese_df_protein$SMD_binned <- ifelse(obese_df_protein$SMD < -1, "<-1",
                                              ifelse(obese_df_protein$SMD > 1, ">1", NA))
            
            # 处理其余部分，将其分成50个等距离的bin
            bins <- 50
            bin_breaks <- seq(-1, 1, length.out = bins + 1)
            obese_df_protein$SMD_binned[is.na(obese_df_protein$SMD_binned)] <- cut(obese_df_protein$SMD[is.na(obese_df_protein$SMD_binned)], 
                                                                           breaks = bin_breaks, include.lowest = TRUE)
            
            # 计算每个bin中的计数
            bin_counts <- table(obese_df_protein$SMD_binned)
            
            # 将数据转换为数据框以便使用ggplot2绘图
            bin_data <- as.data.frame(bin_counts)
            colnames(bin_data) <- c("bin", "count")
            bin_data$bin <- as.character(bin_data$bin)
            
            
            bin_data$bin <- factor(bin_data$bin,levels = c("<-1",1:bins,">1"))
            
            bin_data <- bin_data[order(bin_data$bin),]
            
            
            bin_data$start <- c(-Inf,seq(-1,1,(1--1)/bins),Inf)[1:52]
            bin_data$end <- c(-Inf,seq(-1,1,(1--1)/bins),Inf)[2:53]
            
            bin_data$group <- "out"
            bin_data$group[bin_data$start < obese_df_protein$SMD[rownames(obese_df_protein)==gene] & bin_data$end > obese_df_protein$SMD[rownames(obese_df_protein)==gene]] <- "in"
            bin_data$group <- factor(bin_data$group,levels = c("in","out"))
            
            
            bin_data$classfication <- "no"
            
            bin_data$classfication[bin_data$end<=-0.2|bin_data$start>=0.2] <- "small"
            bin_data$classfication[bin_data$end<=-0.5|bin_data$start>=0.5] <- "medium"
            bin_data$classfication[bin_data$end<=-0.8|bin_data$start>=0.8] <- "large"
            
            bin_data$classfication <- factor(bin_data$classfication,levels = c("no", "small", "medium", "large"))
            
            bin_data$rank <- 1:nrow(bin_data)
            
            
            p <- ggbarplot(bin_data, x = "rank", y = "count", fill = "group", color = NA, xlab = "SMD", ylab = "Count") +
              scale_fill_manual(values = c("in"="red","out"="grey"),name="Group") +
              theme(legend.position = "right",
                    text = element_text(family = "Red Hat Display"),
                    legend.text = element_text(size = 10, family = "Red Hat Display"),
                    axis.line.x = element_line(),  # 启用x轴线
                    axis.ticks.x = element_line(),  # 启用x轴刻度
                    axis.text.x = element_text()) +  # 启用x轴文本
              guides(fill = guide_legend(nrow = 4), override.aes = list(size=4)) +
              scale_y_continuous(breaks = c(0,50,100,150,200), labels = c("     0","      50","     100","     150","     200"))
            
            legend1 <- as_ggplot(get_legend(p))
            
            group_colors <- c("no" = "#7dbfb3", "small" = "#8ca5a5", "medium" = "#4f736a", "large" = "#1a4659", "in"="red", "out"="grey")
            annotation_bar <- data.frame(
              xmin = 0:51 + 0.5,
              xmax = 1:52 + 0.5,
              ymin = -Inf,
              ymax = -max(bin_data$count) * 0.02,
              group = bin_data$classfication
            )
            
            
            # 在条形图下方添加注释条
            p <- p + geom_rect(data = annotation_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), inherit.aes = FALSE) +
              scale_fill_manual(values = group_colors, name = "Effect Size") +
              theme(legend.position = "right",
                    text = element_text(family = "Red Hat Display"))+
              coord_cartesian(ylim = c(-max(bin_data$count) * 0.02, NA))
            
            # legend2 <- as_ggplot(get_legend(p + geom_rect(data = annotation_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), inherit.aes = FALSE) +
            #                                   scale_fill_manual(values = c("no" = "#7dbfb3", "small" = "#8ca5a5", "medium" = "#4f736a", "large" = "#1a4659"), name = "Effect Size") +
            #                                   coord_cartesian(ylim = c(-max(bin_data$count) * 0.02, NA))))
            
            legend2 <- as_ggplot(get_legend(ggbarplot(data = data.frame(value=c(1,2,3,4),group=factor(c("no", "small","medium","large"), levels = c("no", "small","medium","large"))),x = "group", y = "value",fill = "group", color = NA) +
                                              scale_fill_manual(values = c("no" = "#7dbfb3", "small" = "#8ca5a5", "medium" = "#4f736a", "large" = "#1a4659"), name = "Effect Size") +
                                              theme(legend.position = "right",
                                                    text = element_text(family = "Red Hat Display"),
                                                    legend.text = element_text(size = 10, family = "Red Hat Display")) +
                                              guides(fill = guide_legend(nrow = 4), override.aes = list(size=4)))
            )
            
            
            
            # 添加x轴标题和标签
            p <- p + 
              scale_x_continuous(breaks = c(1,6.5,13.5,20.5,26.5,32.5,39.5,46.5,52), labels = c("<-1",-0.8,-0.5,-0.2,0,0.2,0.5,0.8,">1"), name = "SMD", expand = c(0.02, 0)) +
              theme(legend.position = "none",
                    # legend.title = element_blank()
                    # axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5)
              )
            
            
            
            # 假设 p 是您的主图
            # 创建图例，这里是简化的例子，需要您根据实际情况调整
            legend1_grob <- ggplotGrob(legend1)
            legend2_grob <- ggplotGrob(legend2)
            
            # 调整图例尺寸，确保它们适合放置在图中
            legend1_grob$widths <- unit(rep(0.5, length(legend1_grob$widths)), "npc")
            legend2_grob$widths <- unit(rep(0.5, length(legend2_grob$widths)), "npc")
            
            
            p <- p + 
              # annotation_custom(grob = legend1_grob, 
              #                   xmin = layer_scales(p, 1)$x$get_limits()[1] + (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.0, 
              #                   xmax = layer_scales(p, 1)$x$get_limits()[1] + (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.2,
              #                   ymax = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.1,
              #                   ymin = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.2) +
              annotation_custom(grob = legend2_grob, 
                                xmin = layer_scales(p, 1)$x$get_limits()[2] - (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.3, 
                                xmax = layer_scales(p, 1)$x$get_limits()[2] - (layer_scales(p, 1)$x$get_limits()[2]-layer_scales(p, 1)$x$get_limits()[1])*0.0, 
                                ymax = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.1,
                                ymin = layer_scales(p, 1)$y$get_limits()[2] - (layer_scales(p, 1)$y$get_limits()[2]-layer_scales(p, 1)$y$get_limits()[1])*0.2)
            
            p
            
          })

          fluidRow(
            column(width = 8,
                   class = "responsive-col",
                   style = 'padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                   "Transcriptomics",
                   plotOutput("plot1", height = paste0(cohort_number()*15+135,"px"))
            ),
            column(width = 4,
                   class = "responsive-col",
                   style = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                   br(),
                   plotOutput("plot2", height = paste0(cohort_number()*15+135,"px"))
            ),
            column(width = 12,
                   style = 'padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                   br()
            ),
            column(width = 8,
                   class = "responsive-col",
                   style='padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                   "Proteomics",
                   plotOutput("plot3", height = paste0(cohort_number()*15+135,"px"))
            )
            ,
            column(width = 4,
                   class = "responsive-col",
                   style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                   br(),
                   plotOutput("plot4", height = paste0(cohort_number()*15+135,"px"))
            )
          )
        } else {
          fluidRow(
            column(width = 8,
                   class = "responsive-col",
                   style='padding-left:0px; padding-right:15px; padding-top:0px; padding-bottom:0px',
                   plotOutput("plot1", height = paste0(cohort_number()*15+135,"px"))
            ),
            column(width = 4,
                   class = "responsive-col",
                   style='padding-left:0px; padding-right:0px; padding-top:0px; padding-bottom:0px',
                   plotOutput("plot2", height = paste0(cohort_number()*15+135,"px"))
            )
          )
          
        }
                
      } else {
        h1(paste0(gene_reactive(),' is not in this module'),align = "center")
      }
    }
  )
})


