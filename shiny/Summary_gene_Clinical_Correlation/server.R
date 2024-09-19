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
      h1("undefined", align = "center")
    } else {
      if (gene_reactive() %in% rownames(meta_r)) {
        
        gene_protein <- intersect(gene_reactive(), rownames(protein_r))
        
        if (length(gene_protein)==0) {
          output$plot <- renderPlot({
            
            draw_plot <- function(data, max_p, min_p) {
              if(nrow(data) == 0) return(NULL)
              ggbarplot(data, x = "trait_full", y = "r", fill = "n", color = NA, rotate = TRUE) +
                ylim(-max(abs(data$r))*1.3,max(abs(data$r))*1.3) +
                geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
                geom_text(aes(label = p_stars, y = ifelse(r > 0, r + max(abs(data$r))*0.12, r - max(abs(data$r))*0.12)), 
                          vjust = 0.5, color = "black", hjust = 0.5) +
                theme(legend.position = "none",
                      text = element_text(family = "Red Hat Display"),
                      # axis.text.x = element_blank(),  # 隐藏 X 轴文本
                      # axis.line.x = element_blank(),   # 隐藏 X 轴线
                      # axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                      # axis.title.x = element_blank(),  # 隐藏 X 轴标题
                      
                      # axis.text.y = element_blank(),  # 隐藏 Y 轴文本
                      axis.line.y = element_blank(),  # 隐藏 Y 轴线
                      axis.ticks.y = element_blank(),  # 隐藏 Y 轴刻度
                      axis.title.y = element_blank()  # 隐藏 Y 轴标题
                ) +
                # labs(y = "Spearman correlation", fill = "- Log10 P-value") +
                labs(y = "Spearman correlation", fill = "Number of cohort") +
                scale_fill_gradientn(colours = colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)[26:51], limit=c(min_p, max_p))  
                # theme(aspect.ratio = 1)
              # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51])
            }
            
            cor_df <- list()
            cor_df$r <- meta_r[gene_reactive(),,drop=F]
            cor_df$p <- meta_p[gene_reactive(),,drop=F]
            cor_df$n <- meta_n[gene_reactive(),,drop=F]
            cor_df$p <- (-log10(cor_df$p))
            
            
            
            temp <- data.frame(trait=colnames(cor_df$r),
                               trait_full=names(trait_vector)[match(colnames(cor_df$r), trait_vector)],
                               r=cor_df$r[1,],
                               p=cor_df$p[1,],
                               n=cor_df$n[1,],
                               p_stars=ifelse(cor_df$p[1,] < -log10(0.05), "ns", ifelse(cor_df$p[1,] < -log10(0.01), "*", ifelse(cor_df$p[1,] < -log10(0.001), "**", ifelse(cor_df$p[1,] < -log10(0.0001), "***", "****")))))
            
            temp <- temp[order(temp$r,decreasing = F),] %>% na.omit()
            
            temp$type[temp$trait %in% c("BMI", "HOMA", "Age", "WHR", "Hip", "Waist")] <- "anthropometrical"
            temp$type[temp$trait %in% c("Glucose", "Insulin", "TG", "Chol", "HDL", "LDL", "CRP", "Hba1c")] <- "circulating"
            temp$type[temp$trait %in% c("LEP_protein", "TNF_protein", "MCP1_protein", "cell_vol", "basal_TG", "iso_TG", "iso_basal")] <- "WAT"
            
            
            # max_p <- max(temp$p, na.rm = TRUE)
            # min_p <- min(temp$p, na.rm = TRUE)
            
            max_p <- max(temp$n, na.rm = TRUE)
            min_p <- min(temp$n, na.rm = TRUE)
            
            
            
            
            p_list <- list(
              draw_plot(temp[temp$type == "anthropometrical", ], max_p, min_p),
              draw_plot(temp[temp$type == "circulating", ], max_p, min_p),
              draw_plot(temp[temp$type == "WAT", ], max_p, min_p)
            )
            
            # 移除空图形
            p_list <- p_list[!sapply(p_list, is.null)]
            
            # 设置最后一个非空图形的图例位置为"right"
            # p_list[[length(p_list)]] <- p_list[[length(p_list)]] + theme(legend.position = "right")
            
            legend <- as_ggplot(get_legend(p_list[[length(p_list)]] + theme(legend.position = "right")))
            
            final_plot <- plot_grid(plotlist = p_list, align = "hv", nrow = 1, rel_heights = c(1,1,1))
            
            # egg::ggarrange(plots = p_list,nrow = 1)
            ggpubr::ggarrange(final_plot, legend,nrow = 1, widths = c(10,1))
            
            
          })
        } else {
          output$plot <- renderPlot({
            
            draw_plot <- function(data, max_p, min_p) {
              if(nrow(data) == 0) return(ggplot() + theme_void()+ theme(aspect.ratio=0.5, text = element_text(family = "Red Hat Display")))
              ggbarplot(data, x = "trait_full", y = "r", fill = "n", color = NA, rotate = TRUE) +
                ylim(-max(abs(data$r))*1.35,max(abs(data$r))*1.35) +
                geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
                geom_text(aes(label = p_stars, y = ifelse(r > 0, r + max(abs(data$r))*0.12, r - max(abs(data$r))*0.12)), 
                          vjust = 0.5, color = "black", hjust = 0.5) +
                theme(legend.position = "none",
                      text = element_text(family = "Red Hat Display"),
                      # axis.text.x = element_blank(),  # 隐藏 X 轴文本
                      # axis.line.x = element_blank(),   # 隐藏 X 轴线
                      # axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                      # axis.title.x = element_blank(),  # 隐藏 X 轴标题
                      
                      # axis.text.y = element_blank(),  # 隐藏 Y 轴文本
                      axis.line.y = element_blank(),  # 隐藏 Y 轴线
                      axis.ticks.y = element_blank(),  # 隐藏 Y 轴刻度
                      axis.title.y = element_blank()  # 隐藏 Y 轴标题
                ) +
                # labs(y = "Spearman correlation", fill = "- Log10 P-value") +
                labs(y = "Spearman correlation", fill = "Number of cohort") +
                scale_fill_gradientn(colours = colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)[26:51], limit=c(min_p, max_p))  
              # theme(aspect.ratio = 0.5)
              # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51])
            }
            
            
            cor_df <- list()
            cor_df$r <- meta_r[gene_reactive(),,drop=F]
            cor_df$p <- meta_p[gene_reactive(),,drop=F]
            cor_df$n <- meta_n[gene_reactive(),,drop=F]
            cor_df$p <- (-log10(cor_df$p))
            
            
            
            temp_RNA <- data.frame(trait=colnames(cor_df$r),
                                   trait_full=names(trait_vector)[match(colnames(cor_df$r), trait_vector)],
                                   r=cor_df$r[1,] %>% as.vector(),
                                   p=cor_df$p[1,] %>% as.vector(),
                                   n=cor_df$n[1,] %>% as.vector(),
                                   p_stars=ifelse(cor_df$p[1,] < -log10(0.05), "ns", ifelse(cor_df$p[1,] < -log10(0.01), "*", ifelse(cor_df$p[1,] < -log10(0.001), "**", ifelse(cor_df$p[1,] < -log10(0.0001), "***", "****")))))
            
            temp_RNA <- temp_RNA[order(temp_RNA$r,decreasing = F),] %>% na.omit()
            
            temp_RNA$type[temp_RNA$trait %in% c("BMI", "HOMA", "Age", "WHR", "Hip", "Waist")] <- "anthropometrical"
            temp_RNA$type[temp_RNA$trait %in% c("Glucose", "Insulin", "TG", "Chol", "HDL", "LDL", "CRP", "Hba1c")] <- "circulating"
            temp_RNA$type[temp_RNA$trait %in% c("LEP_protein", "TNF_protein", "MCP1_protein", "cell_vol", "basal_TG", "iso_TG", "iso_basal")] <- "WAT"
            
            
            cor_df <- list()
            cor_df$r <- protein_r[gene_protein,,drop=F]
            cor_df$p <- protein_p[gene_protein,,drop=F]
            cor_df$n <- protein_n[gene_protein,,drop=F]
            cor_df$p <- (-log10(cor_df$p))
            
            
            
            temp_protein <- data.frame(trait=colnames(cor_df$r),
                                       trait_full=names(trait_vector)[match(colnames(cor_df$r), trait_vector)],
                                       r=cor_df$r[1,] %>% as.vector(),
                                       p=cor_df$p[1,] %>% as.vector(),
                                       n=cor_df$n[1,] %>% as.vector(),
                                       p_stars=ifelse(cor_df$p[1,] < -log10(0.05), "ns", ifelse(cor_df$p[1,] < -log10(0.01), "*", ifelse(cor_df$p[1,] < -log10(0.001), "**", ifelse(cor_df$p[1,] < -log10(0.0001), "***", "****")))))
            
            temp_protein <- temp_protein[order(temp_protein$r,decreasing = F),] %>% na.omit()
            
            temp_protein$type[temp_protein$trait %in% c("BMI", "HOMA", "Age", "WHR", "Hip", "Waist")] <- "anthropometrical"
            temp_protein$type[temp_protein$trait %in% c("Glucose", "Insulin", "TG", "Chol", "HDL", "LDL", "CRP", "Hba1c")] <- "circulating"
            temp_protein$type[temp_protein$trait %in% c("LEP_protein", "TNF_protein", "MCP1_protein", "cell_vol", "basal_TG", "iso_TG", "iso_basal")] <- "WAT"
            
            
            
            
            max_p <- max(c(temp_RNA$n,temp_protein$n), na.rm = TRUE)
            min_p <- min(c(temp_RNA$n,temp_protein$n), na.rm = TRUE)
            
            
            
            
            
            p1 <- draw_plot(temp_RNA[temp_RNA$type == "anthropometrical", ], max_p, min_p) + ggtitle("Transcriptomics") + theme(plot.title.position = "plot")
            p2 <- draw_plot(temp_RNA[temp_RNA$type == "circulating", ], max_p, min_p) + theme(plot.title.position = "plot")
            p3 <- draw_plot(temp_RNA[temp_RNA$type == "WAT", ], max_p, min_p) + theme(plot.title.position = "plot")
            p4 <- draw_plot(temp_protein[temp_protein$type == "anthropometrical", ], max_p, min_p) + ggtitle("Proteomics") + theme(plot.title.position = "plot")
            p5 <- draw_plot(temp_protein[temp_protein$type == "circulating", ], max_p, min_p) 
            p6 <- draw_plot(temp_protein[temp_protein$type == "WAT", ], max_p, min_p) 
            
            
            
            
            # # 移除空图形
            # p_list <- p_list[!sapply(p_list, is.null)]
            # 
            # # 设置最后一个非空图形的图例位置为"right"
            # p_list[[length(p_list)]] <- p_list[[length(p_list)]] + theme(legend.position = "right")
            
            legend <- as_ggplot(get_legend(ggbarplot(temp_RNA, x = "trait_full", y = "r", fill = "n", color = NA, rotate = TRUE)+labs(y = "Spearman correlation", fill = "Number of cohort") +
                                             scale_fill_gradientn(colours = colorRampPalette(colors = c("#1C4759","#EEF2F2", "#E2C744"))(51)[26:51], limit=c(min_p, max_p)) +
                                             theme(legend.position = "right",
                                                   text = element_text(family = "Red Hat Display"))))
            
            # final_plot <- ggpubr::ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6),align = "hv")
            final_plot <- plot_grid(plotlist = list(p1,p2,p3, NULL, NULL, NULL, p4, p5, p6),align = "hv",rel_heights = c(10,1,10))
            
            # final_plot <- as_ggplot(egg::ggarrange(plots = list(p1,p2,p3,p4,p5,p6),
            #                                        align = "hv",nrow = 2,
            #                                        # labels = c("Transcriptomics", "", "", "Proteomics\n\n", "\n\n", "\n\n"),
            #                                        label.args = list(gp=gpar(font=1,fontfamily="Red Hat Display",fontsize=15), x=unit(0,"line"), hjust=0)))
            
            ggpubr::ggarrange(final_plot,legend,nrow = 1,widths = c(10,1))
            
          })
        }
        
        fluidRow(
          if (length(gene_protein)==0) {
            plotOutput("plot",height="190px")
          } else {
            plotOutput("plot",height="380px")
          }
          
        )
        
      } else {
        h1(paste0(gene_reactive(),' is not in this module'),align = "center")
      }
    }
  )
})
