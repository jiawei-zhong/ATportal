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
      if (gene_reactive() %in% c(rownames(DEOSHeset_WL),rownames(POeset_WL), rownames(GSE141221_diogenes1_WL), rownames(GSE95640_diogenes2_WL))) {
        
        output$plot <- renderPlot({
          if (gene_reactive() %in% rownames(DEOSHeset_WL)) {
            temp_pdat <- pData(DEOSHeset_WL)
            temp_pdat$expression <- exprs(DEOSHeset_WL)[gene_reactive(),]
            temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
            temp_pdat <- temp_pdat[order(temp_pdat$subject),]
            temp_pdat$expression <- scale_values(temp_pdat$expression)
            p1 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",fill = "time_point",xlab = "years after bariatric surgery",ylab = paste0(gene_reactive(), " expression"), title = "Kerr, A. (2020)")+
              # geom_line(aes(group=subject), linetype = "dashed") +
              stat_compare_means(comparisons = list(c("0","2"),c("2","5"),c("0","5")), paired = TRUE, method = "wilcox.test") +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme_set(theme_bw()) +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
                    axis.line.y=element_line(linetype=1,color="black",size=0.9),
                    axis.title.y = element_text(size = 13),
                    axis.text.x = element_text(face = "plain",size = 12,colour = "black"),
                    axis.text.y = element_text(face = "plain",size = 12,colour = "black")) + NoLegend() + theme(aspect.ratio=1, text = element_text(family = "Red Hat Display")) +
              scale_fill_manual(values = c("#1a4659", "#e2c744", "#ba5c28"))
          } else {
            p1 <- ggplot() + theme_void()+ theme(aspect.ratio=1, text = element_text(family = "Red Hat Display"))
          }
          
          if (gene_reactive() %in% rownames(POeset_WL)) {
            temp_pdat <- pData(POeset_WL)
            temp_pdat$expression <- exprs(POeset_WL)[gene_reactive(),]
            temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
            temp_pdat <- temp_pdat[order(temp_pdat$subject),]
            temp_pdat$expression <- scale_values(temp_pdat$expression)
            p2 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",fill = "time_point",xlab = "years after bariatric surgery",ylab = paste0(gene_reactive(), " expression"), title = "Petrus, P. (2018)") +
              # geom_line(aes(group=subject), linetype = "dashed") +
              stat_compare_means(comparisons = list(c("0","2")), paired = TRUE, method = "wilcox.test") +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme_set(theme_bw()) +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
                    axis.line.y=element_line(linetype=1,color="black",size=0.9),
                    axis.title.y = element_text(size = 13),
                    axis.text.x = element_text(face = "plain",size = 12,colour = "black"),
                    axis.text.y = element_text(face = "plain",size = 12,colour = "black")) + NoLegend() + theme(aspect.ratio=1, text = element_text(family = "Red Hat Display")) +
              scale_fill_manual(values = c("#1a4659", "#e2c744", "#ba5c28"))
          } else {
            p2 <- ggplot() + theme_void()+ theme(aspect.ratio=1, text = element_text(family = "Red Hat Display"))
          }

          if (gene_reactive() %in% rownames(GSE141221_diogenes1_WL)) {
            temp_pdat <- pData(GSE141221_diogenes1_WL)
            temp_pdat$expression <- exprs(GSE141221_diogenes1_WL)[gene_reactive(),]
            temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
            temp_pdat <- temp_pdat[order(temp_pdat$subject),]
            temp_pdat$expression <- scale_values(temp_pdat$expression)
            p3 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",fill = "time_point",xlab = "weeks after diet",ylab = paste0(gene_reactive(), " expression"), title = "Imbert, A. (2022)") +
              # geom_line(aes(group=subject), linetype = "dashed") +
              stat_compare_means(comparisons = list(c("0","8")), paired = TRUE, method = "wilcox.test") +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme_set(theme_bw()) +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
                    axis.line.y=element_line(linetype=1,color="black",size=0.9),
                    axis.title.y = element_text(size = 13),
                    axis.text.x = element_text(face = "plain",size = 12,colour = "black"),
                    axis.text.y = element_text(face = "plain",size = 12,colour = "black")) + NoLegend() + theme(aspect.ratio=1, text = element_text(family = "Red Hat Display")) +
              scale_fill_manual(values = c("#1a4659", "#e2c744", "#ba5c28"))
          } else {
            p3 <- ggplot() + theme_void()+ theme(aspect.ratio=1, text = element_text(family = "Red Hat Display"))
          }

          if (gene_reactive() %in% rownames(GSE95640_diogenes2_WL)) {
            temp_pdat <- pData(GSE95640_diogenes2_WL)
            temp_pdat$expression <- exprs(GSE95640_diogenes2_WL)[gene_reactive(),]
            temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
            temp_pdat <- temp_pdat[order(temp_pdat$subject),]
            temp_pdat$expression <- scale_values(temp_pdat$expression)
            p4 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",fill = "time_point",xlab = "weeks after diet",ylab = paste0(gene_reactive(), " expression"), title = "Armenise, C. (2017)") +
              # geom_line(aes(group=subject), linetype = "dashed") +
              stat_compare_means(comparisons = list(c("0","8")), paired = TRUE, method = "wilcox.test") +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme_set(theme_bw()) +
              theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
              theme(panel.border = element_blank()) +
              theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
                    axis.line.y=element_line(linetype=1,color="black",size=0.9),
                    axis.title.y = element_text(size = 13),
                    axis.text.x = element_text(face = "plain",size = 12,colour = "black"),
                    axis.text.y = element_text(face = "plain",size = 12,colour = "black")) + NoLegend() + theme(aspect.ratio=1, text = element_text(family = "Red Hat Display")) +
              scale_fill_manual(values = c("#1a4659", "#e2c744", "#ba5c28"))
          } else {
            p4 <- ggplot() + theme_void()+ theme(aspect.ratio=1, text = element_text(family = "Red Hat Display"))
          }
          
          # if (gene_reactive() %in% rownames(ADIPOSTRESS_WL)) {
          #   temp_pdat <- pData(ADIPOSTRESS_WL)
          #   temp_pdat$expression <- exprs(ADIPOSTRESS_WL)[gene_reactive(),]
          #   temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
          #   temp_pdat <- temp_pdat[order(temp_pdat$subject),]
          #   p3 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",add = "point",fill = "time_point",palette = "ucscgb",xlab = "weeks after diet",facet.by = "group",ylab = paste0(gene_reactive(), " expression"), title = "ADIPOSTRESS weight loss") +
          #     geom_line(aes(group=subject), linetype = "dashed") +
          #     stat_compare_means(comparisons = list(c("0","5"),c("5","9"),c("0","9")), paired = TRUE, method = "wilcox.test") +
          #     theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
          #     theme(panel.border = element_blank()) +
          #     theme_set(theme_bw()) +
          #     theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
          #     theme(panel.border = element_blank()) +
          #     theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
          #           axis.line.y=element_line(linetype=1,color="black",size=0.9),
          #           axis.title.y = element_text(size = 13),
          #           axis.text.x = element_text(face = "plain",size = 12,colour = "black"),
          #           axis.text.y = element_text(face = "plain",size = 12,colour = "black")) + NoLegend() + theme(aspect.ratio=2)
          # } else {
          #   p3 <- ggplot() + theme_void()+ theme(aspect.ratio=1)
          # }
          
          egg::ggarrange(p1, p2, p3, p4, nrow = 1)
          
        })
        
        plotOutput("plot", height = "400px")
        
      } else {
        h1(paste0(gene_reactive(),' is not in this module'),align = "center")
      }
    }
  )
})

