library(shiny)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(viridis)
library(DT)
library(Biobase)
library(gplots)
library(dplyr)
library(cowplot)
library(fresh)

#
shinyServer(function(input, output) {
    
    v <- reactiveValues(#data_ins=NULL,
                        data_hyp2 = NULL,
                        data_tnf=NULL,
                        data_gls=NULL)
    data <- reactiveValues()
    observeEvent(input$start, {
        gene_split <- strsplit(input$gene_id, split = "\n")  
        v$data_hyp2 <- data_hyp_tt[data_hyp_tt$GENE %in% gene_split[[1]],]
        #v$data_ins <- data_ins_tt[data_ins_tt$gene %in% gene_split[[1]],]
        v$data_tnf <- data_tnf_tt[data_tnf_tt$gene %in% gene_split[[1]],]
        v$data_gls <- data_gls_tt[data_gls_tt$gene %in% gene_split[[1]],]
    })


#outputs hypoxia
    output$DT_hypoxia <- DT::renderDataTable({
        if (is.null(v$data_hyp2)) return()
        datatable(v$data_hyp2,extensions = 'Buttons', 
        options = list(dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  %>% 
          formatSignif(columns = c('P.Value', 'adj.P.Val'), digits = 3) %>%
          formatRound(columns=c('logFC', 'AveExpr',"t","B"), digits=3)
    })
    
    output$heatmap_hypoxia <- renderPlot({
      if (is.null(v$data_hyp2)) return()
      if (input$mode == "normalized") {
        data$p5 <- pheatmap(data_hyp[rownames(v$data_hyp2),], scale = "row", cluster_cols = FALSE, color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                            border_color = NA, main = "Hypoxia", labels_row = v$data_hyp2$GENE)
        dev.off() 
        print(data$p5)
        
      } else {
        data$p5 <- pheatmap(data_hyp[rownames(v$data_hyp2),], cluster_cols = FALSE, color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                            border_color = NA, main = "Hypoxia", labels_row = v$data_hyp2$GENE)
        dev.off() 
        print(data$p5)
      }
    })
    
    output$volcano_hypoxia <- renderPlotly({
        if (is.null(v$data_hyp2)) return() 
        p <- ggplot(data=data_hyp_tt, aes(x=logFC, y=-log10(adj.P.Val), label=GENE)) +
            geom_point(color="#d2d3d4") + geom_point(data = v$data_hyp2, aes(x=logFC, y=-log10(adj.P.Val)), color="#57ac82") +
            theme_classic() +
            #geom_text_repel(data = v$data_hyp2) +
            geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
        ggplotly(p)
    })

    #outputs TNF    
    output$heatmap_tnf <- renderPlot({
      if (is.null(v$data_tnf)) return()
      if (input$mode == "normalized") {
      data$p7 <- pheatmap(data_tnf[rownames(data_tnf) %in% rownames(v$data_tnf),1:6], scale = "row", color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                          border_color = NA, main = "TNF stimulation", labels_row =  data_tnf[rownames(data_tnf) %in% rownames(v$data_tnf),7])
      dev.off()
      print(data$p7)
    } else {
      data$p7 <- pheatmap(data_tnf[rownames(data_tnf) %in% rownames(v$data_tnf),1:6], color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                          border_color = NA, main = "TNF stimulation", labels_row =  data_tnf[rownames(data_tnf) %in% rownames(v$data_tnf),7])
      dev.off() 
      print(data$p7)
    }
    })
    
    output$volcano_tnf <- renderPlotly({
      if (is.null(v$data_tnf)) return() 
      p <- ggplot(data=data_tnf_tt, aes(x=log2FoldChange, y=-log10(padj), label=gene)) +
        geom_point(color="#d2d3d4") + geom_point(data = v$data_tnf, aes(x=log2FoldChange, y=-log10(padj)), color="#57ac82") +
        theme_classic() +
        geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
      ggplotly(p)
    })
    output$DT_tnf <- DT::renderDataTable({
      if (is.null(v$data_tnf)) return()
      datatable(v$data_tnf,extensions = 'Buttons', 
                options = list(dom = 'Bfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  %>% 
        formatSignif(columns = c('pvalue', 'padj','bonf'), digits = 3) %>%
        formatRound(columns=c("baseMean","log2FoldChange","lfcSE",'pvalue', 'padj','bonf'), digits=3)
      
    })
    #outputs siGLS    
    output$heatmap_gls <- renderPlot({
      if (is.null(v$data_gls)) return()
      data$p8 <- pheatmap(data_gls[rownames(data_gls) %in% rownames(v$data_gls),1:6], scale = "row", color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                          border_color = NA, main = "INS stimulation", labels_row =  data_gls[rownames(data_gls) %in% rownames(v$data_gls),7])
      dev.off()
      print(data$p8)
    })
    
    output$volcano_gls <- renderPlotly({
      if (is.null(v$data_gls)) return() 
      p <- ggplot(data=data_gls_tt, aes(x=log2FoldChange, y=-log10(padj), label=gene)) +
        geom_point(color="#d2d3d4") + geom_point(data = v$data_gls, aes(x=log2FoldChange, y=-log10(padj)), color="#57ac82") +
        theme_classic() +
        geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
      ggplotly(p)
    })
    output$DT_gls <- DT::renderDataTable({
      if (is.null(v$data_gls)) return()
      datatable(v$data_gls,extensions = 'Buttons', 
                options = list(dom = 'Bfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  %>% 
        formatSignif(columns = c('pvalue', 'padj','bonf'), digits = 3) %>%
        formatRound(columns=c("baseMean","log2FoldChange","lfcSE",'pvalue', 'padj','bonf'), digits=3)
      
    })
    

    })


        

