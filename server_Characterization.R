v <- reactiveValues(data=NULL, data2=NULL, dataNA=NULL, data_fa = NULL)
data <- reactiveValues()
observeEvent(input$start, {
  gene_split <- strsplit(input$gene_id, split = "\n")  
  v$data <- data_line[data_line$gene %in% gene_split[[1]],]
  v$data2 <- data_SVF[data_SVF$gene %in% gene_split[[1]],]
  v$dataNA <- data_novo[data_novo$gene %in% gene_split[[1]],]
  v$data_fa <- data_fantom[rownames(data_fantom) %in% gene_split[[1]],]
})
#NIGA timecourse    
output$linechart <- renderPlotly({
  if (is.null(v$data)) return()
  if (input$mode == "normalized") {
    data$p <- ggplot(data = v$data, aes(x=log(timepoint), y=normalized, group=gene, color=gene))+ 
      geom_point()+ geom_line()+ theme_classic() +
      scale_color_brewer(type = "qual", palette = input$discrete_col)
    ggplotly(data$p)    
  } else {
    data$p<-  ggplot(data = v$data, aes(x=log(timepoint), y=mean_abund, group=gene, color=gene))+ geom_point()+ geom_line()+ theme_classic()+
      scale_color_brewer(type = "qual", palette = input$discrete_col)
    ggplotly(data$p)    
  }
})
output$download_pdf1 <- downloadHandler(
  filename = function() {
    paste('timecourse', Sys.Date(), '.pdf', sep='')
  },
  content = function(file) {
    ggsave(file, plot = data$p, device = "pdf")
  })
output$download_gg1 <- downloadHandler(
  filename = function() {
    paste('timecourse', Sys.Date(), '.rds', sep='')
  },
  content = function(file) {
    saveRDS(file, object = data$p)
  }) 

output$heatmap <- renderPlot({
  if (is.null(v$data)) return()
  if (input$mode == "normalized") {
    data$p2 <- pheatmap((data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]), scale = "row", cluster_cols = FALSE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "Time course NIGA cells")
    dev.off() #somehow needed to plot correctly without resizing.
    print(data$p2)
    
  } else {
    data$p2 <- pheatmap((data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]), cluster_cols = FALSE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "Time course NIGA cells")
    dev.off() #somehow needed to plot correctly without resizing.
    print(data$p2)
  }
})

output$download_pdf2 <- downloadHandler(
  filename = function() {
    paste('timecourse_heat', Sys.Date(), '.pdf', sep='')
  },
  content = function(file) {
    ggsave(file, plot = data$p2, device = "pdf")
  })
output$download_gg2 <- downloadHandler(
  filename = function() {
    paste('timecourse_heat', Sys.Date(), '.rds', sep='')
  },
  content = function(file) {
    saveRDS(file, object = data$p2)
  }) 

output$DT_NIGA <- DT::renderDataTable({
  if (is.null(v$data)) return()
  v$data},extensions = 'Buttons', 
  options = list(dom = 'Bfrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))
)

#SVF
output$linechart2 <- renderPlotly({
  if (is.null(v$data2)) return()
  if (input$mode == "normalized") {
    data$p3 <- ggplot(data = v$data2, aes(x=timepoint, y=normalized, group=gene, color=gene))+ 
      geom_point()+ geom_line()+ theme_classic() +
      scale_color_brewer(type = "qual", palette = input$discrete_col)
    ggplotly(data$p3)    
  } else {
    data$p3<-  ggplot(data = v$data2, aes(x=timepoint, y=mean_abund, group=gene, color=gene))+ geom_point()+ geom_line()+ theme_classic() +
      scale_color_brewer(type = "qual", palette = input$discrete_col)
    ggplotly(data$p3)   
  }
})
output$heatmap2 <- renderPlot({
  if (is.null(v$data2)) return()
  if (input$mode == "normalized") {
    data$p4 <- pheatmap((data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]), scale = "row", cluster_cols = FALSE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "Time course SVF cells")
    dev.off() 
    print(data$p4)
    
  } else {
    data$p4 <- pheatmap((data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]), cluster_cols = FALSE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "Time course SVF cells")
    dev.off() 
    print(data$p4)
  }
})
output$DT_SVF<- DT::renderDataTable({
  if (is.null(v$data2)) return()
  v$data},extensions = 'Buttons', 
  options = list(dom = 'Bfrtip',
                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))
)

#novo    
output$heatmap_novo <-  renderPlot({
  if (is.null(v$dataNA)) return()
  if (input$mode == "normalized") {
    data$p5 <- pheatmap(data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8], scale = "row", cluster_cols = TRUE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "FACS sorted WAT celltypes", labels_row = data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,9])
    dev.off() 
    print(data$p5)
    
  } else {
    data$p5 <- pheatmap(data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8], cluster_cols = TRUE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "FACS sorted WAT celltypes", labels_row = data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,9])
    dev.off() 
    print(data$p5)
  }
})

output$table_novo <- DT::renderDataTable({
  if (is.null(v$dataNA)) return()
  datatable(v$dataNA,extensions = 'Buttons', 
            options = list(dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  %>% 
    formatRound(columns=c('mean', 'z_norm'), digits=3)
})

#FANTOM
output$heatmap_fantom <-  renderPlot({
  if (is.null(v$data_fa)) return()
  if (input$mode == "normalized") {
    data$p6 <- pheatmap(v$data_fa[,2:4], scale = "row", cluster_cols = TRUE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "Comparison FANTOM data")
    dev.off() 
    print(data$p6)
    
  } else {
    data$p6 <- pheatmap(v$data_fa[,2:4], cluster_cols = TRUE, color = farben[[input$gradient_col]], clustering_method = "ward.D",
                        border_color = NA, main = "Comparison FANTOM data")
    dev.off() 
    print(data$p6)
  }
})

output$table_fantom <- DT::renderDataTable({
  if (is.null(v$data_fa)) return()
  datatable(v$data_fa,extensions = 'Buttons', 
            options = list(dom = 'Bfrtip',
                           buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))  %>% 
    formatRound(columns=c( "medianratio_day12overallfantom" , "medianratio_matureoverday12" ,  "medianratio_matureoverallfantom"), digits=3)
})