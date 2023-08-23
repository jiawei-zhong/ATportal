
observeEvent(input$ref_depot_bx, {
  if (input$ref_depot_bx == input$qry_depot_bx) {
    updateSelectInput(session = session, inputId = 'qry_depot_bx',label = 'Query depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$ref_depot_bx,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})
observeEvent(input$qry_depot_bx, {
  if (input$ref_depot_bx == input$qry_depot_bx) {
    updateSelectInput(session = session, inputId = 'ref_depot_bx',label = 'Reference depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$qry_depot_bx,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})


plot_temp_bx <- eventReactive(input$SearchButton_bx, {
  
  P_values <- list()
  for (i in input$cohort_bx) {
    
    Eset <- readRDS(paste0('./data/Depots/',i,'.rds'))
    genes <- Eset@featureData@data[['Genes']]
    
    if (i == 'Schleinitz') {Eset@phenoData@data[['ID']] <- gsub('DS |DS_|_.*|-.*','',Eset@phenoData@data[['ID']])}
    
    Eset.pData <- pData(Eset)
    tissues <- Eset.pData$Tissue
    
    temp <- switch('Healthy' %in% input$disease_status_bx,T=double(length = length(Eset@phenoData@data[['ID']])) + 1, F=double(length = length(Eset@phenoData@data[['ID']])))
    for (j in list('Cancer','Obese')) {
      if (is.null(Eset@phenoData@data[[j]])) {
        next
      } else {
        if (j %in% input$disease_status_bx) {
          temp[Eset@phenoData@data[[j]] != 0] <- 1
        } else {
          temp[Eset@phenoData@data[[j]] != 0] <- 0
        }
      }
    }
    
    temp <- temp > 0
    
    gene_index <- which(genes == input$gene_bx)
    
    temp_1 <- which(grepl(input$qry_depot_bx,tissues) & temp)
    temp_2 <- which(grepl(input$ref_depot_bx,tissues) & temp)
    
    keep <- c(Eset.pData$ID[temp_1],Eset.pData$ID[temp_2])[duplicated(c(Eset.pData$ID[temp_1],Eset.pData$ID[temp_2]))]
    temp_1 <- temp_1[Eset.pData$ID[temp_1] %in% keep]
    temp_2 <- temp_2[Eset.pData$ID[temp_2] %in% keep]
    
    
    FC_iter <- exprs(Eset)[gene_index, temp_1[match(Eset.pData$ID[temp_1], Eset.pData$ID[temp_2])]] -
      exprs(Eset)[gene_index, temp_2]
    
    
    if (i == input$cohort_bx[1]) {
      df.iter <- data.frame(a = FC_iter, b = rep(i,length(FC_iter)))
    } else {
      temp <- data.frame(a = FC_iter, b = rep(i,length(FC_iter)))
      df.iter <- rbind(df.iter, temp)
    }
  }
  
  ggplot(df.iter, aes(x = b, y = a, fill = b))+ 
    geom_boxplot(alpha=0.6)+ 
    theme_bw()+ ylab('log2 FC')+ xlab('')+ 
    guides(fill=guide_legend(title='Cohort'))
  
})

output$graph_bx <- renderPlotly({

  ggplotly(plot_temp_bx())
  
})

output$download_boxplot_pdf <- downloadHandler(
  filename = function() {paste0('boxplot-depot_',Sys.Date(),'.pdf')},
  content = function(file) {ggsave(file,plot_temp_bx())},
  contentType = 'application/pdf')
