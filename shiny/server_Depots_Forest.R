
plot_temp_fr <- eventReactive(input$SearchBotton_fr, {
  
  if (input$ref_depot_fr == input$qry_depot_fr) {
    stop('Reference depot is the same as the query depot.')
  }
  
  cor <- list()
  for (i in 1:length(input$cohort_fr)) {
    Eset <- readRDS(paste0('./data/Depots/',input$cohort_fr[i],'.rds'))
    
    genes <- Eset@featureData@data[['Genes']]
    
    Eset.pData <- pData(Eset)
    tissues <- Eset.pData$Tissue
    
    if (!input$gene_fr %in% genes) {
      next
    } else if (!input$ref_depot_fr %in% tissues | !input$qry_depot_fr %in% tissues) {
      next
    }
    
    keep <- switch('Healthy' %in% input$disease_status_fr,T=double(length = length(Eset@phenoData@data[['ID']])) + 1, F=double(length = length(Eset@phenoData@data[['ID']])))
    for (j in list('Cancer','Obese')) {
      if (is.null(Eset@phenoData@data[[j]])) {
        next
      } else {
        if (j %in% input$disease_status_fr) {
          keep[Eset@phenoData@data[[j]] != 0] <- 1
        } else {
          keep[Eset@phenoData@data[[j]] != 0] <- 0
        }
      }
    }
    
    keep <- keep > 0
    
    depots_regex <- paste0(input$ref_depot_fr,'|',input$qry_depot_fr)
    depot_nums <- gsub(paste0(input$ref_depot_fr,'.*'), 1, tissues)
    depot_nums <- gsub(paste0(input$qry_depot_fr,'.*'), 2, depot_nums)
    
    
    gene_indices <- which(genes == input$gene_fr)
    depot_indices <- which(grepl(depots_regex, tissues) & keep)
    
    if (length(depot_indices) < 4) {next}
    
    temp_1 <- depot_nums[depot_indices]
    temp_2 <- exprs(Eset)[gene_indices, depot_indices]
    
    temp <- Hmisc::rcorr(temp_2,temp_1,type = input$cor_method_fr)
    cor$r <- c(cor$r,temp$r[1,2])
    cor$n <- c(cor$n,temp$n[1,2])
    cor$ID <- c(cor$ID,input$cohort_fr[i])
    
  }
  
  meta::metacor(cor = cor$r,n = cor$n, studlab = cor$ID)
  
})

forest <- reactive({meta::forest(plot_temp_fr(), test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4)})

output$graph_fr <- renderPlot({
  
  forest()

})

output$download_forest_pdf <- downloadHandler(
  
  filename = function() {paste0('forest-plot',Sys.Date(),'.pdf')},
  content = function(file) {pdf(file,width=10,height=6); meta::forest(plot_temp_fr(), test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4); dev.off()},
  contentType = 'application/pdf'
  
)
