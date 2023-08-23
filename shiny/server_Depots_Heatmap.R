
output$ref_depot_ht <- renderUI(
  selectInput(inputId = 'ref_depot_ht',label = 'Reference depot.',
              choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
              selected = 'SAT Abdomen',multiple = F))
output$qry_depot_ht <- renderUI(
  selectInput(inputId = 'qry_depot_ht',label = 'Query depot.',
              choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
              selected = 'VAT Omentum',multiple = F))

observeEvent(input$ref_depot_ht, {
  if (input$ref_depot_ht == input$qry_depot_ht) {
    updateSelectInput(session = session, inputId = 'qry_depot_ht',label = 'Query depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$ref_depot_ht,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})
observeEvent(input$qry_depot_ht, {
  if (input$ref_depot_ht == input$qry_depot_ht) {
    updateSelectInput(session = session, inputId = 'ref_depot_ht',label = 'Reference depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$qry_depot_ht,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})



DE_holder_ht <- reactiveValues(H1=NULL,H2=NULL)

DiffExp_ht <- eventReactive(input$SearchButton_ht, {
  
  
  genes <- unlist(strsplit(input$gene_lst_ht,split=', |,| |\n|\t'))
  
  Mat.r.temp <- matrix(0, nrow = length(genes), ncol = length(input$cohort_ht))
  Mat.P.temp <- matrix(1, nrow = length(genes), ncol = length(input$cohort_ht))
  rownames(Mat.r.temp) <- genes; rownames(Mat.P.temp) <- genes
  colnames(Mat.r.temp) <- input$cohort_ht; colnames(Mat.P.temp) <- input$cohort_ht
  
  counter <- 1
  for (i in input$cohort_ht) {
    Eset <- readRDS(paste0('./data/Depots/',i,'.rds'))
    
    Eset.pData <- pData(Eset)
    tissues <- Eset.pData$Tissue
    
    keep <- switch('Healthy' %in% input$disease_status_ht, T=double(length = length(Eset@phenoData@data[['ID']])) + 1, F=double(length = length(Eset@phenoData@data[['ID']])))
    for (j in list('Cancer','Obese')) {
      if (is.null(Eset@phenoData@data[[j]])) {
        next
      } else {
        if (j %in% input$disease_status_ht) {
          keep[Eset@phenoData@data[[j]] != 0] <- 1
        } else {
          keep[Eset@phenoData@data[[j]] != 0] <- 0
        }
      }
    }
    
    keep <- keep > 0
    
    depots_regex <- paste0(input$ref_depot_ht,'|',input$qry_depot_ht)
    depot_nums <- gsub(paste0(input$ref_depot_ht,'.*'), 1, tissues)
    depot_nums <- gsub(paste0(input$qry_depot_ht,'.*'), 2, depot_nums)
    
    gene_indices <- match(genes, Eset@featureData@data[['Genes']])
    depot_indices <- which(grepl(depots_regex, tissues) & keep)
    
    
    #if (length(depot_indices) < 4) { next }
    
    
    temp_1 <- depot_nums[depot_indices]
    temp_2 <- exprs(Eset)[gene_indices, depot_indices]
    
    temp_rcorr <- NULL
    tryCatch({temp_rcorr <- apply(temp_2, 1, function(x) {Hmisc::rcorr(x,temp_1,type = input$cor_method_ht)})},
             error = function(cond) {  })
    if (is.null(temp_rcorr)) { next }
    Mat.r.temp[, grep(i,colnames(Mat.r.temp))] <- unlist(lapply(temp_rcorr, function(x) {x$r[1,2]}))
    Mat.P.temp[, grep(i,colnames(Mat.P.temp))] <- unlist(lapply(temp_rcorr, function(x) {x$P[1,2]}))
    
  }
  
  list(Mat.r.temp, Mat.P.temp)
  
})

observeEvent(input$SearchButton_ht, {temp <- DiffExp_ht(); DE_holder_ht$H1 <- temp[[1]]; DE_holder_ht$H2 <- temp[[2]]})

H1 <- eventReactive(input$SearchButton_ht, {
  
  breaks_lst_r <- seq(-1,1,by=0.001)
  hmcols_r <- colorRampPalette(brewer.prgn(12))(length(breaks_lst_r) - 1)
  
  pheatmap(DE_holder_ht$H1, color=hmcols_r, breaks = breaks_lst_r,
           cluster_rows = T, cluster_cols = T, main = 'Correlation',
           treeheight_col = 0, treeheight_row = 0, border_color = NA)
  
})

H2 <- eventReactive(input$SearchButton_ht, {
  
  breaks_lst_P <- seq(0,1,by=0.001)
  hmcols_P <- colorRampPalette(rev(brewer.greys(12)), bias=3)(length(breaks_lst_P) - 1)
  
  temp_row <- H1()$tree_row$order
  temp_col <- H1()$tree_col$order
  
  pheatmap(DE_holder_ht$H2[temp_row,temp_col], 
           color=hmcols_P, breaks = breaks_lst_P,
           cluster_rows = F, cluster_cols = F, main = 'Significance',
           legend_breaks = c(0,0.1,1), border_color = NA)
  
})

#matrix(c(1,2,3,4,5,6,7,8,9,10),nrow=5,ncol=2)



output$H1_heatmap_ht <- renderPlot({ temp <- H1(); dev.off(); print(temp) })
output$H2_heatmap_ht <- renderPlot({ temp <- H2(); dev.off(); print(temp) })


output$download_H1_pdf <- downloadHandler(
  filename = function() {paste0('heatmap-corr_',Sys.Date(),'.pdf')},
  content = function(file) {ggsave(file,H1())},
  contentType = 'application/pdf')
output$download_H2_pdf <- downloadHandler(
  filename = function() {paste0('heatmap-P_',Sys.Date(),'.pdf')},
  content = function(file) {ggsave(file,H2())},
  contentType = 'application/pdf')
