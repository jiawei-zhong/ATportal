
output$cohort_vc <- renderUI(
  selectInput(inputId = 'cohort_vc',label = 'Select cohort to visualize.',
              choices = list('EMIF'='EMIF','Krieg et al.'='Krieg','Schleinitz et al.'='Schleinitz',
                             'Keller et al'='Keller','Mazaki et al'='Mazaki',
                             'Hoggard et al'='Hoggard'),
              selected = list('EMIF'='EMIF'),multiple = F))
output$ref_depot_vc <- renderUI(
  selectInput(inputId = 'ref_depot_vc',label = 'Reference depot.',
              choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
              selected = 'SAT Abdomen',multiple = F))
output$qry_depot_vc <- renderUI(
  selectInput(inputId = 'qry_depot_vc',label = 'Query depot.',
              choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
              selected = 'VAT Omentum',multiple = F))

observeEvent(input$ref_depot_vc, {
  if (input$ref_depot_vc == input$qry_depot_vc) {
    updateSelectInput(session = session, inputId = 'qry_depot_vc',label = 'Query depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$ref_depot_vc,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})
observeEvent(input$qry_depot_vc, {
  if (input$ref_depot_vc == input$qry_depot_vc) {
    updateSelectInput(session = session, inputId = 'ref_depot_vc',label = 'Reference depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$qry_depot_vc,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})

output$FC_threshold_vc <- renderUI(
  sliderInput(inputId = 'logFC_threshold_vc',label = 'Absolute fold change cut-off.',
              min = 0, max = 10, value = 1, step = 0.05))
output$P_value_threshold_vc <- renderUI(
  sliderInput(inputId = 'P_threshold_vc',label = 'P-value cut-off.',
              min = 0, max = 1, value = 0.05, step = 0.01))


DE_holder_vc <- reactiveValues(var=NULL)

DiffExp_vc <- eventReactive(input$SearchButton_vc, {
  
  Eset <- readRDS(paste0('./data/Depots/',input$cohort_vc,'.rds'))
  genes <- Eset@featureData@data[['Genes']]
  tissues <- Eset@phenoData@data[['Tissue']]
  
  keep <- switch('Healthy' %in% input$disease_status_vc,T=double(length = length(Eset@phenoData@data[['ID']])) + 1, F=double(length = length(Eset@phenoData@data[['ID']])))
  for (j in list('Cancer','Obese')) {
    if (is.null(Eset@phenoData@data[[j]])) {
      next
    } else {
      if (j %in% input$disease_status_vc) {
        keep[Eset@phenoData@data[[j]] != 0] <- 1
      } else {
        keep[Eset@phenoData@data[[j]] != 0] <- 0
      }
    }
  }
  
  keep <- keep > 0
  
  temp <- double(length = length(tissues))
  temp[grepl(paste0('^',input$ref_depot_vc,'$'),tissues) & keep] <- input$ref_depot_vc
  temp[grepl(paste0('^',input$qry_depot_vc,'$'),tissues) & keep] <- input$qry_depot_vc
  
  if (sum(temp != 0) < 3) { return(NULL) } 
  
  test <- mod.t.test(exprs(Eset)[,grep(paste0(input$ref_depot_vc,'|',input$qry_depot_vc),tissues)],group=factor(temp[grep(paste0(input$ref_depot_vc,'|',input$qry_depot_vc),tissues)]),adjust.method=input$p.adjust_vc)
  logFC <- test[, grep(input$qry_depot_vc,colnames(test))] - test[, grep(input$ref_depot_vc,colnames(test))]
  
  data.frame(z = genes, x = logFC, y = test$adj.p.value, 
             a = test[, grep(input$ref_depot_vc,colnames(test))],
             b = test[, grep(input$qry_depot_vc,colnames(test))])
  
})


observeEvent(input$SearchButton_vc, {DE_holder_vc$var <- DiffExp_vc()})
observeEvent(input$p.adjust_vc, {DE_holder_vc$var <- DiffExp_vc()})
observeEvent(input$ref_depot_vc, {DE_holder_vc$var <- DiffExp_vc()})

DE_genes_vc <- reactive({(abs(as.numeric(DE_holder_vc$var$x)) >= input$logFC_threshold_vc) & (DE_holder_vc$var$y <= input$P_threshold_vc)})

DE_genes_vc_DT <- reactiveValues(var=NULL)

observeEvent(input$SearchButton_vc, {
  if (!is.null(DE_holder_vc$var)) {
    if (input$show_DE_vc) {
      if (sum(DE_genes_vc()) == 0) {
        temp <- NULL
      } else {
        temp <- dplyr::filter(DE_holder_vc$var, DE_genes_vc()) %>% dplyr::arrange(.,desc(x))
      }
    } else {
      temp <- dplyr::arrange(DE_holder_vc$var,desc(x))
    }
    DE_genes_vc_DT$var <- data.frame(row.names=temp$z,Ref_Depot_Norm=round(temp$a,2),
                                  Qry_Depot_Norm=round(temp$b,2),log2FC=round(temp$x,2),
                                  adj.P_value=round(temp$y,5))
  } 
})

observeEvent(input$show_DE_vc, {
  if (!is.null(DE_holder_vc$var)) {
    if (input$show_DE_vc) {
      temp <- dplyr::filter(DE_holder_vc$var, DE_genes_vc()) %>% dplyr::arrange(.,desc(x))
    } else {
      temp <- dplyr::arrange(DE_holder_vc$var,desc(x))
    }
    DE_genes_vc_DT$var <- data.frame(row.names=temp$z,Ref_Depot=round(temp$a,2),
                                  Qry_Depot=round(temp$b,2),log2FC=round(temp$x,2),
                                  adj.P_value=round(temp$y,5))
  } 
})

output$tab_DT_vc <- renderDataTable(DE_genes_vc_DT$var)


vline <- function(x = 0, color = "grey") {
  return(list(type="line",y0=0,y1=1,yref="paper",x0=x,x1=x,line=list(color=color,dash='dot')))
}

hline <- function(y = 0, color = "grey") {
  return(list(type="line",x0=0,x1=1,xref="paper",y0=y,y1=y,line=list(color=color,dash='dot')))
}

volcano <- reactive({
  
  plot_ly(type='scatter',mode='markers',color = DE_genes_vc(), colors=c('grey70','forestgreen')) %>% 
    add_trace(x = DE_holder_vc$var$x, y = -log10(as.numeric(DE_holder_vc$var$y)), 
              text= DE_holder_vc$var$z,alpha=0.6, 
              hovertemplate=paste('<b>%{text}</b>','<br>log2 FC: %{x:.2f}','<br>-log10 FDR: %{y:.2f}','<extra></extra>')) %>%
    layout(shapes=list(hline(input$P_threshold_vc),vline(-input$logFC_threshold_vc),vline(input$logFC_threshold_vc)),
           legend = list(title = list(text = "<b>DE</b>")),
           xaxis = list(title = "log2 FC"), yaxis = list(title = "-log10 P.adj"))

})

output$graph_vc <- renderPlotly({ volcano() })


output$download_volcano_pdf <- downloadHandler(
  filename = function() {paste0('volcano-plot_',Sys.Date(),'.pdf')},
  content = function(file) {export(volcano(), file=file)},
  contentType = 'application/pdf')
output$download_DT_vc_csv <- downloadHandler(
  filename = function() {paste0('volcano-plot data table_',Sys.Date(),'.csv')},
  content = function(file) {write.csv(DE_genes_vc_DT$var,file)},
  contentType = 'application/csv')
output$download_DT_vc_xlsx <- downloadHandler(
  filename = function() {paste0('volcano-plot data table_',Sys.Date(),'.xlsx')},
  content = function(file) {write_xlsx(DE_genes_vc_DT$var,file)},
  contentType = 'application/xlsx')
