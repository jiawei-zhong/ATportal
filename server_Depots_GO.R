output$cohort_go <- renderUI(
  selectInput(inputId = 'cohort_go',label = 'Select cohort to visualize.',
              choices = list('EMIF'='EMIF','Krieg et al.'='Krieg','Schleinitz et al.'='Schleinitz',
                             'Keller et al'='Keller','Mazaki et al'='Mazaki',
                             'Hoggard et al'='Hoggard'),
              selected = list('EMIF'='EMIF'),multiple = F))
output$ref_depot_go <- renderUI(
  selectInput(inputId = 'ref_depot_go',label = 'Reference depot.',
              choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
              selected = 'SAT Abdomen',multiple = F))
output$qry_depot_go <- renderUI(
  selectInput(inputId = 'qry_depot_go',label = 'Query depot.',
              choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
              selected = 'VAT Omentum',multiple = F))


observeEvent(input$ref_depot_go, {
  if (input$ref_depot_go == input$qry_depot_go) {
    updateSelectInput(session = session, inputId = 'qry_depot_go',label = 'Query depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$ref_depot_go,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})
observeEvent(input$qry_depot_go, {
  if (input$ref_depot_go == input$qry_depot_go) {
    updateSelectInput(session = session, inputId = 'ref_depot_go',label = 'Reference depot.',
                      choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                      selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$qry_depot_go,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])}})


depot_factor_go <- reactiveValues( var=1 )

observeEvent(input$ref_depot_go, {depot_factor_go$var <- switch(input$ref_depot_go,'SAT Abdomen'=1,'VAT Omentum'=-1)})

GO_term <- reactiveValues(var=NULL)

GO_DiffExp <- eventReactive(input$SearchButton_go, {
  
  temp <- switch(input$cor_method_go,'pearson'=data_GO_NES[[1]],'spearman'=data_GO_NES[[2]]) * depot_factor_go$var
  temp2 <- switch(input$cor_method_go,'pearson'=data_GO_P[[1]],'spearman'=data_GO_P[[2]])
  temp3 <- switch(input$cor_method_go,'pearson'=data_GO_Size[[1]],'spearman'=data_GO_Size[[2]])
  
  df_temp <- data.frame(pathway = gsub('_',' ',gsub('GOBP_','',rownames(temp))),
                        NES = temp[, grep(input$cohort_go,colnames(temp))], 
                        Enrichment = ifelse(temp[, grep(input$cohort_go,colnames(temp))] >= 0, "Up-regulated", "Down-regulated"),
                        adj.P_val = p.adjust(temp2[, grep(input$cohort_go,colnames(temp2))],method = input$p.adjust_go),
                        size = temp3[, grep(input$cohort_go,colnames(temp))])
  df_temp <- df_temp[complete.cases(df_temp), ] %>% dplyr::arrange(.,desc(NES))
  df_temp$Enrichment[df_temp$adj.P_val >= 0.05] <- '-'
  
  df_temp
  
})

observeEvent(input$SearchButton_go, {GO_term$var <- GO_DiffExp()})
observeEvent(input$cor_method_go, {GO_term$var <- GO_DiffExp()})
observeEvent(input$cohort_go, {GO_term$var <- GO_DiffExp()})


DE_GO <- reactive({ GO_term$var$adj.P_val < input$P_threshold_go })


DE_GO_DT <- reactiveValues( var=NULL )

observeEvent(input$SearchButton_go, {
  if (!is.null(GO_term$var)) {
    if (input$show_DE_go) {
      if (sum(DE_GO()) == 0) {
        temp <- NULL
      } else {
        temp <- dplyr::filter(GO_term$var, DE_GO()) #%>% dplyr::arrange(.,desc(NES))
      }
    } else {
      temp <- GO_term$var #dplyr::arrange(GO_term$var,desc(NES))
    }
    DE_GO_DT$var <- data.frame(Pathway=temp$pathway,NES=round(as.numeric(temp$NES),2),
                               adj.P_value=round(as.numeric(temp$adj.P_val),5),
                               size=temp$size)
  }
})

observeEvent(input$show_DE_go, {
  if (!is.null(GO_term$var)) {
    if (input$show_DE_go) {
      temp <- dplyr::filter(GO_term$var, DE_GO()) #%>% dplyr::arrange(.,desc(NES))
    } else {
      temp <- GO_term$var #dplyr::arrange(GO_term$var,desc(NES))
    }
    DE_GO_DT$var <- data.frame(Pathway=temp$pathway,NES=round(as.numeric(temp$NES),2),
                               adj.P_value=round(as.numeric(temp$adj.P_val),5),
                               size=temp$size)
  }
})

output$DT_go <- renderDataTable(rbind(head(DE_GO_DT$var, n = 2),tail(DE_GO_DT$var, n = 2 )))
output$tab_DT_go <- renderDataTable(DE_GO_DT$var)


enrichment_plot <- eventReactive(input$SearchButton_go, {
  
  colos <- setNames(c("darkgreen", "purple","grey"),c("Up-regulated", "Down-regulated","-"))
  GO_term_temp <- rbind(head(GO_term$var, n = 10),tail(GO_term$var, n = 10 ))
  ggplot(GO_term_temp, aes(reorder(GO_term_temp$pathway,GO_term_temp$NES), GO_term_temp$NES))+
    geom_point(alpha=0.6, aes(fill = GO_term_temp$Enrichment, size = abs(as.numeric(GO_term_temp$size))), shape=21)+
    guides(size=guide_legend(title='Pathway Size'),fill=guide_legend(title='Enrichment'))+
    scale_fill_manual(values = colos )+
    scale_size_continuous(range = c(4,10))+
    geom_hline(yintercept = 0)+
    coord_flip()+
    labs(x="", y="Normalized Enrichment Score")+ 
    theme_bw()+
    ylim(c(GO_term_temp$NES[length(GO_term_temp$NES)]-0.5, GO_term_temp$NES[1]+0.5))
})

output$graph_go <- renderPlot({
  
  enrichment_plot()
  
})


output$download_GO_enrch <- downloadHandler(
  filename = function() {paste0('GO_enriched_',Sys.Date(),'.pdf')},
  content = function(file) {ggsave(file,enrichment_plot())},
  contentType = 'application/pdf')
output$download_DT_go_csv <- downloadHandler(
  filename = function() {paste0('GO data table_',Sys.Date(),'.csv')},
  content = function(file) {write.csv(DE_GO_DT$var,file)},
  contentType = 'application/csv')
output$download_DT_go_xlsx <- downloadHandler(
  filename = function() {paste0('GO data table_',Sys.Date(),'.xlsx')},
  content = function(file) {write_xlsx(DE_GO_DT$var,file)},
  contentType = 'application/xlsx')