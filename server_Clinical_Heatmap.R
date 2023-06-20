## heatmap



# observe({updatePrettyCheckboxGroup(session = session, 
#                                    inputId = 'clinical_measure_hm',  
#                                    choices = clinical_measure_vector,
#                                    selected = if (input$select_all_hm) {clinical_measure_vector},
#                                    inline = T )})

observeEvent(input$cohort_hm, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_hm,'.rds'))
  pData <- pData(eset)
  pData <- pData[,colSums(is.na(pData))<nrow(pData)]
  
  
  
  updatePickerInput(session = session,
                    inputId = 'clinical_measure_hm',
                    label = 'Clinical measures',
                    # selected = '',
                    choices = intersect(clinical_measure_vector,colnames(pData))
                    # width = "300px",
                    # options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                    # multiple = TRUE
  )
  
})

observeEvent(input$clinical_measure_hm, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_hm,'.rds'))
  pData <- pData(eset)
  pData <- pData[,colSums(is.na(pData))<nrow(pData)]
  
  
  
  updatePickerInput(session = session,
                    inputId = 'clinical_measure_adjust_hm',
                    label = 'Adjustment',
                    # selected = '',
                    choices = setdiff(intersect(clinical_measure_vector,colnames(pData)),input$clinical_measure_hm)
                    # width = "300px",
                    # options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                    # multiple = TRUE
  )
})




plotdata_hm <- eventReactive(input$SearchButton_hm, {
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_hm,'.rds'))
  exprs <- exprs(eset)
  if (nchar(input$gene_hm) > 0) {
    exprs <- exprs[rownames(exprs) %in% unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_hm),split='\n|\t|,')),]
  } else {
    exprs <- exprs[1:30,]
  }
  
  
  
  
  
  if (length(input$clinical_measure_adjust_hm)==0) {
    pData <- pData(eset)
    pData <- pData[,colnames(pData) %in% input$clinical_measure_hm]
    pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
    cor_df <- t(rbind(exprs,t(pData)))
    cor_df <- psych::corr.test(cor_df, method = input$correlation_method_hm, adjust="none")
    cor_df$r <- (cor_df$r)[rownames(exprs),colnames(pData)]
    cor_df$p <- (cor_df$p)[rownames(exprs),colnames(pData)]
    
    
    
  } else {
    
    pData <- pData(eset)
    pData <- pData[,colnames(pData) %in% c(input$clinical_measure_hm, input$clinical_measure_adjust_hm)]
    pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
    cor_df <- t(rbind(exprs,t(pData)))
    temp <- list()
    temp$r <- matrix(1, nrow = length(rownames(exprs)), ncol = length(input$clinical_measure_hm))
    rownames(temp$r) <- rownames(exprs)
    colnames(temp$r) <- input$clinical_measure_hm
    temp$p <- temp$r
    
    for (gene in rownames(exprs)) {
      for (measure in input$clinical_measure_hm) {
        cor_df_temp <- cor_df[,c(gene,measure,input$clinical_measure_adjust_hm)] %>% na.omit()
        cor_df_temp <- pcor.test(x = cor_df_temp[,gene], y = cor_df_temp[,measure], z = cor_df_temp[,input$clinical_measure_adjust_hm], method = input$correlation_method_hm)
        temp$r[gene,measure] <- cor_df_temp$estimate
        temp$p[gene,measure] <- cor_df_temp$p.value
        
      }
    }
    
    cor_df <- temp
    
  }
  
  cor_df$p <- (-log10(cor_df$p))
  
  # cor_df$n <- (cor_df$n)[rownames(exprs),colnames(pData)]
  
  row_order <- hclust(d = dist(cor_df$r))$order
  column_order <- hclust(d = dist(t(cor_df$r)))$order
  # cor_df$order_row <- row_order
  # cor_df$order_col <- column_order
  cor_df$r <- cor_df$r[row_order,column_order]
  cor_df$p <- cor_df$p[row_order,column_order]
  
  
  
  p1 <- pheatmap::pheatmap(mat = cor_df$r,
                           cluster_cols = T,
                           cluster_rows = T,
                           cellwidth = 20,
                           cellheight = 20,
                           border_color = NA,
                           treeheight_col = 0,
                           treeheight_row = 30, 
                           fontsize = 15,
                           main = "Heatmap of correlation",
                           silent = T,
                           show_colnames = T,
                           color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
  )
  
  p2 <- pheatmap::pheatmap(mat = cor_df$p,
                           cluster_cols = F,
                           cluster_rows = F,
                           display_numbers = ifelse(test = cor_df$p > (-log10(0.001)),
                                                    yes =  "***",
                                                    no =  ifelse(test = cor_df$p > (-log10(0.01)),
                                                                 yes =  "**",
                                                                 no =  ifelse(test = cor_df$p > (-log10(0.05)),
                                                                              yes =  "*",
                                                                              no =  ""))),
                           cellwidth = 20,
                           cellheight = 20,
                           border_color = NA,
                           treeheight_col = 0,
                           treeheight_row = 0,
                           fontsize = 15,
                           main = "Heatmap of P-value (-log10)",
                           silent = T,
                           show_colnames = T
  )
  
  
  cor_df$heatmap <- plot_grid(as.ggplot(p1), as.ggplot(p2),ncol=2)
  
  cor_df
})

output$heat_map <- renderPlot({
  
  plotdata_hm()$heatmap
  
})




output$r_hm <- DT::renderDataTable(server = FALSE,{return(round(plotdata_hm()$r,digits = 4))},
                                   extensions = c('Buttons'),
                                   options = list(scrollX = TRUE,
                                                  pageLength = 10,
                                                  lengthMenu = c(10, 25, 50, 100),
                                                  dom = 'Blfrtip',
                                                  buttons = list(
                                                    list(extend = "csv", text = "Download Current Page", filename = "page",
                                                         exportOptions = list(modifier = list(page = "current"))),
                                                    list(extend = "csv", text = "Download Full matrix", filename = "data",
                                                         exportOptions = list(modifier = list(page = "all"))))))
output$p_hm <- DT::renderDataTable(server = FALSE,{return(round(plotdata_hm()$p,digits = 4))},
                                   extensions = c('Buttons'),
                                   options = list(scrollX = TRUE,
                                                  pageLength = 10,
                                                  lengthMenu = c(10, 25, 50, 100),
                                                  dom = 'Blfrtip',
                                                  buttons = list(
                                                    list(extend = "csv", text = "Download Current Page", filename = "page",
                                                         exportOptions = list(modifier = list(page = "current"))),
                                                    list(extend = "csv", text = "Download Full matrix", filename = "data",
                                                         exportOptions = list(modifier = list(page = "all"))))))

# output$p <- DT::renderDataTable({return(plotdata()$p)},extensions = c('Buttons'),options = list(scrollX = TRUE))
# output$download_r <- downloadHandler(filename = function(){"matrix_r.rho.csv"}, 
#                                      content = function(fname){write.csv(plotdata()$r, fname, quote = F, row.names = T)})
# output$download_P <- downloadHandler(filename = function(){"matrix_P.csv"}, 
#                                      content = function(fname){write.csv(plotdata()$p, fname, quote = F, row.names = T)})

output$downloadPlot_hm <- downloadHandler(filename = function() {"heatmap.pdf"},
                                          content = function(fname) {ggsave(filename = fname, plot = plotdata_hm()$heatmap, height = 40, width = 20, units = "in")})



# output$ui_plot_hm <- renderUI({ if (input$SearchButton_hm) {
#   list(box(width=12,
#            style=paste0('height: ',min(30*20+120,nrow(plotdata_hm()$r)*20+120),'px; overflow-x: scroll; overflow-y: scroll;'),
#            plotOutput(outputId = "heat_map", height = paste0(nrow(plotdata_hm()$r)*20+100,'px'), width = paste0(ncol(plotdata_hm()$r)*20+900,'px'))),
#        column(12,downloadBttn(outputId = 'downloadPlot_hm', label = "Download heatmap", style= "simple", color = "primary", size = "sm"))
#   )}
# })




output$ui_all_hm <- renderUI({if (input$SearchButton_hm) {
  tabsetPanel(id = "inTabset",
              tabPanel(title = "plot",
                       column(12,list(box(width=12,
                                          style=paste0('height: ',min(30*20+120,nrow(plotdata_hm()$r)*20+120),'px; overflow-x: scroll; overflow-y: scroll;'),
                                          plotOutput(outputId = "heat_map", height = paste0(nrow(plotdata_hm()$r)*20+100,'px'), width = paste0(ncol(plotdata_hm()$r)*20+900,'px'))),
                                      column(12,downloadBttn(outputId = 'downloadPlot_hm', label = "Download heatmap", style= "simple", color = "primary", size = "sm"))
                       ))),
              tabPanel(title = "r/rho matrix",
                       column(12,DT::dataTableOutput("r_hm"))),
              tabPanel(title = "P-value matrix (-log10)", 
                       column(12,DT::dataTableOutput("p_hm"))))
}})

# output$ui_table <- renderUI({if (input$SearchButton) {
#   tabsetPanel(id = "inTabset",
#               tabPanel(title = "r/rho matrix",
#                        column(12,DT::dataTableOutput("r")),
#                        column(12,downloadButton('download_r',"Download r/rho matrix"))),
#               tabPanel(title = "P-value matrix (-log10)", 
#                        column(12,DT::dataTableOutput("P")),
#                        column(12,downloadButton('download_P',"Download P-value matrix"))))
# }})