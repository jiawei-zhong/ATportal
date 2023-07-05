## heatmap



# observe({updatePrettyCheckboxGroup(session = session, 
#                                    inputId = 'clinical_measure_sp',  
#                                    choices = clinical_measure_vector,
#                                    selected = if (input$select_all_sp) {clinical_measure_vector},
#                                    inline = T )})

observeEvent(input$cohort_sp, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_sp,'.rds'))
  pData <- pData(eset)
  pData <- pData[,colSums(is.na(pData))<nrow(pData)]
  
  
  
  updatePickerInput(session = session,
                    inputId = 'clinical_measure_sp',
                    label = 'Clinical measures',
                    # selected = '',
                    choices = intersect(clinical_measure_vector,colnames(pData))
                    # width = "300px",
                    # options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                    # multiple = TRUE
  )
  
})

observeEvent(input$clinical_measure_sp, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_sp,'.rds'))
  pData <- pData(eset)
  pData <- pData[,colSums(is.na(pData))<nrow(pData)]
  
  
  
  updatePickerInput(session = session,
                    inputId = 'clinical_measure_adjust_sp',
                    label = 'Adjustment',
                    # selected = '',
                    choices = setdiff(intersect(clinical_measure_vector,colnames(pData)),input$clinical_measure_sp)
                    # width = "300px",
                    # options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                    # multiple = TRUE
  )
})




plotdata_sp <- eventReactive(input$SearchButton_sp, {

  
  if (length(input$clinical_measure_adjust_sp)==0) {
    
    eset <- readRDS(paste0('./data/Clinical/',input$cohort_sp,'.rds'))
    exprs <- exprs(eset)
    if (nchar(input$gene_sp) > 0) {
      exprs <- exprs[rownames(exprs) %in% unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_sp),split='\n|\t|,')),]
    }
    
    pData <- pData(eset)
    pData <- pData[,colnames(pData) %in% input$clinical_measure_sp,drop=FALSE]
    # pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
    # cor_df <- t(rbind(exprs,t(pData)))
    # cor_df <- psych::corr.test(cor_df, method = input$correlation_method_sp, adjust="none")
    # cor_df$r <- (cor_df$r)[rownames(exprs),colnames(pData),drop=FALSE]
    # cor_df$p <- (cor_df$p)[rownames(exprs),colnames(pData),drop=FALSE]
    
    cor_df <- psych::corr.test(x = t(exprs),y = pData, method = input$correlation_method_sp, adjust="none")
    cor_df <- data.frame(r=cor_df$r[,1],p=cor_df$p[,1],gene=rownames(cor_df$r),mean_exp=rowMeans(exprs))
    
    cor_df <- cor_df[order(cor_df$r,decreasing = T),]
    cor_df$gene <- factor(cor_df$gene,levels = cor_df$gene)
    # cor_df <- cor_df[c(1:30,(nrow(cor_df)-29):nrow(cor_df)),]
    
  } else {
    
    eset <- readRDS(paste0('./data/Clinical/',input$cohort_sp,'.rds'))
    exprs <- exprs(eset)
    if (nchar(input$gene_sp) > 0) {
      exprs <- exprs[rownames(exprs) %in% unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_sp),split='\n|\t|,')),]
    } else {
      exprs <- exprs[1:30,]
    }
    
    pData <- pData(eset)
    pData <- pData[,colnames(pData) %in% c(input$clinical_measure_sp, input$clinical_measure_adjust_sp)]
    pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
    cor_df <- t(rbind(exprs,t(pData)))
    temp <- list()
    temp$r <- matrix(1, nrow = length(rownames(exprs)), ncol = length(input$clinical_measure_sp))
    rownames(temp$r) <- rownames(exprs)
    colnames(temp$r) <- input$clinical_measure_sp
    temp$p <- temp$r
    
    for (gene in rownames(exprs)) {
      for (measure in input$clinical_measure_sp) {
        cor_df_temp <- cor_df[,c(gene,measure,input$clinical_measure_adjust_sp)] %>% na.omit()
        cor_df_temp <- pcor.test(x = cor_df_temp[,gene], y = cor_df_temp[,measure], z = cor_df_temp[,input$clinical_measure_adjust_sp], method = input$correlation_method_sp)
        temp$r[gene,measure] <- cor_df_temp$estimate
        temp$p[gene,measure] <- cor_df_temp$p.value
        
      }
    }
    
    cor_df <- temp
    
    cor_df <- data.frame(r=cor_df$r[,1],p=cor_df$p[,1],gene=rownames(cor_df$r),mean_exp=rowMeans(exprs))
    
    cor_df <- cor_df[order(cor_df$r,decreasing = T),]
    cor_df$gene <- factor(cor_df$gene,levels = cor_df$gene)
  }
  
  cor_df$p <- (-log10(cor_df$p))
  
  cor_df_temp <- cor_df[abs(cor_df$r)>=as.numeric(input$correlation_cutoff_sp),]
  
  if (nrow(cor_df_temp)>60) {
    cor_df_temp <- cor_df_temp[c(1:30,(nrow(cor_df_temp)-29):nrow(cor_df_temp)),]
  } else {
    cor_df_temp <- cor_df_temp
  }
  
  temp <- list()
  
  temp$df <- cor_df
  
  temp$plot <- ggscatter(cor_df_temp, x = "gene", y = "r", color = "p", size = "mean_exp")+
    theme(legend.position = "right",
          axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) + 
    labs(x = "Genes", y = paste0(ifelse(test = input$correlation_method_sp=="spearman",yes = "Spearman",no = "Pearson")," correlation"), color = "- Log10 P-value", size = "Mean expression") + 
    scale_color_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[26:51])
  
  
  temp
  

  
})

output$scatter <- renderPlot({
  
  plotdata_sp()$plot
  
})


output$df_sp <- DT::renderDataTable(server = FALSE,{return(plotdata_sp()$df %>% mutate_if(is.numeric, round,4))},
                                    extensions = c('Buttons'),
                                    options = list(scrollX = TRUE,
                                                   pageLength = 10,
                                                   lengthMenu = c(10, 25, 50, 100),
                                                   dom = 'Blfrtip',
                                                   buttons = list(
                                                     list(extend = "csv", text = "Download Current Page", filename = "page",
                                                          exportOptions = list(modifier = list(page = "current"))),
                                                     list(extend = "csv", text = "Download Full data frame", filename = "data",
                                                          exportOptions = list(modifier = list(page = "all"))))))



output$downloadPlot_sp <- downloadHandler(filename = function() {"scatter.pdf"},
                                          content = function(fname) {ggsave(filename = fname, plot = plotdata_sp()$plot, height = 10, width = 20, units = "in")})




output$ui_sp <- renderUI({if (input$SearchButton_sp) {
  tabsetPanel(id = "inTabset_sp",
              tabPanel(title = "Plot",
                       column(12,plotOutput(outputId = "scatter", height = '600px')),
                       if (nrow(plotdata_sp()$plot$data)==60) {
                         column(12,paste0('please note that only top 60 genes were visualized'))
                       },
                       column(12,downloadBttn(outputId = 'downloadPlot_sp', label = "Download scatter plot", style= "simple", color = "primary", size = "sm"))),
              tabPanel(title = "Data frame",
                       column(12,DT::dataTableOutput("df_sp"))))
  }})






