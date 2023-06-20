## HOMA boxplot


plotdata_box_HOMA <- eventReactive(input$SearchButton_box_HOMA, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_box_HOMA,'.rds'))
  pData <- pData(eset)
  pData <- pData[,colSums(is.na(pData))<nrow(pData)]
  pData$expression <- exprs(eset)[input$gene_HOMA,]
  pData <- pData[!is.na(pData$HOMA),]
  pData <- pData[order(pData$HOMA,decreasing = F),]
  pData$HOMA_catelogy <- "1"
  quantile <- quantile(pData$HOMA,probs = seq(0,1,1/as.numeric(input$quantile_HOMA)))
  for (i in 1:as.numeric(input$quantile_HOMA )) {
    pData$HOMA_catelogy[pData$HOMA >= quantile[i] & pData$HOMA <= quantile[i+1]] <- paste0(signif(quantile[i],digits = 3),'~',signif(quantile[i+1],digits = 3))
  }
  
  
  compare_list <- list()
  
  for (i in 1:(length(unique(pData$HOMA_catelogy))-1)) {
    for (j in (i+1):length(unique(pData$HOMA_catelogy))) {
      compare_list <- append(compare_list,list(unique(pData$HOMA_catelogy)[c(i,j)]))
    }
  }
  
  
  ggboxplot(data = pData,x = "HOMA_catelogy",y = "expression",fill = "HOMA_catelogy",palette = "ucscgb",xlab = "HOMA catelogies",ylab = paste0(input$gene_HOMA, " expression"), title = input$cohort_box_HOMA)+
    NoLegend() +
    theme(aspect.ratio=1,
          axis.text.x = element_text(face = "plain",size = 11,angle = 45, hjust = 1, vjust = 1)) +
    stat_compare_means(comparisons = compare_list, paired = F, method = input$statistical_test_cor_HOMA)
  
})

output$boxplot_HOMA <- renderPlot({
  
  plotdata_box_HOMA()
  
})


output$downloadPlot_box_HOMA <- downloadHandler(filename = function() {"boxplot.pdf"},
                                                content = function(fname) {
                                                  ggsave(filename = fname, plot = plotdata_box_HOMA(), height = 10, width = 10, units = "in")})

output$ui_boxplot_HOMA <- renderUI({ if (input$SearchButton_box_HOMA) {
  
  list(column(12,plotOutput(outputId = "boxplot_HOMA", height = '600px')),
       column(12,downloadBttn(outputId = 'downloadPlot_box_HOMA', label = "Download boxplot", style= "simple", color = "primary", size = "sm")))
  
}
})