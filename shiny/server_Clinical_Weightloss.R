## WL boxplot

plotdata_wl <- eventReactive(input$SearchButton_wl, {
  
  DEOSH_WL <- readRDS("./data/Clinical/DEOSHeset_WL.rds")
  temp_pdat <- pData(DEOSH_WL)
  temp_pdat$expression <- exprs(DEOSH_WL)[input$gene_wl,]
  
  # DEOSH_baseline <- readRDS("./data/Clinical/DEOSHeset_Baseline.rds")
  # pData <- pData(DEOSH_baseline)
  # pData <- pData[pData$DNA %in% temp_pdat$subject,]
  # pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
  # colnames(pData)[1] <- 'subject'
  # temp_pdat <- merge(temp_pdat,pData,by="subject",all=T)
  
  p1 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",add = "point",fill = "time_point",palette = "ucscgb",xlab = "years after bariatric surgery",ylab = paste0(input$gene_wl, " expression"), title = "DEOSH weight loss")+
    geom_line(aes(group=subject), linetype = "dashed") + NoLegend() + theme(aspect.ratio=1) +
    stat_compare_means(comparisons = list(c("0","2"),c("2","5"),c("0","5")), paired = TRUE, method = input$statistical_test_wl)
  
  
  PO_WL <- readRDS("./data/Clinical/POeset_WL.rds")
  temp_pdat <- pData(PO_WL)
  temp_pdat$expression <- exprs(PO_WL)[input$gene_wl,]
  p2 <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",add = "point",fill = "time_point",palette = "ucscgb",xlab = "years after bariatric surgery",ylab = paste0(input$gene_wl, " expression"), title = "PO weight loss") +
    stat_compare_means(comparisons = list(c("0","2")), paired = TRUE, method = input$statistical_test_wl) +
    geom_line(aes(group=subject), linetype = "dashed") + NoLegend() + theme(aspect.ratio=1)
  
  p1+p2
})

output$boxplot_wl <- renderPlot({
  
  plotdata_wl()
  
})



output$downloadPlot_wl <- downloadHandler(filename = function() {"boxplot.pdf"},
                                          content = function(fname) {
                                            ggsave(filename = fname, plot = plotdata_wl(), height = 10, width = 10, units = "in")})

output$ui_boxplot_wl <- renderUI({ if (input$SearchButton_wl) {
  
  list(column(12,plotOutput(outputId = "boxplot_wl", height = '600px')),
       column(12,downloadBttn(outputId = 'downloadPlot_wl', label = "Download boxplot", style= "simple", color = "primary", size = "sm")))
  
}
})
