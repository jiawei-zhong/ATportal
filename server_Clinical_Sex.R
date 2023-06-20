## Sex boxplot

plotdata_sex <- eventReactive(input$SearchButton_sex, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_sex,'.rds'))
  pData <- pData(eset)
  pData$expression <- exprs(eset)[input$gene_sex,]
  pData$Gender[pData$Gender=="m"] <- "Male"
  pData$Gender[pData$Gender=="f"] <- "Female"
  pData$Gender <- factor(pData$Gender,levels = c("Male","Female"))
  
  
  ggboxplot(data = pData,x = "Gender",y = "expression",fill = "Gender",palette = "ucscgb",xlab = "Gender",ylab = paste0(input$gene_sex, " expression"), title = input$cohort_sex) +
    stat_compare_means(comparisons = list(c("Male","Female")), paired = FALSE, method = input$statistical_test_sex) + 
    NoLegend() + 
    theme(aspect.ratio=1)
  
  
})

output$boxplot_sex <- renderPlot({
  
  plotdata_sex()
  
})



output$downloadPlot_sex <- downloadHandler(filename = function() {"boxplot.pdf"},
                                           content = function(fname) {
                                             ggsave(filename = fname, plot = plotdata_sex(), height = 10, width = 10, units = "in")})

output$ui_boxplot_sex <- renderUI({ if (input$SearchButton_sex) {
  
  list(column(12,plotOutput(outputId = "boxplot_sex", height = '600px')),
       column(12,downloadBttn(outputId = 'downloadPlot_sex', label = "Download boxplot", style= "simple", color = "primary", size = "sm")))
  
}
})