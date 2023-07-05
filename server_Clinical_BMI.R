## BMI boxplot


observeEvent(input$cohort_box_BMI, {
  
  
  BMI_catelogy <- readRDS("./data/Clinical/cohort_BMI_list.rds")
  BMI_catelogy <- BMI_catelogy[[input$cohort_box_BMI]]
  name <- BMI_catelogy
  name[name=="Underweight"] <- "Underweight (BMI < 18.5)"
  name[name=="Healthy Weight"] <- "Healthy Weight (BMI: 18.5 ~ 25)"
  name[name=="Overweight"] <- "Overweight (BMI: 25 ~ 30)"
  name[name=="Obese I"] <- "Obese I (BMI: 30 ~ 35)"
  name[name=="Obese II"] <- "Obese II (BMI: 35 ~ 40)"
  name[name=="Obese III"] <- "Obese III (BMI > 40)"
  name[name=="Obese"] <- "Obese (BMI > 30)"
  names(BMI_catelogy) <- name
  
  updatePickerInput(session = session,
                    inputId = 'group_box_BMI',
                    label = 'Choose BMI group',
                    selected = NA,
                    choices = BMI_catelogy
                    # width = "300px",
                    # options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                    # multiple = TRUE
  )
  
  
  
})



plotdata_box_BMI <- eventReactive(input$SearchButton_box_BMI, {
  
  eset <- readRDS(paste0('./data/Clinical/',input$cohort_box_BMI,'.rds'))
  pData <- pData(eset)
  pData <- pData[,colSums(is.na(pData))<nrow(pData)]
  pData$expression <- exprs(eset)[input$gene_BMI,]
  
  BMI_catelogy <- readRDS("./data/Clinical/cohort_BMI_list.rds")
  BMI_catelogy <- BMI_catelogy[[input$cohort_box_BMI]]
  
  if ("Obese" %in% BMI_catelogy) {
    pData_new <- pData[pData$BMI_catelogy %in% c("Obese I","Obese II","Obese III"),]
    pData_new$BMI_catelogy <- "Obese"
    pData <- rbind(pData,pData_new)
  }
  
  BMI_catelogy <- input$group_box_BMI
  pData <- pData[pData$BMI_catelogy %in% input$group_box_BMI,]
  
  
  compare_list <- list()
  
  for (i in 1:(length(BMI_catelogy)-1)) {
    for (j in (i+1):length(BMI_catelogy)) {
      
      if ((!BMI_catelogy[j]=="Obese") | (!BMI_catelogy[i] %in% c("Obese I","Obese II","Obese III"))) {
        compare_list <- append(compare_list,list(BMI_catelogy[c(i,j)]))
      }
      
      
    }
  }
  
  
  
  
  
  p <- ggboxplot(data = pData,x = "BMI_catelogy",y = "expression",fill = "BMI_catelogy",palette = "ucscgb",xlab = "BMI catelogies",ylab = paste0(input$gene_BMI, " expression"), title = input$cohort_box_BMI)+
    NoLegend() +
    theme(aspect.ratio=1,
          axis.text.x = element_text(face = "plain",size = 11,angle = 45, hjust = 1, vjust = 1)) +
    stat_compare_means(comparisons = compare_list, paired = F, method = input$statistical_test_cor_BMI)
  
  p
  
})

output$boxplot_BMI <- renderPlot({
  
  plotdata_box_BMI()
  
})



output$downloadPlot_box_BMI <- downloadHandler(filename = function() {"boxplot.pdf"},
                                               content = function(fname) {
                                                 ggsave(filename = fname, plot = plotdata_box_BMI(), height = 10, width = 10, units = "in")})

output$ui_boxplot_BMI <- renderUI({ if (input$SearchButton_box_BMI) {
  
  list(column(12,plotOutput(outputId = "boxplot_BMI", height = '600px')),
       column(12,downloadBttn(outputId = 'downloadPlot_box_BMI', label = "Download boxplot", style= "simple", color = "primary", size = "sm")))
  
}
})