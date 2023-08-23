## Sex forest plot

plotdata_sex <- eventReactive(input$SearchButton_sex, {
  
  
  male_list <- list()
  female_list <- list()
  missing_gene_cohort <- c()

  for (i in 1:length(input$cohort_sex)) {
    eset <- readRDS(paste0('./data/Clinical/',input$cohort_sex[i],'.rds'))
    
    
    if (!input$gene_sex %in% rownames(exprs(eset))) {
      missing_gene_cohort <- c(missing_gene_cohort,input$cohort_sex[i])
    } else {
      
      male_list$n <- c(male_list$n, length(exprs(eset)[input$gene_sex,colnames(eset)[eset$Gender=="m"]]))
      male_list$mean <- c(male_list$mean, mean(exprs(eset)[input$gene_sex,colnames(eset)[eset$Gender=="m"]]))
      male_list$sd <- c(male_list$sd, sd(exprs(eset)[input$gene_sex,colnames(eset)[eset$Gender=="m"]]))
      male_list$ID <- c(male_list$ID, names(cohort_name[cohort_name %in% input$cohort_sex[i]]))
      
      female_list$n <- c(female_list$n, length(exprs(eset)[input$gene_sex,colnames(eset)[eset$Gender=="f"]]))
      female_list$mean <- c(female_list$mean, mean(exprs(eset)[input$gene_sex,colnames(eset)[eset$Gender=="f"]]))
      female_list$sd <- c(female_list$sd, sd(exprs(eset)[input$gene_sex,colnames(eset)[eset$Gender=="f"]]))
      female_list$ID <- c(female_list$ID, names(cohort_name[cohort_name %in% input$cohort_sex[i]]))
      

        
      
    }
    
    
  }
  
  
  
  test <- list()
  test$plot <- meta::metacont(n.e = male_list$n,
                              mean.e = male_list$mean,
                              sd.e = male_list$sd,
                              n.c = female_list$n,
                              mean.c = female_list$mean,
                              sd.c = female_list$sd,
                              studlab = male_list$ID,
                              label.e = "Male",
                              label.c = "Female",
                              label.right = "Male-specific",
                              label.left = "Female-Specific",
                              sm = "SMD",
                              method.smd = "Hedges")
  
  
  test$gene <- input$gene_sex
  test$missing_gene_cohort <- names(cohort_name[cohort_name %in% missing_gene_cohort])

  
  test
  
})


output$forestplot_sex <- renderPlot({
  
  meta::forest(plotdata_sex()$plot, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4)
  
})

output$downloadPlot_sex <- downloadHandler(filename = function() {"forestplot.pdf"},
                                           content = function(fname) {
                                             pdf(fname, width=10, height=6)
                                             meta::forest(plotdata_sex()$plot, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4)
                                             dev.off()
                                           })



output$ui_forestplot_sex <- renderUI({ if (input$SearchButton_sex) {
  
  list(column(12,plotOutput(outputId = "forestplot_sex", height = '600px')),
       if (!length(plotdata_sex()$missing_gene_cohort)==0) {
         column(12,paste0('please note that ',plotdata_sex()$gene,' is not in the cohort: ',paste0(plotdata_sex()$missing_gene_cohort,collapse = ", ")))
       },
       column(12,downloadBttn(outputId = 'downloadPlot_sex', label = "Download forest plot", style= "simple", color = "primary", size = "sm")))
  
}
})
