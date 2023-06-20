## forest plot

observeEvent(input$clinical_measure_fp, {
  
  updatePickerInput(session = session,
                    inputId = 'clinical_measure_adjust_fp',
                    label = 'Adjustment',
                    # selected = '',
                    choices = setdiff(clinical_measure_vector,input$clinical_measure_fp)
                    # width = "300px",
                    # options = list(size = 10,`live-search` = TRUE,`actions-box` = TRUE),
                    # multiple = TRUE
  )
})

plotdata_fp <- eventReactive(input$SearchButton_fp, {
  
  
  cor <- list()
  missing_gene_cohort <- c()
  missing_measure_cohort <- c()
  missing_measureadjust_cohort <- c()
  for (i in 1:length(input$cohort_fp)) {
    eset <- readRDS(paste0('./data/Clinical/',input$cohort_fp[i],'.rds'))
    
    
    if (!input$gene_fp %in% rownames(exprs(eset))) {
      missing_gene_cohort <- c(missing_gene_cohort,input$cohort_fp[i])
    } else {
      if (all(is.na(pData(eset)[,input$clinical_measure_fp]))) {
        missing_measure_cohort <- c(missing_measure_cohort,input$cohort_fp[i])
      } else {
        if ((length(input$clinical_measure_adjust_fp)==0)) {
          temp <- Hmisc::rcorr(exprs(eset)[input$gene_fp,],pData(eset)[,input$clinical_measure_fp],type = input$correlation_method_fp)
          cor$r <- c(cor$r,temp$r[1,2])
          cor$n <- c(cor$n,temp$n[1,2])
          cor$ID <- c(cor$ID,names(cohort_forest[cohort_forest %in% input$cohort_fp[i]]))
        } else {
          pData <- pData(eset)
          pData <- pData[,colnames(pData) %in% c(input$clinical_measure_adjust_fp,input$clinical_measure_fp)]
          pData <- pData[,colSums(is.na(pData))<nrow(pData)]  # remove columns where ALL valued are NA
          
          if (length(intersect(input$clinical_measure_adjust_fp,colnames(pData)))==0) {
            
            missing_measureadjust_cohort <- c(missing_measureadjust_cohort,paste0(names(cohort_forest[cohort_forest %in% input$cohort_fp[i]]),': ',paste0(input$clinical_measure_adjust_fp,collapse = ", ")))
            
            temp <- Hmisc::rcorr(exprs(eset)[input$gene_fp,],pData(eset)[,input$clinical_measure_fp],type = input$correlation_method_fp)
            cor$r <- c(cor$r,temp$r[1,2])
            cor$n <- c(cor$n,temp$n[1,2])
            cor$ID <- c(cor$ID,names(cohort_forest[cohort_forest %in% input$cohort_fp[i]]))
          } else {
            if (!length(setdiff(input$clinical_measure_adjust_fp,colnames(pData)))==0) {
              missing_measureadjust_cohort <- c(missing_measureadjust_cohort,paste0(names(cohort_forest[cohort_forest %in% input$cohort_fp[i]]),': ',paste0(setdiff(input$clinical_measure_adjust_fp,colnames(pData)),collapse = ", ")))
            }
            cor_df <- t(rbind(eset@assayData$exprs[input$gene_fp,],t(pData)))
            colnames(cor_df)[1] <- input$gene_fp
            cor_df_temp <- cor_df[,c(input$gene_fp,input$clinical_measure_fp,intersect(input$clinical_measure_adjust_fp,colnames(pData)))] %>% na.omit()
            cor_df_temp <- pcor.test(x = cor_df_temp[,input$gene_fp], y = cor_df_temp[,input$clinical_measure_fp], z = cor_df_temp[,intersect(input$clinical_measure_adjust_fp,colnames(pData))], method = input$correlation_method_fp)
            cor$r <- c(cor$r, cor_df_temp$estimate)
            cor$n <- c(cor$n,cor_df_temp$n)
            cor$ID <- c(cor$ID,names(cohort_forest[cohort_forest %in% input$cohort_fp[i]]))
          }
          
        }
        
        
        
        
      }
    }
    
    
  }
  
  
  
  test <- list()
  test$plot <- meta::metacor(cor = cor$r,n = cor$n, studlab = cor$ID)
  test$measure <- input$clinical_measure_fp
  test$gene <- input$gene_fp
  test$missing_gene_cohort <- names(cohort_forest[cohort_forest %in% missing_gene_cohort])
  test$missing_measure_cohort <- names(cohort_forest[cohort_forest %in% missing_measure_cohort])
  test$missing_measureadjust_cohort <- missing_measureadjust_cohort
  
  test
  
})


output$forestplot <- renderPlot({
  
  meta::forest(plotdata_fp()$plot)
  
})

output$downloadPlot_fp <- downloadHandler(filename = function() {"forestplot.pdf"},
                                          content = function(fname) {
                                            pdf(fname, width=10, height=6)
                                            meta::forest(plotdata_fp()$plot)
                                            dev.off()
                                          })



output$ui_forestplot <- renderUI({ if (input$SearchButton_fp) {
  
  list(column(12,plotOutput(outputId = "forestplot", height = '600px')),
       if (!length(plotdata_fp()$missing_gene_cohort)==0) {
         column(12,paste0('please note that ',plotdata_fp()$gene,' is not in the cohort: ',paste0(plotdata_fp()$missing_gene_cohort,collapse = ", ")))
       },
       if (!length(plotdata_fp()$missing_measure_cohort)==0) {
         column(12,paste0('please note that ',plotdata_fp()$measure,' is not in the cohort: ',paste0(plotdata_fp()$missing_measure_cohort,collapse = ", ")))
       },
       if (!length(plotdata_fp()$missing_measureadjust_cohort)==0) {
         column(12,paste0('please note that these measures are not in the cohorts: ',paste0(plotdata_fp()$missing_measureadjust_cohort,collapse = "; "),'. Be careful for the adjustment selection.'))
       },
       column(12,downloadBttn(outputId = 'downloadPlot_fp', label = "Download forest plot", style= "simple", color = "primary", size = "sm")))
  
}
})
