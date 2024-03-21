suppressMessages(library(shiny))
suppressMessages(library(shinyhelper))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))  # for box()
suppressMessages(library(shinycustomloader))
suppressMessages(library(shinyjs))
suppressMessages(library(fresh))
suppressMessages(library(DT))
suppressMessages(library(meta))
suppressMessages(library(ppcor))
suppressMessages(library(Hmisc))
suppressMessages(library(psych))
suppressMessages(library(MKmisc))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Biobase))
suppressMessages(library(grid))
suppressMessages(library(plotly))
suppressMessages(library(pheatmap))
suppressMessages(library(pals))
suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(cowplot))
suppressMessages(library(ggplotify))
suppressMessages(library(ggpubr))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(scCustomize))

shinyServer(function(input, output, session) {
  
  observeEvent(input$depot_ref_dc, {
    if (input$depot_ref_dc == input$depot_qry_dc) {
      updatePickerInput(session = session, inputId = 'depot_qry_dc',label = 'Adipose Depot 2',
                        choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                        selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$depot_ref_dc,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])
    }
  })
  observeEvent(input$depot_qry_dc, {
    if (input$depot_ref_dc == input$depot_qry_dc) {
      updatePickerInput(session = session, inputId = 'depot_ref_dc',label = 'Adipose Depot 1',
                        choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                        selected = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum')[!grepl(input$depot_qry_dc,list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'))])
    }
  })
  
  
  plotdata_dc <- eventReactive(input$start, {
    
    gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_dc),split='\n|\t|,'))
    
    if (length(input$cohort_dc)>1) {
      
      if (length(gene) == 1) {
        
        res_out <- c()
        for (i in input$cohort_dc) {
          
          Eset <- get(i)
          
          genes <- Eset@featureData@data[["Genes"]]
          
          Eset.pData <- pData(Eset)
          tissues <- Eset.pData$Tissue
          
          
          if (i == 'Schleinitz_et_al') {
            Eset.pData$ID <- gsub('DS_','DS ',Eset.pData$ID)
            Eset.pData$ID <- gsub('_.*','',Eset.pData$ID)
          }
          
          if (!input$gene_dc %in% genes) {
            next
          } else if (!input$depot_ref_dc %in% tissues | !input$depot_qry_dc %in% tissues) {
            next
          }
          
          
          depots_regex <- paste0(input$depot_ref_dc,'|',input$depot_qry_dc)
          depot_nums <- gsub(paste0(input$depot_ref_dc,'.*'), 1, tissues)
          depot_nums <- gsub(paste0(input$depot_qry_dc,'.*'), 2, depot_nums)
          
          gene_indices <- which(genes == input$gene_dc)
          depot_indices <- which(grepl(depots_regex, tissues))
          
          if (length(depot_indices) < 4) {print('next'); next}
          
          
          iter_ID_Ref <- gsub('_sc|_om','',Eset.pData$ID[depot_nums == '1'])
          iter_ID_Qry <- gsub('_sc|_om','',Eset.pData$ID[depot_nums == '2'])
          
          iter_Ref <- exprs(Eset)[gene_indices, (depot_nums == '1')]
          iter_Qry <- exprs(Eset)[gene_indices, (depot_nums == '2')]
          
          iter_Ref <- iter_Ref[iter_ID_Ref %in% iter_ID_Qry]
          iter_Qry <- iter_Qry[iter_ID_Qry %in% iter_ID_Ref]
          
          iter_Qry <- iter_Qry[match(iter_ID_Qry[iter_ID_Qry %in% iter_ID_Ref],iter_ID_Ref[iter_ID_Ref %in% iter_ID_Qry])]
          
          
          if (input$effect_size == 'correlation') {
            iter_res <- cor.test(x=iter_Ref, y=iter_Qry, exact=F, method=input$correlation_method)
            res_out$r <- c(res_out$r,iter_res$estimate)
            res_out$n <- c(res_out$n,length(iter_Ref))
          } else if (input$effect_size == 'SMD') {
            
            res_out$n <- c(res_out$n,length(iter_Ref))
            res_out$sd <- c(res_out$sd,sd(log2(iter_Qry)-log2(iter_Ref)))
            res_out$m <- c(res_out$m,mean(log2(iter_Qry)-log2(iter_Ref))/sd(log2(iter_Qry)-log2(iter_Ref)))
            
          }
          res_out$ID <- c(res_out$ID,i)
          
          
        }
        
        if (input$effect_size == 'correlation') {
          
          meta.temp <- meta::metacor(cor = res_out$r, n = res_out$n, studlab = res_out$ID, label.left = input$depot_ref_dc, label.right = input$depot_qry_dc)
          
          meta.temp$label.left <- input$depot_ref_dc; meta.temp$label.right <- input$depot_qry_dc
          
        } else if (input$effect_size == 'SMD') {
          meta.temp <- meta::metamean(n = res_out$n, mean = res_out$m, sd = res_out$sd, studlab = res_out$ID, null.effect = 0)
          meta.temp$sm <- 'SMD'
          
          if (all(res_out$m > 0)) {
            meta.temp$label.left <- ''; meta.temp$label.right <- paste0(input$depot_qry_dc,'-',input$depot_ref_dc)
          } else if (all(res_out$m < 0)) {
            meta.temp$label.left <- paste0(input$depot_qry_dc,'-',input$depot_ref_dc); meta.temp$label.right <- ''
          } else {
            meta.temp$label.left <- input$depot_ref_dc; meta.temp$label.right <- input$depot_qry_dc
          }
          
        }
        
        list(meta.temp,1,1,input$cohort_dc)
        
      } else if (length(gene) > 1) {
        
        Mat.r.temp <- matrix(NA, nrow = length(gene), ncol = length(input$cohort_dc))
        Mat.P.temp <- matrix(NA, nrow = length(gene), ncol = length(input$cohort_dc))
        rownames(Mat.r.temp) <- gene; rownames(Mat.P.temp) <- gene
        colnames(Mat.r.temp) <- input$cohort_dc; colnames(Mat.P.temp) <- input$cohort_dc
        
        
        for (i in input$cohort_dc) {
          
          
          Eset <- get(i)
          
          temp <<- Eset
          
          Eset.genes <- Eset@featureData@data[["Genes"]]
          iter_genes <- gene[gene %in% Eset.genes]
          
          Eset.pData <- pData(Eset)
          tissues <- Eset.pData$Tissue
          
          if (i == 'Schleinitz_et_al') {
            Eset.pData$ID <- gsub('DS_','DS ',Eset.pData$ID)
            Eset.pData$ID <- gsub('_.*','',Eset.pData$ID)
          }
          
          
          depots_regex <- paste0(input$depot_ref_dc,'|',input$depot_qry_dc)
          depot_nums <- gsub(paste0(input$depot_ref_dc,'.*'), 1, tissues)
          depot_nums <- gsub(paste0(input$depot_qry_dc,'.*'), 2, depot_nums)
          
          gene_indices <- match(iter_genes, Eset@featureData@data[['Genes']])
          depot_indices <- which(grepl(depots_regex, tissues))
          
          if (length(depot_indices) < 4) {print('next'); next}
          
          
          iter_ID_Ref <- gsub('_sc|_om','',Eset.pData$ID[depot_nums == '1'])
          iter_ID_Qry <- gsub('_sc|_om','',Eset.pData$ID[depot_nums == '2'])
          
          iter_Ref <- exprs(Eset)[gene_indices, (depot_nums == '1')]
          iter_Qry <- exprs(Eset)[gene_indices, (depot_nums == '2')]
          
          iter_Ref <- iter_Ref[ , iter_ID_Ref %in% iter_ID_Qry]
          iter_Qry <- iter_Qry[ , iter_ID_Qry %in% iter_ID_Ref]
          
          iter_Qry <- iter_Qry[ , match(iter_ID_Qry[iter_ID_Qry %in% iter_ID_Ref],iter_ID_Ref[iter_ID_Ref %in% iter_ID_Qry])]
          
          temp <- lapply(1:length(gene_indices), function(x) {cor.test(x=iter_Ref[x,], y=iter_Qry[x,], exact=F, method=input$correlation_method, alternative='two.sided')})
          
          temp_1 <- unlist(lapply(1:length(temp), function(x) {temp[[x]][['estimate']]}))
          temp_2 <- unlist(lapply(1:length(temp), function(x) {temp[[x]][['p.value']]}))
          
          
          Mat.r.temp[match(iter_genes, rownames(Mat.r.temp)),which(colnames(Mat.r.temp) == i)] <- temp_1
          Mat.P.temp[match(iter_genes, rownames(Mat.P.temp)),which(colnames(Mat.P.temp) == i)] <- temp_2
          
          
          
        }
        
        Mat.P.temp <- -log10(Mat.P.temp)
        Mat.P.temp[Mat.P.temp == Inf] <- 1e6
        
        list('heatmap',Mat.r.temp,Mat.P.temp,input$cohort_dc)
        
      }
    } else if (length(input$cohort_dc) == 1) {
      
      Eset <- get(input$cohort_dc)
      
      genes <- Eset@featureData@data[["Genes"]]
      
      Eset.pData <- pData(Eset)
      tissues <- Eset.pData$Tissue
      
      if (input$cohort_dc == 'Schleinitz_et_al') {
        Eset.pData$ID <- gsub('DS_','DS ',Eset.pData$ID)
        Eset.pData$ID <- gsub('_.*','',Eset.pData$ID)
      }
      
      if (!input$gene_dc %in% genes) {
        next
      } else if (!input$depot_ref_dc %in% tissues | !input$depot_qry_dc %in% tissues) {
        next
      }
      
      depots_regex <- paste0(input$depot_ref_dc,'|',input$depot_qry_dc)
      depot_nums <- gsub(paste0(input$depot_ref_dc,'.*'), 1, tissues)
      depot_nums <- gsub(paste0(input$depot_qry_dc,'.*'), 2, depot_nums)
      
      gene_indices <- which(genes == input$gene_dc)
      depot_indices <- which(grepl(depots_regex, tissues))
      
      if (length(depot_indices) < 4) {print('next'); next}
      
      iter_ID_Ref <- gsub('_sc|_om','',Eset.pData$ID[depot_nums == '1'])
      iter_ID_Qry <- gsub('_sc|_om','',Eset.pData$ID[depot_nums == '2'])
      
      iter_Ref <- exprs(Eset)[gene_indices, (depot_nums == '1')]
      iter_Qry <- exprs(Eset)[gene_indices, (depot_nums == '2')]
      
      iter_Ref <- iter_Ref[iter_ID_Ref %in% iter_ID_Qry]
      iter_Qry <- iter_Qry[iter_ID_Qry %in% iter_ID_Ref]
      
      iter_Qry <- iter_Qry[match(iter_ID_Qry[iter_ID_Qry %in% iter_ID_Ref],iter_ID_Ref[iter_ID_Ref %in% iter_ID_Qry])]
      
      
      log2_Exprs_Qry <- log2(iter_Qry)
      log2_Exprs_Ref <- log2(iter_Ref)
      iter_test <- t.test(log2_Exprs_Qry-log2_Exprs_Ref,alternative = 'two.sided')
      iter_df <- data.frame(group1=input$depot_ref_dc,group2=input$depot_qry_dc,p=round(iter_test[["p.value"]],4),y.position=(max(c(log2_Exprs_Ref,log2_Exprs_Qry))*1.05))
      
      df.temp <- data.frame(Exprs = c(log2_Exprs_Ref,log2_Exprs_Qry),
                            Depots = c(Eset.pData$Tissue[depot_nums == '1'][iter_ID_Ref %in% iter_ID_Qry],Eset.pData$Tissue[depot_nums == '2'][iter_ID_Qry %in% iter_ID_Ref]))
      
      df.temp$Depots <- factor(df.temp$Depots, levels=c(unique(Eset.pData$Tissue[depot_nums == '1']),unique(Eset.pData$Tissue[depot_nums == '2'])))
      
      list('boxplot',df.temp,iter_df,input$cohort_dc)
    }
  })
  
  
  dat_tab <- eventReactive(input$start, {
    temp <- plotdata_dc()
    if (class(temp[[1]])[1] == "metacor") {
      data.frame(Cohort = temp[[1]][['studlab']], Correlation = round(temp[[1]][['cor']], 2), Subjects = temp[[1]][['n']])
    } else if (class(temp[[1]])[1] == "metamean") {
      data.frame(Cohort = temp[[1]][['studlab']], mean = round(temp[[1]][['mean']],2), SD = round(temp[[1]][['sd']], 2), Subjects = temp[[1]][['n']])
    } else if (temp[[1]] == 'boxplot') {
      data.frame(Depot <- temp[[2]]$Depots, log2_FC <- round(temp[[2]]$Exprs, 2))
    } else if (temp[[1]] == 'heatmap') {
      data.frame(Cohort <- rep(colnames(temp[[2]]), nrow(temp[[2]])), Effect_Size <- c(temp[[2]]), P_val <- c(temp[[3]]))
    } 
    
  })
  
  output$data_table_dc <- DT::renderDataTable(dat_tab())
  
  
  output$plot_dc <- renderPlot({
    if (any(plotdata_dc()[[1]] == "boxplot")) {
      temp_plot <- plotdata_dc()
      
      ggpaired(data = temp_plot[[2]], x = "Depots",y = "Exprs",fill = "Depots",palette = input$discrete_col,xlab = "Depots",ylab = paste0(input$gene_dc, " expression"), title = input$cohort_dc, linetype = 'dashed') +
        #stat_compare_means(comparisons = list(c(unique(Eset.pData$Tissue[depot_nums == '1']),unique(Eset.pData$Tissue[depot_nums == '2']))), data=data.frame(iter_Qry-iter_Ref), method = 'wilcox.test') + 
        stat_pvalue_manual(temp_plot[[3]], label = 'p') + 
        theme(aspect.ratio=1, legend.position = 'none')
      
    } else if ((class(plotdata_dc()[[1]])[1] == "metacor") | (class(plotdata_dc()[[1]])[1] == "metamean")) {
      meta::forest(plotdata_dc()[[1]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA")
      
      test <<- plotdata_dc()
      
    } else if (plotdata_dc()[1] == "heatmap") {
      temp_plot <- plotdata_dc()
      
      p1 <- pheatmap::pheatmap(
        temp_plot[[2]],
        cluster_rows = T,
        cluster_cols = T,
        clustering_method = input$cluster_method,
        cellwidth = 20,
        cellheight = 20,
        border_color = NA,
        treeheight_col = 0,
        treeheight_row = 30, 
        fontsize = 15,
        main = "Correlation Heatmap",
        silent = T,
        show_colnames = T,
        color = colorRampPalette(brewer.pal(n = 7, name = input$gradient_col))(51)) # colorRampPalette(brewer.pal(n = 7, name = input$gradient_col))(51))
      
      p2 <- pheatmap::pheatmap(
        temp_plot[[3]][p1[["tree_row"]][["order"]],p1[["tree_col"]][["order"]]],
        cluster_cols = F,
        cluster_rows = F,
        display_numbers = ifelse(
          test = temp_plot[[3]][p1[["tree_row"]][["order"]],p1[["tree_col"]][["order"]]] > (-log10(0.001)),
          yes =  "***",
          no =  ifelse(test = temp_plot[[3]][p1[["tree_row"]][["order"]],p1[["tree_col"]][["order"]]] > (-log10(0.01)),
                       yes =  "**",
                       no =  ifelse(test = temp_plot[[3]][p1[["tree_row"]][["order"]],p1[["tree_col"]][["order"]]] > (-log10(0.05)),
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
        show_colnames = T,
        breaks = seq(0, 10, by = 0.1))
      
      plot_grid(as.ggplot(p1),as.ggplot(p2),ncol=2)
    }
  })

  output$ui_citation <- renderUI({if (input$start) {
    
    citation_temp <- citation[citation$Cohort_ID %in% plotdata_dc()[[4]],3:8,drop=F] %>% unique()
    
    for (i in 1:nrow(citation_temp)) {
      if (i==1) {
        citation_list <- list()
      }
      
      if (is.na(citation_temp$PMID[i])) {
        PMID_tag <- "not available"
      } else {
        PMID_tag <- as.character(a(href=citation_temp$PMID_link[i], citation_temp$PMID[i]))
      }
      
      if (is.na(citation_temp$data[i])) {
        data_tag <- "not available"
      } else {
        data_tag <- as.character(a(href=citation_temp$data_link[i], citation_temp$data[i]))
      }
      
      citation_list <- c(citation_list,paste0(citation_temp$citation_form[i],'[PMID: '), PMID_tag, "; data: ", data_tag, "]; ")
    }
    
    citation_list <- citation_list[-length(citation_list)]
    
    citation_list <- c(citation_list, "]")
    
    citation_list <- c("If you want to use this figure in your publication, please cite: Zhong et al. (portal) and ", citation_list, ".")
    
    return(HTML(paste0(paste0(citation_list),collapse = "")))
  }})


})