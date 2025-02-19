suppressMessages(library(shiny))
suppressMessages(library(shinyhelper))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard)) 
suppressMessages(library(shinycustomloader))
suppressMessages(library(shinyjs))
suppressMessages(library(fresh))
suppressMessages(library(DT))
suppressMessages(library(meta))
suppressMessages(library(metafor))
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
      updatePickerInput(session = session, inputId = "depot_qry_dc",label = "Adipose Depot 2",
                        choices = list("SAT Abdomen"="SAT Abdomen","VAT Omentum"="VAT Omentum"),
                        selected = list("SAT Abdomen"="SAT Abdomen","VAT Omentum"="VAT Omentum")[!grepl(input$depot_ref_dc,list("SAT Abdomen"="SAT Abdomen","VAT Omentum"="VAT Omentum"))])
    }
  })
  observeEvent(input$depot_qry_dc, {
    if (input$depot_ref_dc == input$depot_qry_dc) {
      updatePickerInput(session = session, inputId = "depot_ref_dc",label = "Adipose Depot 1",
                        choices = list("SAT Abdomen"="SAT Abdomen","VAT Omentum"="VAT Omentum"),
                        selected = list("SAT Abdomen"="SAT Abdomen","VAT Omentum"="VAT Omentum")[!grepl(input$depot_qry_dc,list("SAT Abdomen"="SAT Abdomen","VAT Omentum"="VAT Omentum"))])
    }
  })
  observeEvent(input$navbar_id_depot, {
    if (input$navbar_id_depot == "Transcriptomics") {
      updatePickerInput(session = session, inputId = "cohort_dc", label = "Cohort",
                        choices = list("Keller, M. (2017)"="Keller_et_al",
                                       "Krieg, L. (2021)"="Krieg_et_al",
                                       "Arner, P. (2016)"="EMIFeset",
                                       "Schleinitz, D. (2020)"="Schleinitz_et_al",
                                       "Mazaki-Tovi, S. (2016)"="Mazaki_et_al",
                                       "Hoggard, N. (2012)"="Hoggard_et_al",
                                       "MacLaren, RE. (2010)"="GSE15524eset",
                                       "Salcedo-Tacuma, D. (2022)"="GSE188799eset",
                                       "Hardy, OT. (2011)"="GSE20950eset",
                                       "Du Plessis, J. (2015)"="GSE58979eset"))
    } else if (input$navbar_id_depot == "Proteomics") {
      updatePickerInput(session = session, inputId = "cohort_dc", label = "Cohort",
                        choices = list("Hruska, P. (2021)"="Hruska_2021",
                                       "Hruska, P. (2023)"="Hruska_2023",
                                       "Krieg, L. (2021)"="Krieg_et_al_proteomics"))
    }
  })
  
  val <- reactiveValues(
    gene=NULL
  )
  
  observeEvent(input$gene_dc, {val$gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_dc),split="\n|\t|,"))})

  output$ui_effect <- renderUI({
    if (length(input$cohort_dc) > 1 & length(val$gene) == 1) {
      return(list(prettyRadioButtons(inputId = "effect_size",
                                     label = "Effect Size",
                                     selected = "SMD",
                                     choices = c("Correlation" = "correlation", "Standard Mean Difference" = "SMD"),
                                     inline = TRUE))
             )
    } else if (length(input$cohort_dc) > 1 & length(val$gene) > 1) {
      return(list(prettyRadioButtons(inputId = "effect_size",
                                     label = "Effect Size",
                                     selected = "SMD",
                                     choices = c("Correlation" = "correlation", "Standard Mean Difference" = "SMD", "log2 Fold Change" = "log2_FC"),
                                     inline = TRUE))
             )
    } else {
      return(NULL)
    }
  })
  
  output$ui_corr_method <- renderUI({
    if (length(input$cohort_dc) > 1) {
      return(list(prettyRadioButtons(inputId = "correlation_method",
                                label = "Correlation method",
                                selected = "spearman",
                                choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                                inline = TRUE))
             )
    } else {
      return(NULL)
    }
  })
  
  
  plotdata_dc <- eventReactive(input$start, {
    #gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_dc),split="\n|\t|,"))
    
    gene <- val$gene
    
    validate(
      need(length(input$cohort_dc)>0 & length(gene)!=0, "Please select at least one cohort and specify at least one gene.")
    )
    
    iter_temp <- vector(mode="integer",length=length(gene))
    for (i in input$cohort_dc) {
      Eset <- get(i)
      # genes <- Eset@featureData@data[["Genes"]]
      genes <- rownames(Eset) #Jiawei edited
      iter_temp[gene %in% genes] <- iter_temp[gene %in% genes] + 1
    }
    
    validate(
      need(sum(iter_temp) > 0, "Specified gene name/names was/were not found in any of the selected cohorts.")
    )
    
    gene <- gene[iter_temp > 0]
    
    if (length(input$cohort_dc)>1) {
      
      if (length(gene) == 1) {
        
        validate(
          need(input$effect_size != "log2_FC", "Log2 fold change as an effect size measurement is not supported for forest-plots.")
        )
        
        res_out <- c()
        epsilon <- 1e-10  # 添加一个小常数
        for (i in input$cohort_dc) {
          
          Eset <- get(i)
          
          # genes <- Eset@featureData@data[["Genes"]]
          genes <- rownames(Eset) #Jiawei edited
          
          Eset.pData <- pData(Eset)
          tissues <- Eset.pData$Tissue
          
          
          # if (i == 'Schleinitz_et_al') {
          #   Eset.pData$ID <- gsub('DS_','DS ',Eset.pData$ID)
          #   Eset.pData$ID <- gsub('_.*','',Eset.pData$ID)
          # }
          
          if (!gene %in% genes) {
            next
          } else if (!input$depot_ref_dc %in% tissues | !input$depot_qry_dc %in% tissues) {
            next
          }
          
          
          depots_regex <- paste0(input$depot_ref_dc,"|",input$depot_qry_dc)
          depot_nums <- gsub(paste0(input$depot_ref_dc,".*"), 1, tissues)
          depot_nums <- gsub(paste0(input$depot_qry_dc,".*"), 2, depot_nums)
          
          gene_indices <- which(genes == gene)
          depot_indices <- which(grepl(depots_regex, tissues))
          
          if (length(depot_indices) < 4) {print("next"); next}
          
          
          iter_ID_Ref <- gsub("_sc|_om","",Eset.pData$ID[depot_nums == "1"])
          iter_ID_Qry <- gsub("_sc|_om","",Eset.pData$ID[depot_nums == "2"])
          
          iter_Ref <- exprs(Eset)[gene_indices, (depot_nums == "1")]
          iter_Qry <- exprs(Eset)[gene_indices, (depot_nums == "2")]
          
          iter_Ref <- iter_Ref[iter_ID_Ref %in% iter_ID_Qry]
          iter_Qry <- iter_Qry[iter_ID_Qry %in% iter_ID_Ref]
          
          iter_Qry <- iter_Qry[match(iter_ID_Qry[iter_ID_Qry %in% iter_ID_Ref],iter_ID_Ref[iter_ID_Ref %in% iter_ID_Qry])]
          
          
          if (input$effect_size == "correlation") {
            iter_res <- cor.test(x=iter_Ref, y=iter_Qry, exact=F, method=input$correlation_method)
            res_out$r <- c(res_out$r,iter_res$estimate)
            res_out$n <- c(res_out$n,length(iter_Ref))
          } else if (input$effect_size == "SMD") {
            
            res_out$n <- c(res_out$n, length(iter_Ref))
            res_out$sd <- c(res_out$sd, sd(iter_Qry - iter_Ref))
            res_out$m <- c(res_out$m, mean(iter_Qry - iter_Ref))
            
            res_out$Qry_m <- c(res_out$Qry_m, mean(iter_Qry))
            res_out$Qry_sd <- c(res_out$Qry_sd, sd(iter_Qry))
            
            res_out$Ref_m <- c(res_out$Ref_m, mean(iter_Ref))
            res_out$Ref_sd <- c(res_out$Ref_sd, sd(iter_Ref))
            
          }
          res_out$ID <- c(res_out$ID, names(cohort_name[cohort_name==i]))
          res_out$Qry_n <- res_out$n
          res_out$Ref_n <- res_out$n
          
          
          res_out$sd[res_out$sd == 0] <- epsilon
          res_out$Qry_sd[res_out$Qry_sd == 0] <- epsilon
          res_out$Ref_sd[res_out$Ref_sd == 0] <- epsilon

          
        }
        
        if (input$effect_size == "correlation") {
          
          meta.temp <- meta::metacor(cor = res_out$r, n = res_out$n, studlab = res_out$ID, label.left = input$depot_ref_dc, label.right = input$depot_qry_dc)
          
          meta.temp$label.left <- input$depot_ref_dc; meta.temp$label.right <- input$depot_qry_dc
          
        } else if (input$effect_size == "SMD") {
          
          
          res_out_df <- as.data.frame(res_out)
          
          
          temp <- tryCatch({
            result <- metafor::escalc(measure = "SMN",
                             mi = m, 
                             ni = n,
                             sdi = sd,
                             data = res_out_df)
            result
          }, error = function(e) {
            tryCatch({
              result <- metafor::escalc(measure = "SMN",
                               mi = m, 
                               ni = n,
                               sdi = sd,
                               data = res_out_df,
                               control = list(maxiter = 1000))
              result
            }, error = function(e) {
              tryCatch({
                result <- metafor::escalc(measure = "SMN",
                                 mi = m, 
                                 ni = n,
                                 sdi = sd,
                                 data = res_out_df,
                                 control = list(stepadj = 0.5))
                result
              }, error = function(e) {
                NULL
                # Handle the error if needed
                # You can also leave this part empty to just skip the current iteration
              })
            })
          })
          
          
          meta.temp <- tryCatch({
            result <- meta::metagen(TE = temp$yi,seTE = temp$vi^.5, method.tau = 'REML',
                                    method.tau.ci = 'BJ', comb.random = TRUE, comb.fixed = TRUE,
                                    label.c = input$depot_ref_dc,label.e = input$depot_qry_dc,label.left = input$depot_ref_dc,label.right = input$depot_qry_dc,
                                    sm = 'SMD')
          }, error = function(e) {
            tryCatch({
              result <- meta::metagen(TE = temp$yi,seTE = temp$vi^.5, method.tau = 'REML',
                                      method.tau.ci = 'BJ', comb.random = TRUE, comb.fixed = TRUE,
                                      label.c = input$depot_ref_dc,label.e = input$depot_qry_dc,label.left = input$depot_ref_dc,label.right = input$depot_qry_dc,
                                      sm = 'SMD',control = list(maxiter = 1000))
            }, error = function(e) {
              tryCatch({
                result <- meta::metagen(TE = temp$yi,seTE = temp$vi^.5, method.tau = 'REML',
                                        method.tau.ci = 'BJ', comb.random = TRUE, comb.fixed = TRUE,
                                        label.c = input$depot_ref_dc,label.e = input$depot_qry_dc,label.left = input$depot_ref_dc,label.right = input$depot_qry_dc,
                                        sm = 'SMD',control = list(stepadj = 0.5))
              }, error = function(e) {
                NULL
              })
            })
          })
          
          if (!is.null(meta.temp)) {
            class(meta.temp) <- c("metacont","meta")
            
            
            meta.temp$mean.e <- res_out_df$Qry_m
            meta.temp$mean.c <- res_out_df$Ref_m
            
            
            meta.temp$sd.e <- res_out_df$Qry_sd
            meta.temp$sd.c <- res_out_df$Ref_sd
            
            
            meta.temp$n.e <- res_out_df$Qry_n
            meta.temp$n.c <- res_out_df$Ref_n
            
            meta.temp$studlab <- res_out_df$ID
            
          }
          
          
          # meta.temp <- meta::metamean(n = res_out$n, mean = res_out$m, sd = res_out$sd, studlab = res_out$ID, null.effect = 0)
          # meta.temp$sm <- "SMD"
          
          # if (all(res_out$m > 0)) {
          #   meta.temp$label.left <- ""; meta.temp$label.right <- paste0(input$depot_qry_dc,"-",input$depot_ref_dc)
          # } else if (all(res_out$m < 0)) {
          #   meta.temp$label.left <- paste0(input$depot_qry_dc,"-",input$depot_ref_dc); meta.temp$label.right <- ""
          # } else {
          #   meta.temp$label.left <- input$depot_ref_dc; meta.temp$label.right <- input$depot_qry_dc
          # }
        }
        
        if (input$effect_size == "correlation") {
          list(plottype="forest", plot=meta.temp, cohort=res_out$ID)
        } else {
          list(plottype="forest", plot=meta.temp, cohort=res_out$ID, df=res_out_df)
        }
        
      } else if (length(gene) > 1) {
        
        Mat.r.temp <- matrix( nrow = length(gene), ncol = length(input$cohort_dc))
        Mat.P.temp <- matrix( nrow = length(gene), ncol = length(input$cohort_dc))
        rownames(Mat.r.temp) <- gene; rownames(Mat.P.temp) <- gene
        colnames(Mat.r.temp) <- input$cohort_dc; colnames(Mat.P.temp) <- input$cohort_dc
        
        for (i in input$cohort_dc) {
          
          
          Eset <- get(i)
          
          temp <<- Eset
          
          # Eset.genes <- Eset@featureData@data[["Genes"]]
          Eset.genes <- rownames(Eset) #Jiawei edited
          
          iter_genes <- gene[gene %in% Eset.genes]
          
          Eset.pData <- pData(Eset)
          tissues <- Eset.pData$Tissue
          
          # if (i == 'Schleinitz_et_al') {
          #   Eset.pData$ID <- gsub('DS_','DS ',Eset.pData$ID)
          #   Eset.pData$ID <- gsub('_.*','',Eset.pData$ID)
          # }
          
          
          depots_regex <- paste0(input$depot_ref_dc,"|",input$depot_qry_dc)
          depot_nums <- gsub(paste0(input$depot_ref_dc,".*"), 1, tissues)
          depot_nums <- gsub(paste0(input$depot_qry_dc,".*"), 2, depot_nums)
          
          # gene_indices <- match(iter_genes, Eset@featureData@data[["Genes"]])
          gene_indices <- match(iter_genes, rownames(Eset)) # Jiawei edited 
          depot_indices <- which(grepl(depots_regex, tissues))
          
          if (length(depot_indices) < 4) {print("next"); next}
          
          
          iter_ID_Ref <- gsub("_sc|_om","",Eset.pData$ID[depot_nums == "1"])
          iter_ID_Qry <- gsub("_sc|_om","",Eset.pData$ID[depot_nums == "2"])
          
          iter_Ref <- exprs(Eset)[gene_indices, (depot_nums == "1"),drop=F]
          iter_Qry <- exprs(Eset)[gene_indices, (depot_nums == "2"),drop=F]
          
          iter_Ref <- iter_Ref[ , iter_ID_Ref %in% iter_ID_Qry,drop=F]
          iter_Qry <- iter_Qry[ , iter_ID_Qry %in% iter_ID_Ref,drop=F]
          
          iter_Qry <- iter_Qry[ , match(iter_ID_Qry[iter_ID_Qry %in% iter_ID_Ref],iter_ID_Ref[iter_ID_Ref %in% iter_ID_Qry]),drop=F]
          
          
          if (input$effect_size == "correlation") {
            
            temp <- lapply(1:length(gene_indices), function(x) {cor.test(x=iter_Ref[x,], y=iter_Qry[x,], exact=F, method=input$correlation_method, alternative="two.sided")})
            
            temp_1 <- unlist(lapply(1:length(temp), function(x) {temp[[x]][["estimate"]]}))
            temp_2 <- unlist(lapply(1:length(temp), function(x) {temp[[x]][["p.value"]]}))
            
          } else if (input$effect_size == "SMD") {
            
            temp_1 <- unlist(lapply(1:length(gene_indices), function(x) {mean(iter_Qry[x,]-iter_Ref[x,])/sd(iter_Qry[x,]-iter_Ref[x,])}))
            temp_2 <- unlist(lapply(1:length(gene_indices), function(x) {t.test(iter_Qry[x,]-iter_Ref[x,],mu=0)[["p.value"]]}))
            
          } else if (input$effect_size == "log2_FC") {
            
            temp_1 <- unlist(lapply(1:length(gene_indices), function(x) {mean(iter_Qry[x,]-iter_Ref[x,])}))
            temp_2 <- unlist(lapply(1:length(gene_indices), function(x) {t.test(iter_Qry[x,]-iter_Ref[x,],mu=0)[["p.value"]]}))
            
          }
          
          Mat.r.temp[match(iter_genes, rownames(Mat.r.temp)),which(colnames(Mat.r.temp) == i)] <- temp_1
          Mat.P.temp[match(iter_genes, rownames(Mat.P.temp)),which(colnames(Mat.P.temp) == i)] <- temp_2
          
        }
        
        colnames(Mat.r.temp) <- names(cohort_name[match(colnames(Mat.r.temp),cohort_name)])
        colnames(Mat.P.temp) <- names(cohort_name[match(colnames(Mat.P.temp),cohort_name)])
        
        Mat.P.temp <- Mat.P.temp[rowSums(is.na(Mat.r.temp)) != ncol(Mat.r.temp), ]
        Mat.r.temp <- Mat.r.temp[rowSums(is.na(Mat.r.temp)) != ncol(Mat.r.temp), ]
        
        
        Mat.P.temp <- -log10(Mat.P.temp)
        Mat.P.temp[Mat.P.temp == Inf] <- 1e6
        
        list(plottype="heatmap",round(Mat.r.temp,5),round(Mat.P.temp,5), colnames(Mat.r.temp))
        
      }
    } else if (length(input$cohort_dc) == 1) {
      
      validate(
        need(length(gene) == 1, "Please enter a single gene to generate a boxplot.")
      )
      
      validate(
        need(input$effect_size == "log2_FC", "Correlation as an effect size measurement is not supported for forest-plots.")
      )
      
      Eset <- get(input$cohort_dc)
      
      # genes <- Eset@featureData@data[["Genes"]]
      genes <- rownames(Eset) #Jiawei edited
      
      Eset.pData <- pData(Eset)
      tissues <- Eset.pData$Tissue
      
      # if (input$cohort_dc == 'Schleinitz_et_al') {
      #   Eset.pData$ID <- gsub('DS_','DS ',Eset.pData$ID)
      #   Eset.pData$ID <- gsub('_.*','',Eset.pData$ID)
      # }
      
      
      depots_regex <- paste0(input$depot_ref_dc,"|",input$depot_qry_dc)
      depot_nums <- gsub(paste0(input$depot_ref_dc,".*"), 1, tissues)
      depot_nums <- gsub(paste0(input$depot_qry_dc,".*"), 2, depot_nums)
      
      gene_indices <- which(genes == input$gene_dc)
      depot_indices <- which(grepl(depots_regex, tissues))
      
      if (length(depot_indices) < 4) {print("next"); next}
      
      iter_ID_Ref <- gsub("_sc|_om","",Eset.pData$ID[depot_nums == "1"])
      iter_ID_Qry <- gsub("_sc|_om","",Eset.pData$ID[depot_nums == "2"])
      
      iter_Ref <- exprs(Eset)[gene_indices, (depot_nums == "1")]
      iter_Qry <- exprs(Eset)[gene_indices, (depot_nums == "2")]
      
      iter_Ref <- iter_Ref[iter_ID_Ref %in% iter_ID_Qry]
      iter_Qry <- iter_Qry[iter_ID_Qry %in% iter_ID_Ref]
      
      iter_Qry <- iter_Qry[match(iter_ID_Qry[iter_ID_Qry %in% iter_ID_Ref],iter_ID_Ref[iter_ID_Ref %in% iter_ID_Qry])]
      
      
      log2_Exprs_Qry <- log2(iter_Qry)
      log2_Exprs_Ref <- log2(iter_Ref)
      iter_test <- t.test(log2_Exprs_Qry-log2_Exprs_Ref,alternative = "two.sided")
      iter_df <- data.frame(group1=input$depot_ref_dc,group2=input$depot_qry_dc,p=round(iter_test[["p.value"]],4),y.position=(max(c(log2_Exprs_Ref,log2_Exprs_Qry))*1.05))
      
      df.temp <- data.frame(Exprs = c(log2_Exprs_Ref,log2_Exprs_Qry),
                            Depots = c(Eset.pData$Tissue[depot_nums == "1"][iter_ID_Ref %in% iter_ID_Qry],Eset.pData$Tissue[depot_nums == "2"][iter_ID_Qry %in% iter_ID_Ref]))
      
      df.temp$Depots <- factor(df.temp$Depots, levels=c(unique(Eset.pData$Tissue[depot_nums == "1"]),unique(Eset.pData$Tissue[depot_nums == "2"])))
      
      list(plottype="boxplot",df.temp,iter_df, names(cohort_name[match(input$cohort_dc,cohort_name)]))
      
    }
  })
  
  
  dat_tab <- eventReactive(input$start, {
    temp <- plotdata_dc()
    if (class(temp[[2]])[1] == "metacor") {
      data.frame(Cohort = temp[[2]][["studlab"]], Correlation = round(temp[[2]][["cor"]], 2), Subjects = temp[[2]][["n"]])
    } else if (class(temp[[2]])[1] == "metacont") {
      # data.frame(Cohort = temp[[2]][["studlab"]], SMD = round(temp[[2]][["mean"]],2), SD = round(temp[[2]][["sd"]], 2), Subjects = temp[[2]][["n"]])
      temp[[4]]
    } else if (temp[[1]] == "boxplot") {
      aggregate(.~Depot,data=data.frame(Depot = temp[[2]]$Depots, log2_mean_expression = temp[[2]]$Exprs), FUN = function(x) round(mean(x),2))
    } else if (temp[[1]] == "heatmap") {
      temp1 <- temp[[2]]; colnames(temp1) <- paste("Effect Size:",colnames(temp1))
      temp2 <- temp[[3]]; colnames(temp2) <- paste("P-value:",colnames(temp2))
      cbind(as.data.frame(temp1),as.data.frame(temp2))
      
    } 
    
  })
  
  output$data_table_dc_1 <- DT::renderDataTable(
    {dat_tab()},
    extensions = c("Buttons"),
    options = list(scrollX = TRUE,
                   pageLength = 10,
                   lengthMenu = c(10, 25, 50, 100),
                   dom = "Blfrtip",
                   buttons = c("copy", "csv", "excel", "pdf", "print"))
  )
  output$data_table_dc_2 <- DT::renderDataTable(
    {dat_tab()},
    extensions = c("Buttons"),
    options = list(scrollX = TRUE,
                   pageLength = 10,
                   lengthMenu = c(10, 25, 50, 100),
                   dom = "Blfrtip",
                   buttons = c("copy", "csv", "excel", "pdf", "print"))
  )
  
  
  
  output$plot_dc_1 <- renderPlot({
    
    if (is.null(input$cohort_dc)) {
      NULL
    } else {
      
      if (is.null(plotdata_dc()$plot)) {
        # 如果 plot 对象不存在，则可以选择在控制台打印信息，或者进行其他的错误处理操作
        print("No plot data available. Check the data input and processing steps.")
        return(NULL)  # 返回 NULL，避免尝试绘制不存在的图形
      }
      
      
      if (plotdata_dc()[[1]][1] == "boxplot") {
        temp_plot <- plotdata_dc()
        
        ggpaired(data = temp_plot[[2]], x = "Depots",y = "Exprs",fill = "Depots",palette = input$discrete_col,xlab = "Depots",ylab = paste0(input$gene_dc, " expression"), title = names(cohort_name[grep(input$cohort_dc,cohort_name)]), linetype = "dashed") +
          stat_pvalue_manual(temp_plot[[3]], label = "p") +
          theme(aspect.ratio=1, legend.position = "none")
        
      } else if (plotdata_dc()[[1]][1] == "forest") {
        discrete_color <- portalcol2[1:2]
        # meta::forest(plotdata_dc(), test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA")
        # meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        
        if (is.na(plotdata_dc()[[2]]$pval.Q)) {
          meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        } else {
          if (plotdata_dc()[[2]]$pval.Q<0.05) {
            if (plotdata_dc()[[2]]$TE.random>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          } else {
            if (plotdata_dc()[[2]]$TE.common>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          }
        }
        
      }  else if (plotdata_dc()[[1]][1] == "heatmap") {
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
          main = "Heatmap of Effect Size",
          silent = T,
          show_colnames = T,
          fontfamily= "Red Hat Display",
          color = get(input$gradient_col)) 
        
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
          fontfamily= "Red Hat Display",
          breaks = seq(0, 10, by = 0.1))
        
        plot_grid(as.ggplot(p1),as.ggplot(p2),ncol=2)
      }
    }
    
  })
  
  output$plot_dc_2 <- renderPlot({
    
    if (is.null(input$cohort_dc)) {
      NULL
    } else {
      if (plotdata_dc()[[1]][1] == "boxplot") {
        temp_plot <- plotdata_dc()
        
        ggpaired(data = temp_plot[[2]], x = "Depots",y = "Exprs",fill = "Depots",palette = input$discrete_col,xlab = "Depots",ylab = paste0(input$gene_dc, " expression"), title = input$cohort_dc, linetype = "dashed") +
          stat_pvalue_manual(temp_plot[[3]], label = "p") +
          theme(aspect.ratio=1, legend.position = "none")
        
      } else if (plotdata_dc()[1] == "forest") {
        discrete_color <- portalcol2[1:2]
        # meta::forest(plotdata_dc(), test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA")
        # meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        
        if (is.na(plotdata_dc()[[2]]$pval.Q)) {
          meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        } else {
          if (plotdata_dc()[[2]]$pval.Q<0.05) {
            if (plotdata_dc()[[2]]$TE.random>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          } else {
            if (plotdata_dc()[[2]]$TE.common>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          }
        }
        
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
          main = "Heatmap of Effect Size",
          silent = T,
          show_colnames = T,
          fontfamily= "Red Hat Display",
          color = get(input$gradient_col)) 
        
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
          fontfamily= "Red Hat Display",
          breaks = seq(0, 10, by = 0.1))
        
        plot_grid(as.ggplot(p1),as.ggplot(p2),ncol=2)
      }
    }
  })
  
  output$dc_pdf_dl1 <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      
      if (plotdata_dc()[1] == "boxplot") {
        
        temp_plot <- plotdata_dc()
        
        ggsave(file, plot =
                 ggpaired(data = temp_plot[[2]], x = "Depots",y = "Exprs",fill = "Depots",palette = input$discrete_col,xlab = "Depots",ylab = paste0(input$gene_dc, " expression"), title = input$cohort_dc, linetype = "dashed") +
                 stat_pvalue_manual(temp_plot[[3]], label = "p") +
                 theme(aspect.ratio=1, legend.position = "none"),
               height = 20, width = 20, units = "in", device = cairo_pdf)
        
      } else if (plotdata_dc()[1] == "forest") {
        
        cairo_pdf(file, width=13, height=6)
        
        discrete_color <- portalcol2[1:2]
        # meta::forest(plotdata_dc(), test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA")
        # meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        
        if (is.na(plotdata_dc()[[2]]$pval.Q)) {
          meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        } else {
          if (plotdata_dc()[[2]]$pval.Q<0.05) {
            if (plotdata_dc()[[2]]$TE.random>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          } else {
            if (plotdata_dc()[[2]]$TE.common>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          }
        }
        
        
        
        dev.off()
        
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
          main = "Heatmap of Effect Size",
          silent = T,
          show_colnames = T,
          fontfamily= "Red Hat Display",
          color = get(input$gradient_col)) 
        
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
          fontfamily= "Red Hat Display",
          breaks = seq(0, 10, by = 0.1))
        
        ggsave(file, plot = plot_grid(as.ggplot(p1),as.ggplot(p2),ncol=2), height = 20, width = 20, units = "in", device = cairo_pdf)
      }
    }
  )
  
  output$dc_pdf_dl2 <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), ".pdf", sep="")
    },
    content = function(file) {
      
      if (plotdata_dc()[1] == "boxplot") {
        
        temp_plot <- plotdata_dc()
        
        ggsave(file, plot =
                 ggpaired(data = temp_plot[[2]], x = "Depots",y = "Exprs",fill = "Depots",palette = input$discrete_col,xlab = "Depots",ylab = paste0(input$gene_dc, " expression"), title = input$cohort_dc, linetype = "dashed") +
                 stat_pvalue_manual(temp_plot[[3]], label = "p") +
                 theme(aspect.ratio=1, legend.position = "none"),
               height = 20, width = 20, units = "in", device = cairo_pdf)
        
      } else if (plotdata_dc()[1] == "forest") {
        
        cairo_pdf(file, width=13, height=6)
        
        discrete_color <- portalcol2[1:2]
        # meta::forest(plotdata_dc(), test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, layout = "JAMA")
        # meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        
        if (is.na(plotdata_dc()[[2]]$pval.Q)) {
          meta::forest(plotdata_dc()[[2]], test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        } else {
          if (plotdata_dc()[[2]]$pval.Q<0.05) {
            if (plotdata_dc()[[2]]$TE.random>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          } else {
            if (plotdata_dc()[[2]]$TE.common>0) {
              meta::forest(plotdata_dc()[[2]], sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              meta::forest(plotdata_dc()[[2]], sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            }
          }
        }
        
        
        
        dev.off()
        
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
          main = "Heatmap of Effect Size",
          silent = T,
          show_colnames = T,
          fontfamily= "Red Hat Display",
          color = get(input$gradient_col)) 
        
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
          fontfamily= "Red Hat Display",
          breaks = seq(0, 10, by = 0.1))
        
        ggsave(file, plot = plot_grid(as.ggplot(p1),as.ggplot(p2),ncol=2), height = 20, width = 20, units = "in", device = cairo_pdf)
      }
    }
  )
  
  output$ui_citation <- renderUI({if (input$start) {
    
    cohort_ID <- plotdata_dc()
    cohort_ID <- cohort_ID[[3]]
    cohort_ID <- cohort_name[match(cohort_ID, names(cohort_name))]

    if (length(cohort_ID) > 0) {
      citation_temp <- citation[citation$Cohort_ID %in% cohort_ID,3:9,drop=F] %>% unique()
      
      for (i in 1:nrow(citation_temp)) {
        if (i==1) {
          citation_list <- list()
        }
        
        if (is.na(citation_temp$PMID[i])) {
          PMID_tag <- "not available"
        } else {
          PMID_tag <- as.character(a(href=citation_temp$PMID_link[i], target = "_blank", citation_temp$PMID[i]))
        }
        
        if (is.na(citation_temp$data[i])) {
          data_tag <- "not available"
        } else {
          data_tag <- as.character(a(href=citation_temp$data_link[i], target = "_blank", citation_temp$data[i]))
        }
        
        citation_list <- c(citation_list,paste0(citation_temp$citation_form[i],'[PMID: '), PMID_tag, "; data: ", data_tag, "]; ")
      }
      
      citation_list <- citation_list[-length(citation_list)]
      
      citation_list <- c(citation_list, "]")
      
      citation_list <- c("If you want to use this figure in your publication, please cite: Zhong et al. (portal) and ", citation_list, ".")
      
      return(HTML(paste0(paste0(citation_list),collapse = "")))       
    } else {
      return(NULL)
    }
  }})
  
})
