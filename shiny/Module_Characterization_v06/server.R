library(shiny)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(viridis)
library(DT)
library(Biobase)
library(gplots)
library(dplyr)
library(ggpubr)
library(cowplot)
library(fresh)
library(shinyjs)
library(reshape2)

# Define server logic required to draw plots and tables
shinyServer(function(input, output) {
  
  # Initialize reactive values
  v <- reactiveValues(data=NULL, data2=NULL, dataNA=NULL, data_fa=NULL)
  
  # Reactive expression to load data based on user input
  observeEvent(input$start, {
    gene_split <- strsplit(input$gene_id, split = "\n")
    if (any(rownames(adipogenesis) %in% gene_split[[1]])) {
      df <- data.frame(pData(adipogenesis),t(exprs(adipogenesis)[rownames(adipogenesis) %in% gene_split[[1]],,drop=F]))
      df <- df[!df$timepoint==20160,]
      df$timepoint <- log1p(df$timepoint)
      colnames(df)[3] <- "log_timepoint"
      df <- melt(df, id.vars = c("DNA", "ID", "log_timepoint"), variable.name = "gene", value.name = "expression")
      df$gene <- as.character(df$gene)
      v$data <- df
    }  
    if (any(rownames(SVF) %in% gene_split[[1]])) {
      # v$data2 <- data_SVF[data_SVF$gene %in% gene_split[[1]],]
      df <- data.frame(pData(SVF),t(exprs(SVF)[rownames(SVF) %in% gene_split[[1]],,drop=F]))
      df <- melt(df, id.vars = c("DNA", "ID", "timepoint"), variable.name = "gene", value.name = "expression")
      df$gene <- as.character(df$gene)
      v$data2 <- df
    }
    if (any(data_novo$gene %in% gene_split[[1]])) {
      v$dataNA <- data_novo[data_novo$gene %in% gene_split[[1]],]
    }
    if (any(rownames(FANTOM) %in% gene_split[[1]])) {

      df <- data.frame(pData(FANTOM),t(exprs(FANTOM)[rownames(FANTOM) %in% gene_split[[1]],,drop=F]))
      
      df$classfication[df$ID %in% c("adipose, donor1", "adipose, donor2", "adipose, donor3", "adipose, donor4")] <- "other tissue"
      df$classfication[df$classfication %in% c("muture adipocyte", "NiGa adipocyte")] <- "adipocyte"

      v$data_fa <- df
    }
    
    
    
  })
  
  observeEvent(input$toggle_table, {
    toggleElement(id="gene_table")
  })
  observeEvent(input$toggle_table2, {
    toggleElement(id="gene_table2")
  })
  observeEvent(input$toggle_table3, {
    toggleElement(id="gene_table3")
  })
  observeEvent(input$toggle_table4, {
    toggleElement(id="gene_table4")
  })
  
  # Define reactive expression to generate plot for NIGA timecourse
  color_vector <- reactiveVal(NULL)
  num_unique_genes <- reactiveVal(NULL)
  
  observe({
    if (!is.null(v$data)) {
      num_unique_genes(length(unique(v$data$gene)))
      if (input$discrete_col_sel == "default_2") {
        if (num_unique_genes() > length(portalcol2)) {
          color_vector(portal_col_ex(num_unique_genes()))
        } else {
          color_vector(portalcol2)
        }
      }
      if (input$discrete_col_sel == "CB2") {
        if (num_unique_genes() > 8) {
          color_vector(colorRampPalette(brewer.pal(8, input$discrete_col))(num_unique_genes()))
        } else {
          color_vector(brewer.pal(8, input$discrete_col))
        }
      }
    }
  })
  
  observe({
    if (!is.null(v$data) & input$cell_line == "hAT MSC") {

      df <- v$data
      if (input$mode == "normalized"){
        for (i in unique(df$gene)) {
           df$expression[df$gene==i] <- scale(df$expression[df$gene==i])[,1]
        }
      }

      v$plot <- ggline(data = df, x = "log_timepoint", y = "expression",numeric.x.axis = T, add = c("mean_se"), color = "gene") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])

    }
    if (!is.null(v$data) & input$cell_line == "SGBS") {
      if (input$mode == "normalized"){
        # v$plot <- ggplot(data = data_SGBS_t_z[data_SGBS_t_z$SYMBOL %in% v$data$gene,], aes(x = timepoint, y = value, group = SYMBOL, color = SYMBOL)) + 
        #   geom_point() + stat_summary(fun = "mean", geom = "line", aes(group = SYMBOL)) + theme_classic() +
        #   scale_color_manual(values = color_vector()[1:num_unique_genes()])
          
        v$plot <- ggline(data = data_SGBS_t_z[data_SGBS_t_z$SYMBOL %in% v$data$gene,], x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "SYMBOL") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])

      } else {
        # v$plot <- ggplot(data = data_SGBS_t_raw[data_SGBS_t_raw$SYMBOL %in% v$data$gene,], aes(x = timepoint, y = value, group = SYMBOL, color = SYMBOL)) + 
        #   geom_point() + stat_summary(fun = "mean", geom = "line", aes(group = SYMBOL)) + theme_classic() +
        #   scale_color_manual(values = color_vector()[1:num_unique_genes()])
        
        v$plot <- ggline(data = data_SGBS_t_raw[data_SGBS_t_raw$SYMBOL %in% v$data$gene,], x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "SYMBOL") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])
      }
    }
    
  })
  # proteome timecourse NIGA
  observe({
    if (input$cell_line == "hAT MSC"){
      if (!is.null(v$data) & input$prot_switch == TRUE & any(data_prot_z$gene %in% v$data$gene)) {
        if (input$mode == "normalized"){
          # v$plot_p <- ggplot(data = data_prot_z[data_prot_z$gene %in% v$data$gene,], aes(x = timepoint, y = value, group = gene, color = gene)) + 
          #   geom_point() + stat_summary(fun = "mean", geom = "line", aes(group = gene)) + theme_classic() +
          #   scale_color_manual(values = color_vector()[1:num_unique_genes()])
          v$plot_p <- ggline(data = data_prot_z[data_prot_z$gene %in% v$data$gene,] %>% na.omit(), x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "gene") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])
        } else {
          # v$plot_p <- ggplot(data = data_prot_raw[data_prot_raw$gene %in% v$data$gene,], aes(x = timepoint, y = value, group = gene, color = gene)) + 
          #   geom_point()+ stat_summary(fun = "mean", geom = "line", aes(group = gene)) + theme_classic() +
          #   scale_color_manual(values = color_vector()[1:num_unique_genes()])
          v$plot_p <- ggline(data = data_prot_raw[data_prot_raw$gene %in% v$data$gene,] %>% na.omit(), x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "gene") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])
        }
      } else {
        v$plot_p <- NULL  
      }
    } 
    if (input$cell_line == "SGBS"){
      if (!is.null(v$data) & input$prot_switch == TRUE) {
        if (input$mode == "normalized"){
          # v$plot_p <- ggplot(data = data_SGBS_p_z[data_SGBS_p_z$gene %in% v$data$gene,], aes(x = timepoint, y = value, group = gene, color = gene)) + 
          #   geom_point() + stat_summary(fun = "mean", geom = "line", aes(group = gene)) + theme_classic() +
          #   scale_color_manual(values = color_vector()[1:num_unique_genes()])
          v$plot_p <- ggline(data = data_SGBS_p_z[data_SGBS_p_z$gene %in% v$data$gene,] %>% na.omit(), x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "gene") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])
        } else {
          # v$plot_p <- ggplot(data = data_SGBS_p_raw[data_SGBS_p_raw$gene %in% v$data$gene,], aes(x = timepoint, y = value, group = gene, color = gene)) + 
          #   geom_point()+ stat_summary(fun = "mean", geom = "line", aes(group = gene)) + theme_classic() +
          #   scale_color_manual(values = color_vector()[1:num_unique_genes()])
          v$plot_p <- ggline(data = data_SGBS_p_raw[data_SGBS_p_raw$gene %in% v$data$gene,] %>% na.omit(), x = "timepoint", y = "value",numeric.x.axis = T, add = c("mean_se"), color = "gene") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])
        }
      } else {
        v$plot_p <- NULL  
      }
    }
  })
  
  # Define reactive expression to generate plot for SVF
  observe({
    if (!is.null(v$data2)) {
      df <- v$data2
      if (input$mode == "normalized"){
        for (i in unique(df$gene)) {
           df$expression[df$gene==i] <- scale(df$expression[df$gene==i])[,1]
        }
      }
      v$plot2 <- ggline(data = df, x = "timepoint", y = "expression",numeric.x.axis = F, add = c("mean_se"), color = "gene") +
            theme(text = element_text(family = "Red Hat Display")) +
            scale_color_manual(values = color_vector()[1:num_unique_genes()])

    }
  })
  
  # Define reactive expression to generate heatmap for NIGA timecourse
  observe({
    if (!is.null(v$data) & input$cell_line == "hAT MSC") {
      if (input$mode == "normalized") {
        if (nrow(data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]) == 1) {
          v$heatmap <- pheatmap((data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]), 
                                scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                main = "Time course hADSC differentiation",
                                labels_col = c("0 minutes","15 minutes","30 minutes","45 minutes","60 minutes",
                                               "80 minutes","100 minutes", "2 hours", "2.5 hours", "3 hours",
                                               "12 hours", "24 hours", "48 hours", "4 days", "8 days", "12 days"),silent = T)
        } else {
          v$heatmap <- pheatmap((data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]), 
                                scale = "row", cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                border_color = NA, main = "Time course hADSC differentiation",
                                labels_col = c("0 minutes","15 minutes","30 minutes","45 minutes","60 minutes",
                                               "80 minutes","100 minutes", "2 hours", "2.5 hours", "3 hours",
                                               "12 hours", "24 hours", "48 hours", "4 days", "8 days", "12 days"),silent = T)
        }
      } else {
        if (nrow(data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]) == 1) {
          v$heatmap <- pheatmap((data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]), 
                                cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                main = "Time course hADSC differentiation",
                                labels_col = c("0 minutes","15 minutes","30 minutes","45 minutes","60 minutes",
                                               "80 minutes","100 minutes", "2 hours", "2.5 hours", "3 hours",
                                               "12 hours", "24 hours", "48 hours", "4 days", "8 days", "12 days"),silent = T)
        } else {
          v$heatmap <- pheatmap((data_niga_heat[rownames(data_niga_heat) %in% v$data$gene,]), 
                                cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                border_color = NA, main = "Time course hADSC differentiation",
                                labels_col = c("0 minutes","15 minutes","30 minutes","45 minutes","60 minutes",
                                               "80 minutes","100 minutes", "2 hours", "2.5 hours", "3 hours",
                                               "12 hours", "24 hours", "48 hours", "4 days", "8 days", "12 days"),silent = T)
        }
      }
    }
    if (!is.null(v$data) & input$cell_line == "SGBS"){
      if (input$mode == "normalized") {
        if (nrow(data_SGBS_t_heat[rownames(data_SGBS_t_heat) %in% v$data$gene,]) == 1) {
          v$heatmap <- pheatmap((data_SGBS_t_heat[rownames(data_SGBS_t_heat) %in% v$data$gene,]), 
                                scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                main = "Time course SGBS differentiation",
                                labels_col = c("0 hours","24 hours","48 hours","96 hours","192 hours",
                                               "384 hours"),silent = T)
        } else {
          v$heatmap <- pheatmap((data_SGBS_t_heat[rownames(data_SGBS_t_heat) %in% v$data$gene,]), 
                                scale = "row", cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                border_color = NA, main = "Time course SGBS differentiation",
                                labels_col = c("0 hours","24 hours","48 hours","96 hours","192 hours",
                                               "384 hours"),silent = T)
        }
      } else {
        if (nrow(data_SGBS_t_heat[rownames(data_SGBS_t_heat) %in% v$data$gene,]) == 1) {
          v$heatmap <- pheatmap((data_SGBS_t_heat[rownames(data_SGBS_t_heat) %in% v$data$gene,]), 
                                cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                main = "Time course SGBS differentiation",
                                labels_col = c("0 hours","24 hours","48 hours","96 hours","192 hours",
                                               "384 hours"),silent = T)
        } else {
          v$heatmap <- pheatmap((data_SGBS_t_heat[rownames(data_SGBS_t_heat) %in% v$data$gene,]), 
                                cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                border_color = NA, main = "Time course SGBS differentiation",
                                labels_col = c("0 hours","24 hours","48 hours","96 hours","192 hours",
                                               "384 hours"),silent = T)
        }
      }
    }
  })
  
  observe({
    if (!is.null(v$data) & input$cell_line == "hAT MSC" & input$prot_switch == TRUE  & any(rownames(data_prot_heat) %in% v$data$gene) > 0) {
      if (input$mode == "normalized") {
        if (nrow(data_prot_heat[rownames(data_prot_heat) %in% v$data$gene,]) == 1) {
          v$heatmap_prot <- pheatmap((data_prot_heat[rownames(data_prot_heat) %in% v$data$gene,]), 
                                scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, fontfamily= "Red Hat Display",
                                color = farben[[input$gradient_col]], border_color = NA, 
                                main = "Time course hADSC differentiation",
                                labels_col = c("-2 days","0 days","1 day","2 days","4 days",
                                               "6 days","8 days", "10 days", "12 days"),silent = T)
        } else {
          v$heatmap_prot <- pheatmap((data_prot_heat[rownames(data_prot_heat) %in% v$data$gene,]), 
                                scale = "row", cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                border_color = NA, main = "Time course hADSC differentiation",
                                labels_col = c("-2 days","0 days","1 day","2 days","4 days",
                                               "6 days","8 days", "10 days", "12 days"),silent = T)
        }
      } else {
        if (nrow(data_prot_heat[rownames(data_prot_heat) %in% v$data$gene,]) == 1) {
          v$heatmap_prot <- pheatmap((data_prot_heat[rownames(data_prot_heat) %in% v$data$gene,]), 
                                cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                main = "Time course hADSC differentiation",
                                labels_col = c("-2 days","0 days","1 day","2 days","4 days",
                                               "6 days","8 days", "10 days", "12 days"),silent = T)
        } else {
          v$heatmap_prot <- pheatmap((data_prot_heat[rownames(data_prot_heat) %in% v$data$gene,]), 
                                cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                border_color = NA, main = "Time course hADSC differentiation",
                                labels_col = c("-2 days","0 days","1 day","2 days","4 days",
                                               "6 days","8 days", "10 days", "12 days"),silent = T)
        }
      } 
    } else {
      v$heatmap_prot <- NULL  
    }
    if (!is.null(v$data) & input$cell_line == "SGBS" & input$prot_switch == TRUE) {
      if (input$mode == "normalized") {
        if (nrow(data_SGBS_p_heat[rownames(data_SGBS_p_heat) %in% v$data$gene,]) == 1) {
          v$heatmap_prot <- pheatmap((data_SGBS_p_heat[rownames(data_SGBS_p_heat) %in% v$data$gene,]), 
                                     scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                     main = "Time course hADSC differentiation",
                                     labels_col = c("0 days","2 days","4 days","6 days",
                                                    "10 days","12 days", "14 days", "20 days"),silent = T)
        } else {
          v$heatmap_prot <- pheatmap((data_SGBS_p_heat[rownames(data_SGBS_p_heat) %in% v$data$gene,]), 
                                     scale = "row", cluster_cols = FALSE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                     border_color = NA, main = "Time course hADSC differentiation",
                                     labels_col = c("0 days","2 days","4 days","6 days",
                                                    "10 days","12 days", "14 days", "20 days"),silent = T)
        }
      } else {
        if (nrow(data_SGBS_p_heat[rownames(data_SGBS_p_heat) %in% v$data$gene,]) == 1) {
          v$heatmap_prot <- pheatmap((data_SGBS_p_heat[rownames(data_SGBS_p_heat) %in% v$data$gene,]), 
                                     cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                     main = "Time course hADSC differentiation",
                                     labels_col = c("0 days","2 days","4 days","6 days",
                                                    "10 days","12 days", "14 days", "20 days"),silent = T)
        } else {
          v$heatmap_prot <- pheatmap((data_SGBS_p_heat[rownames(data_SGBS_p_heat) %in% v$data$gene,]), 
                                     cluster_cols = FALSE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                     border_color = NA, main = "Time course hADSC differentiation",
                                     labels_col = c("0 days","2 days","4 days","6 days",
                                                    "10 days","12 days", "14 days", "20 days"),silent = T)
        }
      } 
    } else {
      v$heatmap_prot <- NULL  
    }
  })
  
  # Define reactive expression to generate heatmap for SVF
  observe({
    if (!is.null(v$data2)) {
      if (input$mode == "normalized") {
        if (nrow(data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]) == 1) {
          v$heatmap2 <- pheatmap((data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]), 
                                 cluster_cols = FALSE, cluster_rows= FALSE, scale="row", 
                                 color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                 main = "Time course SVF cells",silent = T)
        } else {
          v$heatmap2 <- pheatmap((data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]), 
                                 scale = "row", cluster_cols = FALSE, 
                                 color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                 border_color = NA, main = "Time course SVF cells",silent = T)
        }
      } else {
        if (nrow(data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]) == 1) {
          v$heatmap2 <- pheatmap((data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]), 
                                 cluster_cols = FALSE, cluster_rows= FALSE, 
                                 color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                 main = "Time course SVF cells",silent = T)
        } else {
          v$heatmap2 <- pheatmap((data_SVF_heat[rownames(data_SVF_heat) %in% v$data2$gene,]), 
                                 cluster_cols = FALSE, 
                                 color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                 border_color = NA, main = "Time course SVF cells",silent = T)
        }
      }
    }
  })
  
  # Define reactive expression to generate heatmap for novo
  observe({
    if (!is.null(v$dataNA)) {
      if (input$mode == "normalized") {
        if (nrow(data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8]) == 1) {
          v$heatmap_novo <- pheatmap((data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8]), 
                                     cluster_cols = FALSE, cluster_rows= FALSE, scale="row", 
                                     color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                     main = "FACS sorted WAT celltypes", 
                                     labels_row = data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,9],silent = T)
        } else {
          v$heatmap_novo <- pheatmap((data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8]), 
                                     scale = "row", cluster_cols = TRUE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                     border_color = NA, main = "FACS sorted WAT celltypes", 
                                     labels_row = data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,9],silent = T)
        }
      } else {
        if (nrow(data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8]) == 1) {
          v$heatmap_novo <- pheatmap((data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8]), 
                                     cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, fontfamily= "Red Hat Display",
                                     main = "FACS sorted WAT celltypes", 
                                     labels_row = data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,9],silent = T)
        } else {
          v$heatmap_novo <- pheatmap((data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,1:8]), 
                                     cluster_cols = TRUE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method, fontfamily= "Red Hat Display",
                                     border_color = NA, main = "FACS sorted WAT celltypes", 
                                     labels_row = data_novo_heat[rownames(data_novo_heat) %in% v$dataNA$transcriptID,9],silent = T)
        }
      }
    }
  })
  
  # Define reactive expression to generate heatmap for FANTOM
  observe({
    if (!is.null(v$data_fa)) {
      if (input$mode == "normalized") {
        scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
        for (i in 4:ncol(v$data_fa)) {
           v$data_fa[,i] <- scale_values(v$data_fa[,i])
        }
      }

      plot_list <- list()
      for (i in 4:ncol(v$data_fa)) {
        df <- v$data_fa[order(v$data_fa[,i],decreasing = T),]

        df$rank <- 1:nrow(df)

          p1 <- ggbarplot(df, x = "rank", y = colnames(df)[i], fill = "#adbec4", color = "#adbec4", ylab = paste0(colnames(df)[i], " expression")) +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  axis.line.x = element_line(),  
                  axis.ticks.x = element_line(),  #
                  axis.text.x = element_text())   


          group_colors <- c("adipocyte" = "#1a4659", "adipose" = "#e2c744", "other tissue" = "grey")
          annotation_bar <- data.frame(
            xmin = 0:203 + 0.5,
            xmax = 1:204 + 0.5,
            ymin = -Inf,
            ymax = -max(df[,i]) * 0.02,
            group = df$classfication
          )

          p1 <- p1 + geom_rect(data = annotation_bar, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = group), inherit.aes = FALSE) +
            scale_fill_manual(values = group_colors) +
            coord_cartesian(ylim = c(-max(df[,i]) * 0.05, NA))

          p1 <- p1 + 
            scale_x_continuous(breaks = seq(0,200,50)+0.5, labels = seq(0,200,50), name = "Rank", expand = c(0.02, 0)) +
            theme(legend.position = c(0.9,0.8),
            legend.title = element_blank())

          plot_list[[length(plot_list) + 1]] <- p1
      }

      v$heatmap_fantom <- plot_grid(plotlist = plot_list, ncol = 1, align = "hv")

    }
  })
  
  # Define output for NIGA timecourse plot
  output$linechart <- renderPlot({
    if (!is.null(v$plot)) v$plot
  })
  
  output$linechart_prot<- renderPlot({
    if (!is.null(v$plot_p)) v$plot_p
  })
  
  # Define output for SVF plot
  output$linechart2 <- renderPlot({
    if (!is.null(v$plot2)) v$plot2
  })
  
  # Define output for NIGA timecourse heatmap
  output$heatmap <- renderPlot({
    if (!is.null(v$heatmap)) v$heatmap
  })
  
  output$heatmap_prot <- renderPlot({
    if (!is.null(v$heatmap_prot)) v$heatmap_prot
  })
  
  # Define output for SVF heatmap
  output$heatmap2 <- renderPlot({
    if (!is.null(v$heatmap2)) v$heatmap2
  })
  
  # Define output for novo heatmap
  output$heatmap_novo <- renderPlot({
    if (!is.null(v$heatmap_novo)) v$heatmap_novo
  })
  
  # Define output for FANTOM heatmap
  output$heatmap_fantom <- renderPlot({
    if (!is.null(v$heatmap_fantom)) v$heatmap_fantom
  })
  
  # Define output for data tables
  output$DT_NIGA <- DT::renderDataTable({
    if (!is.null(v$data)) {
      datatable(v$data, extensions = 'Buttons', 
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )
    }
  })
  
  output$DT_SVF <- DT::renderDataTable({
    if (!is.null(v$data2)) {
      datatable(v$data2, extensions = 'Buttons', 
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )
    }
  })
  
  output$table_novo <- DT::renderDataTable({
    if (!is.null(v$dataNA)) {
      datatable(v$dataNA, extensions = 'Buttons', 
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) ) %>% 
        formatRound(columns=c('mean', 'z_norm'), digits=3)
    }
  })
  
  output$table_fantom <- DT::renderDataTable({
    if (!is.null(v$data_fa)) {
      datatable(v$data_fa[,-1], extensions = 'Buttons', 
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )
    }
  })
  
  # Download handlers for plots
  output$download_pdf1 <- downloadHandler(
    filename = function() { paste('timecourse', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$plot, device = cairo_pdf) }
  )
  
  output$download_pdf2 <- downloadHandler(
    filename = function() { paste('timecourse_heat', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap, device = cairo_pdf) }
  )
  
  output$download_pdf3 <- downloadHandler(
    filename = function() { paste('SVF_line_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$plot2, device = cairo_pdf) }
  )
  
  output$download_pdf4 <- downloadHandler(
    filename = function() { paste('SVF_heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap2, device = cairo_pdf) }
  )
  
  output$download_pdf7 <- downloadHandler(
    filename = function() { paste('WAT_fractionation_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_novo, device = cairo_pdf) }
  )
  
  output$download_pdf8 <- downloadHandler(
    filename = function() { paste('tissue_FANTOM_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_fantom, device = cairo_pdf) }
  )
  
  # Download handlers for RDS objects
  output$download_gg1 <- downloadHandler(
    filename = function() { paste('timecourse', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$plot, file) }
  )
  
  output$download_gg2 <- downloadHandler(
    filename = function() { paste('timecourse_heat', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap, file) }
  )
  
  output$download_gg3 <- downloadHandler(
    filename = function() { paste('SVF_line_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$plot2, file) }
  )
  
  output$download_gg4 <- downloadHandler(
    filename = function() { paste('SVF_heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap2, file) }
  )
  
  output$download_gg7 <- downloadHandler(
    filename = function() { paste('WAT_fractionation_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_novo, file) }
  )
  
  output$download_gg8 <- downloadHandler(
    filename = function() { paste('tissue_FANTOM_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_fantom, file) }
  )
})
