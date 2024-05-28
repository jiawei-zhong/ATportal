library(shiny)
library(pheatmap)
library(ggplotify)
library(ggplot2)
library(ggrepel)
library(viridis)
library(DT)
library(Biobase)
library(gplots)
library(dplyr)
library(cowplot)
library(fresh)
library(shinyjs)

# Define server logic required to draw plots and tables
shinyServer(function(input, output) {
  
  # Initialize reactive values
  v <- reactiveValues(hyp=NULL, tnf=NULL, ins=NULL, gls=NULL, cd248=NULL, aqp7=NULL, c14=NULL, stim=NULL, inflamed=NULL)
  
  # Reactive expression to load data based on user input
  observeEvent(input$start, {
    gene_split <- strsplit(input$gene_id, split = "\n")  
    v$hyp <- data_hyp[data_hyp$GENE %in% gene_split[[1]],]
    v$tnf <- data_tnf[data_tnf$gene %in% gene_split[[1]],]
    v$ins <- data_ins[data_ins$genes %in% gene_split[[1]],]
    v$gls <- data_gls[data_gls$gene %in% gene_split[[1]],]
    v$cd248 <- data_cd248[data_cd248$GENE %in% gene_split[[1]],]
    v$aqp7 <- data_aqp7[data_aqp7$GENE %in% gene_split[[1]],]
    v$c14 <- data_c14[data_c14$GENE %in% gene_split[[1]],]
    v$stim <- data_stim_heat[data_stim_heat$SYMBOL %in% gene_split[[1]],]
    v$inflamed <-  data_infl_heat[data_infl_heat$SYMBOL %in% gene_split[[1]],]
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
  inflamed_select <- reactiveVal(NULL)
  inflamed_select2 <- reactiveVal(NULL)
  
  observe({
    if (!is.null(v$inflamed)){
      if (input$inflamedAT_compare == "24h Treatment vs Control") {
        inflamed_select(c(1,2,5,6,9,10))
        inflamed_select2(data_infl$Compare24h)
      }
      if (input$inflamedAT_compare == "48h Treatment vs Control") {
        inflamed_select(c(3,4,7,8,11,12))
        inflamed_select2(data_infl$Compare48h)
      }
      if (input$inflamedAT_compare == "Treatment 24h vs 48h") {
        inflamed_select(c(2,4,6,8,10,12))
        inflamed_select2(data_infl$CompareActive)
      }
      if (input$inflamedAT_compare == "All Treatment vs Control") {
        inflamed_select(1:12)
        inflamed_select2(data_infl$CompareAll)
      }
    }
  })
  
  observe({
    if (!is.null(v$hyp)) {
      num_unique_genes(length(unique(v$hyp$GENE)))
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
#observe events for correct plot generation  
#Hypoxia 
  #heatmap
  observe({
    if (!is.null(v$hyp)) {
      if (input$mode == "normalized") {
        if (nrow(data_hyp_heat[rownames(data_hyp_heat) %in% rownames(v$hyp),]) == 1) {
          v$heatmap_hyp <- pheatmap((data_hyp_heat[rownames(data_hyp_heat) %in% rownames(v$hyp),]), 
                                scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, 
                                main = "Hypoxia", labels_row = v$hyp$GENE)
        } else {
          v$heatmap_hyp <- pheatmap((data_hyp_heat[rownames(data_hyp_heat) %in% rownames(v$hyp),]), 
                                scale = "row", cluster_cols = TRUE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                border_color = NA, main = "Hypoxia", labels_row = v$hyp$GENE)
        }
      } else {
        if (nrow(data_hyp_heat[rownames(data_hyp_heat) %in% rownames(v$hyp),]) == 1) {
          v$heatmap_hyp <- pheatmap((data_hyp_heat[rownames(data_hyp_heat) %in% rownames(v$hyp),]), 
                                cluster_cols = FALSE, cluster_rows= FALSE, 
                                color = farben[[input$gradient_col]], border_color = NA, 
                                main = "Hypoxia", labels_row = v$hyp$GENE)
        } else {
          v$heatmap_hyp <- pheatmap((data_hyp_heat[rownames(data_hyp_heat) %in% rownames(v$hyp),]), 
                                cluster_cols = FALSE, 
                                color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                border_color = NA, main = "Hypoxia", labels_row = v$hyp$GENE)
        }
      }
    }
  })

  #Volcano
 observe({
   if(!is.null(v$hyp)){
     v$volcano_hyp <- ggplot(data = data_hyp, aes(x=logFC, y=-log10(adj.P.Val), label=GENE)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$hyp, aes(x=logFC, y=-log10(adj.P.Val)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$hyp, aes(label=GENE)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })

 #TNF SGBS AND NIGA
 #heatmap
 observe({
   if(input$cell_line == "hMADS"){
     if (!is.null(v$tnf)) {
       if (input$mode == "normalized") {
         if (nrow(data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,1:6]) == 1) {
           v$heatmap_tnf <- pheatmap((data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,1:6]), 
                                     scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, 
                                     main = "TNF stimulation", labels_row = data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,7])
         } else {
           v$heatmap_tnf <- pheatmap((data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,1:6]), 
                                     scale = "row", cluster_cols = TRUE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                     border_color = NA, main = "TNF stimulation", labels_row = data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,7])
         }
       } else {
         if (nrow(data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,1:6]) == 1) {
           v$heatmap_tnf <- pheatmap((data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,1:6]), 
                                     cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, 
                                     main = "TNF stimulation", labels_row = data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,7])
         } else {
           v$heatmap_tnf <- pheatmap((data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,1:6]), 
                                     cluster_cols = TRUE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                     border_color = NA, main = "TNF stimulation", labels_row = data_tnf_heat[data_tnf_heat$gene %in% v$tnf$gene,7])
         }
       }
     }
   }
   if(input$cell_line == "SGBS"){
     if (!is.null(v$stim)) {
       if (input$mode == "normalized") {
         if (nrow(v$stim) == 1) {
           v$heatmap_tnf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa", arr.ind = TRUE)], 
                                     scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, 
                                     main = "TNFa", labels_row = v$stim$SYMBOL,
                                     labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa"),2])
         } else {
           v$heatmap_tnf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa", arr.ind = TRUE)],
                                     scale = "row", cluster_cols = TRUE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                     border_color = NA, main = "TNFa", labels_row = v$stim$SYMBOL,
                                     labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa"),2])
         }
       } else {
         if (nrow(v$stim) == 1) {
           v$heatmap_tnf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa", arr.ind = TRUE)], 
                                     cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, 
                                     main = "TNFa", labels_row = v$stim$SYMBOL,
                                     labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa"),2])
         } else {
           v$heatmap_tnf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa", arr.ind = TRUE)],
                                     cluster_cols = FALSE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                     border_color = NA, main = "TNFa", labels_row = v$stim$SYMBOL,
                                     labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TNFa"),2])
         }
       }
     }
   }
 })
 
 #Volcano
 observe({
   if(input$cell_line == "hMADS"){
     if(!is.null(v$tnf)){
       v$volcano_tnf <- ggplot(data = data_tnf, aes(x=log2FoldChange, y=-log10(padj), label=gene)) +
         geom_point(color="#EEF2F2") + geom_point(data = v$tnf, aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
         theme_classic() +
         geom_text_repel(data = v$tnf, aes(label=gene)) + #not yet implemented in ggplotly
         geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
     }
   }
   if(input$cell_line == "SGBS"){
     if(!is.null(v$stim)){
       v$volcano_tnf <- ggplot(data = data_stim$TNFa, aes(x=log2FoldChange, y=-log10(padj))) +
         geom_point(color="#EEF2F2") + geom_point(data = data_stim$TNFa[data_stim$TNFa$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
         theme_classic() +
         geom_text_repel(data = data_stim$TNFa[data_stim$TNFa$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
         geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
     }
   }
 })

 #insulin SGBS AND NIGA
 #heatmap
 observe({
   if (input$cell_line2 == "in vivo") {
      if (!is.null(v$ins)) {
         if (input$mode == "normalized") {
           if (nrow(data_ins_heat[rownames(data_ins_heat) %in% v$ins$genes,]) == 1) {
             v$heatmap_ins <- pheatmap((data_ins_heat[rownames(data_ins_heat) %in% v$ins$genes,]), 
                                       scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                       color = farben[[input$gradient_col]], border_color = NA, 
                                       main = "Insulin", labels_row = v$ins$genes)
           } else {
             v$heatmap_ins <- pheatmap((data_ins_heat[rownames(data_ins_heat) %in% v$ins$genes,]), 
                                       scale = "row", cluster_cols = TRUE, 
                                       color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                       border_color = NA, main = "Insulin", labels_row = v$ins$genes)
           }
         } else {
           if (nrow(data_ins_heat[rownames(data_ins_heat) %in% v$ins$genes,]) == 1) {
             v$heatmap_ins <- pheatmap((data_ins_heat[rownames(data_ins_heat) %in% v$ins$genes,]), 
                                       cluster_cols = FALSE, cluster_rows= FALSE, 
                                       color = farben[[input$gradient_col]], border_color = NA, 
                                       main = "Insulin", labels_row = v$ins$genes)
           } else {
             v$heatmap_ins <- pheatmap((data_ins_heat[rownames(data_ins_heat) %in% v$ins$genes,]), 
                                       cluster_cols = FALSE, 
                                       color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                       border_color = NA, main = "Insulin", labels_row = v$ins$genes)
           }
         }
       }
   }
   if (input$cell_line == "SGBS") {
       if (!is.null(v$stim)) {
         if (input$mode == "normalized") {
           if (nrow(v$stim) == 1) {
             v$heatmap_ins <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin", arr.ind = TRUE)], 
                                       scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                       color = farben[[input$gradient_col]], border_color = NA, 
                                       main = "Insulin", labels_row = v$stim$SYMBOL,
                                       labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin"),2])
           } else {
             v$heatmap_ins <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin", arr.ind = TRUE)],
                                       scale = "row", cluster_cols = TRUE, 
                                       color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                       border_color = NA, main = "Insulin", labels_row = v$stim$SYMBOL,
                                       labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin"),2])
           }
         } else {
           if (nrow(v$stim) == 1) {
             v$heatmap_ins <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin", arr.ind = TRUE)], 
                                       cluster_cols = FALSE, cluster_rows= FALSE, 
                                       color = farben[[input$gradient_col]], border_color = NA, 
                                       main = "Insulin", labels_row = v$stim$SYMBOL,
                                       labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin"),2])
           } else {
             v$heatmap_ins <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin", arr.ind = TRUE)],
                                       cluster_cols = FALSE, 
                                       color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                       border_color = NA, main = "Insulin", labels_row = v$stim$SYMBOL,
                                       labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Insulin"),2])
           }
         }
     }
   }
 
 })
 
 #Volcano
 observe({
   if(input$cell_line2 == "in vivo"){
     if(!is.null(v$ins)){
       v$volcano_ins <- ggplot(data = data_ins, aes(x=log2FoldChange, y=-log10(padj), label=genes)) +
         geom_point(color="#EEF2F2") + geom_point(data = v$ins, aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
         theme_classic() +
         geom_text_repel(data = v$ins, aes(label=genes)) + #not yet implemented in ggplotly
         geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
     }
   }
   if(input$cell_line == "SGBS"){
     if(!is.null(v$stim)){
       v$volcano_ins <- ggplot(data = data_stim$Insulin, aes(x=log2FoldChange, y=-log10(padj))) +
         geom_point(color="#EEF2F2") + geom_point(data = data_stim$Insulin[data_stim$Insulin$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
         theme_classic() +
         geom_text_repel(data = data_stim$Insulin[data_stim$Insulin$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
         geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
     }
   }
 })
 
#IL6
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_il6 <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IL6", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "IL6", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IL6"),2])
       } else {
         v$heatmap_il6 <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IL6", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "IL6", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IL6"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_il6 <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IL6", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "IL6", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IL6"),2])
       } else {
         v$heatmap_il6 <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IL6", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "IL6", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IL6"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_il6 <- ggplot(data = data_stim$IL6, aes(x=log2FoldChange, y=-log10(padj), label=SYMBOL)) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$IL6, aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$IL6, aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#TGFB1
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_tgf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "TGFB1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1"),2])
       } else {
         v$heatmap_tgf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "TGFB1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_tgf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "TGFB1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1"),2])
       } else {
         v$heatmap_tgf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "TGFB1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "TGFB1"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_tgf <- ggplot(data = data_stim$TGFB1, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$TGFB1[data_stim$TGFB1$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$TGFB1[data_stim$TGFB1$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#Rosi
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_rosi <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Rosiglitazone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone"),2])
       } else {
         v$heatmap_rosi <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Rosiglitazone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_rosi <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Rosiglitazone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone"),2])
       } else {
         v$heatmap_rosi <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Rosiglitazone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Rosiglitazone"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_rosi <- ggplot(data = data_stim$Rosiglitazone, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Rosiglitazone[data_stim$Rosiglitazone$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Rosiglitazone[data_stim$Rosiglitazone$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#Dexa
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_dexa <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Dexamethasone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone"),2])
       } else {
         v$heatmap_dexa <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Dexamethasone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_dexa <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Dexamethasone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone"),2])
       } else {
         v$heatmap_dexa <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Dexamethasone", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Dexamethasone"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_dexa <- ggplot(data = data_stim$Dexamethasone, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Dexamethasone[data_stim$Dexamethasone$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Dexamethasone[data_stim$Dexamethasone$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#Iso
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_iso <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Isoprenaline", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline"),2])
       } else {
         v$heatmap_iso <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Isoprenaline", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_iso <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Isoprenaline", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline"),2])
       } else {
         v$heatmap_iso <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Isoprenaline", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Isoprenaline"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_iso <- ggplot(data = data_stim$Isoprenaline, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Isoprenaline[data_stim$Isoprenaline$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Isoprenaline[data_stim$Isoprenaline$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#IBMX
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_ibmx <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "IBMX", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX"),2])
       } else {
         v$heatmap_ibmx <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "IBMX", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_ibmx <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "IBMX", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX"),2])
       } else {
         v$heatmap_ibmx <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "IBMX", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IBMX"),2])
       }
     }
   }
 })
#Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_ibmx <- ggplot(data = data_stim$IBMX, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$IBMX[data_stim$IBMX$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$IBMX[data_stim$IBMX$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#Metformin
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_met <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Metformin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin"),2])
       } else {
         v$heatmap_met <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Metformin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_met <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Metformin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin"),2])
       } else {
         v$heatmap_met <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Metformin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Metformin"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_met <- ggplot(data = data_stim$Metformin, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Metformin[data_stim$Metformin$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Metformin[data_stim$Metformin$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#Atorvastatin
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_atorva <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Atorvastatin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin"),2])
       } else {
         v$heatmap_atorva <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Atorvastatin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_atorva <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Atorvastatin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin"),2])
       } else {
         v$heatmap_atorva <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Atorvastatin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Atorvastatin"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_atorva <- ggplot(data = data_stim$Atorvastatin, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Atorvastatin[data_stim$Atorvastatin$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Atorvastatin[data_stim$Atorvastatin$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#SB
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_sb <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "SB203580", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580"),2])
       } else {
         v$heatmap_sb <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "SB203580", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_sb <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "SB203580", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580"),2])
       } else {
         v$heatmap_sb <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "SB203580", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SB203580"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_sb <- ggplot(data = data_stim$SB203580, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$SB203580[data_stim$SB203580$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$SB203580[data_stim$SB203580$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#SP
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_sp <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "SP600125", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125"),2])
       } else {
         v$heatmap_sp <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "SP600125", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_sp <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "SP600125", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125"),2])
       } else {
         v$heatmap_sp <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "SP600125", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "SP600125"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_sp <- ggplot(data = data_stim$SP600125, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$SP600125[data_stim$SP600125$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$SP600125[data_stim$SP600125$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#UO
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_uo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "UO126", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "UO126", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "UO126"),2])
       } else {
         v$heatmap_uo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "UO126", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "UO126", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "UO126"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_uo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "UO126", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "UO126", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "UO126"),2])
       } else {
         v$heatmap_uo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "UO126", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "UO126", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "UO126"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_uo <- ggplot(data = data_stim$UO126, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$UO126[data_stim$UO126$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$UO126[data_stim$UO126$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#leptin
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_lep <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Leptin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin"),2])
       } else {
         v$heatmap_lep <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Leptin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_lep <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Leptin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin"),2])
       } else {
         v$heatmap_lep <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Leptin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Leptin"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_lep <- ggplot(data = data_stim$Leptin, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Leptin[data_stim$Leptin$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Leptin[data_stim$Leptin$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#adiponectin
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_adipo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Adiponectin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin"),2])
       } else {
         v$heatmap_adipo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Adiponectin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_adipo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Adiponectin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin"),2])
       } else {
         v$heatmap_adipo <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Adiponectin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Adiponectin"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_adipo <- ggplot(data = data_stim$Adiponectin, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Adiponectin[data_stim$Adiponectin$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Adiponectin[data_stim$Adiponectin$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#igf1
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_igf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "IGF1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1"),2])
       } else {
         v$heatmap_igf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "IGF1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_igf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "IGF1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1"),2])
       } else {
         v$heatmap_igf <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "IGF1", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "IGF1"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_igf <- ggplot(data = data_stim$IGF1, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$IGF1[data_stim$IGF1$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$IGF1[data_stim$IGF1$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#glucose
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_glucose <- pheatmap(v$stim[,which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Glucose", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose"),2])
       } else {
         v$heatmap_glucose <- pheatmap(v$stim[,which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Glucose", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_glucose <- pheatmap(v$stim[,which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Glucose", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose"),2])
       } else {
         v$heatmap_glucose <- pheatmap(v$stim[,which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Glucose", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control_Glucose" | data_stim_group$group == "Glucose"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_glucose <- ggplot(data = data_stim$Glucose, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Glucose[data_stim$Glucose$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Glucose[data_stim$Glucose$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#lauroyl
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_lauro <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Lauroylcarnitine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine"),2])
       } else {
         v$heatmap_lauro <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Lauroylcarnitine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_lauro <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Lauroylcarnitine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine"),2])
       } else {
         v$heatmap_lauro <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Lauroylcarnitine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Lauroylcarnitine"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_lauro <- ggplot(data = data_stim$Lauroylcarnitine, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Lauroylcarnitine[data_stim$Lauroylcarnitine$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Lauroylcarnitine[data_stim$Lauroylcarnitine$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#decanoyl
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_decano <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Decanoyllcarnithine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine"),2])
       } else {
         v$heatmap_decano <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Decanoyllcarnithine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_decano <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Decanoyllcarnithine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine"),2])
       } else {
         v$heatmap_decano <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Decanoyllcarnithine", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Decanoyllcarnithine"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_decano <- ggplot(data = data_stim$Decanoyllcarnithine, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Decanoyllcarnithine[data_stim$Decanoyllcarnithine$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Decanoyllcarnithine[data_stim$Decanoyllcarnithine$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#retinoic acid
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_ra <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "RetinoicAcid", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid"),2])
       } else {
         v$heatmap_ra <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "RetinoicAcid", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_ra <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "RetinoicAcid", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid"),2])
       } else {
         v$heatmap_ra <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "RetinoicAcid", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "RetinoicAcid"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_ra <- ggplot(data = data_stim$RetinoicAcid, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$RetinoicAcid[data_stim$RetinoicAcid$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$RetinoicAcid[data_stim$RetinoicAcid$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#wortmannin
 #heatmap
 observe({
   if (!is.null(v$stim)) {
     if (input$mode == "normalized") {
       if (nrow(v$stim) == 1) {
         v$heatmap_wort <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin", arr.ind = TRUE)], 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Wortmanin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin"),2])
       } else {
         v$heatmap_wort <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin", arr.ind = TRUE)],
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Wortmanin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin"),2])
       }
     } else {
       if (nrow(v$stim) == 1) {
         v$heatmap_wort <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin", arr.ind = TRUE)], 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "Wortmanin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin"),2])
       } else {
         v$heatmap_wort <- pheatmap(v$stim[,which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin", arr.ind = TRUE)],
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "Wortmanin", labels_row = v$stim$SYMBOL,
                                   labels_col = data_stim_group[which(data_stim_group$group == "Control" | data_stim_group$group == "Wortmanin"),2])
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$stim)){
     v$volcano_wort <- ggplot(data = data_stim$Wortmanin, aes(x=log2FoldChange, y=-log10(padj))) +
       geom_point(color="#EEF2F2") + geom_point(data = data_stim$Wortmanin[data_stim$Wortmanin$SYMBOL %in% v$stim$SYMBOL,], aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = data_stim$Wortmanin[data_stim$Wortmanin$SYMBOL %in% v$stim$SYMBOL,], aes(label=SYMBOL)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#Define conditional UI for knock downs
 output$knockdowns <- renderUI({
   # if specific knock out line selected
   if(input$ko_gene == "GLS") {
     # return normal output with correct description etc. for GLS
     return(list(              includeHTML("htmls/description_gls.html"),
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_gls"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_gls", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_gls", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_gls"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_gls2", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_gls2", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_gls")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_gls.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_gls.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               includeHTML("htmls/reference_gls.html") 
     )
     )
   }
   if(input$ko_gene == "AQP7") {
     # return normal output with correct description etc. for AQP7
     return(list(              includeHTML("htmls/description_aqp7.html"),
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_aqp7"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_aqp71", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_aqp71", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_aqp7"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_aqp72", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_aqp72", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_aqp7")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_aqp7.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_aqp7.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               includeHTML("htmls/reference_aqp7.html") 
     )
     )
   }
      if(input$ko_gene == "CD248") {
     # return normal output with correct description etc. for CD248
     return(list(              includeHTML("htmls/description_cd248.html"),
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_cd248"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_cd2481", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_cd2481", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_cd248"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_cd2482", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_cd2482", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_cd248")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_cd248.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_cd248.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               includeHTML("htmls/reference_cd248.html") 
     )
     )
   }
   if(input$ko_gene == "C14orf180") {
     # return normal output with correct description etc. for C14orf180
     return(list(              includeHTML("htmls/description_C14.html"),
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_c14"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_c141", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_c141", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_c14"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_c142", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_c142", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_c14")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_C14.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_C14.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               includeHTML("htmls/reference_C14.html") 
     )
     )
   } 
   
   else {
     return(list(div("Please select a knock out/down cell line!")))
   }
 })
 #knock downs/ outs
 #gls
 #heatmap
 observe({
   if (!is.null(v$gls)) {
     if (input$mode == "normalized") {
       if (nrow(data_gls_heat[data_gls_heat$gene %in% v$gls$gene,]) == 1) {
         v$heatmap_gls <- pheatmap((data_gls_heat[rownames(data_gls_heat) %in% v$gls$gene,1:6]), 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "gls", labels_row = v$gls$gene)
       } else {
         v$heatmap_gls <- pheatmap((data_gls_heat[data_gls_heat$gene %in% v$gls$gene,1:6]), 
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "gls", labels_row = v$gls$gene)
       }
     } else {
       if (nrow(data_gls_heat[data_gls_heat$gene %in% v$gls$gene,]) == 1) {
         v$heatmap_gls <- pheatmap((data_gls_heat[data_gls_heat$gene %in% v$gls$gene,1:6]), 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "gls", labels_row = v$gls$gene)
       } else {
         v$heatmap_gls <- pheatmap((data_gls_heat[data_gls_heat$gene %in% v$gls$gene,1:6]), 
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "gls", labels_row = v$gls$gene)
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$gls)){
     v$volcano_gls <- ggplot(data = data_gls, aes(x=log2FoldChange, y=-log10(padj), label=gene)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$gls, aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$gls, aes(label=gene)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 }) 
 #aqp7
 #heatmap
 observe({
   if (!is.null(v$aqp7)) {
     if (input$mode == "normalized") {
       if (nrow(data_aqp7_heat[rownames(data_aqp7_heat) %in% v$aqp7$GENE,]) == 1) {
         v$heatmap_aqp7 <- pheatmap((data_aqp7_heat[rownames(data_aqp7_heat) %in% v$aqp7$GENE,]), 
                                   scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "aqp7", labels_row = v$aqp7$gene)
       } else {
         v$heatmap_aqp7 <- pheatmap((data_aqp7_heat[rownames(data_aqp7_heat) %in% v$aqp7$GENE,]), 
                                   scale = "row", cluster_cols = TRUE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "aqp7", labels_row = v$aqp7$gene)
       }
     } else {
       if (nrow(data_aqp7_heat[rownames(data_aqp7_heat) %in% v$aqp7$GENE,]) == 1) {
         v$heatmap_aqp7 <- pheatmap((data_aqp7_heat[rownames(data_aqp7_heat) %in% v$aqp7$GENE,]), 
                                   cluster_cols = FALSE, cluster_rows= FALSE, 
                                   color = farben[[input$gradient_col]], border_color = NA, 
                                   main = "aqp7ulin", labels_row = v$aqp7$gene)
       } else {
         v$heatmap_aqp7 <- pheatmap((data_aqp7_heat[rownames(data_aqp7_heat) %in% v$aqp7$GENE,]), 
                                   cluster_cols = FALSE, 
                                   color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                   border_color = NA, main = "aqp7", labels_row = v$aqp7$gene)
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$aqp7)){
     v$volcano_aqp7 <- ggplot(data = data_aqp7, aes(x=logFC, y=-log10(adj.P.Val), label=GENE)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$aqp7, aes(x=logFC, y=-log10(adj.P.Val)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$aqp7, aes(label=GENE)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 }) 
 #cd248
 #heatmap
 observe({
   if (!is.null(v$cd248)) {
     if (input$mode == "normalized") {
       if (nrow(data_cd248_heat[rownames(data_cd248_heat) %in% v$cd248$GENE,]) == 1) {
         v$heatmap_cd248 <- pheatmap((data_cd248_heat[rownames(data_cd248_heat) %in% v$cd248$GENE,]), 
                                    scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                    color = farben[[input$gradient_col]], border_color = NA, 
                                    main = "cd248", labels_row = v$cd248$gene)
       } else {
         v$heatmap_cd248 <- pheatmap((data_cd248_heat[rownames(data_cd248_heat) %in% v$cd248$GENE,]), 
                                    scale = "row", cluster_cols = TRUE, 
                                    color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                    border_color = NA, main = "cd248", labels_row = v$cd248$gene)
       }
     } else {
       if (nrow(data_cd248_heat[rownames(data_cd248_heat) %in% v$cd248$GENE,]) == 1) {
         v$heatmap_cd248 <- pheatmap((data_cd248_heat[rownames(data_cd248_heat) %in% v$cd248$GENE,]), 
                                    cluster_cols = FALSE, cluster_rows= FALSE, 
                                    color = farben[[input$gradient_col]], border_color = NA, 
                                    main = "cd248ulin", labels_row = v$cd248$gene)
       } else {
         v$heatmap_cd248 <- pheatmap((data_cd248_heat[rownames(data_cd248_heat) %in% v$cd248$GENE,]), 
                                    cluster_cols = FALSE, 
                                    color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                    border_color = NA, main = "cd248", labels_row = v$cd248$gene)
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$cd248)){
     v$volcano_cd248 <- ggplot(data = data_cd248, aes(x=logFC, y=-log10(adj.P.Val), label=GENE)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$cd248, aes(x=logFC, y=-log10(adj.P.Val)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$cd248, aes(label=GENE)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 }) 
 #c14
 #heatmap
 observe({
   if (!is.null(v$c14)) {
     if (input$mode == "normalized") {
       if (nrow(data_c14_heat[rownames(data_c14_heat) %in% v$c14$GENE,]) == 1) {
         v$heatmap_c14 <- pheatmap((data_c14_heat[rownames(data_c14_heat) %in% v$c14$GENE,]), 
                                     scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, 
                                     main = "c14", labels_row = v$c14$gene)
       } else {
         v$heatmap_c14 <- pheatmap((data_c14_heat[rownames(data_c14_heat) %in% v$c14$GENE,]), 
                                     scale = "row", cluster_cols = TRUE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                     border_color = NA, main = "c14", labels_row = v$c14$gene)
       }
     } else {
       if (nrow(data_c14_heat[rownames(data_c14_heat) %in% v$c14$GENE,]) == 1) {
         v$heatmap_c14 <- pheatmap((data_c14_heat[rownames(data_c14_heat) %in% v$c14$GENE,]), 
                                     cluster_cols = FALSE, cluster_rows= FALSE, 
                                     color = farben[[input$gradient_col]], border_color = NA, 
                                     main = "c14ulin", labels_row = v$c14$gene)
       } else {
         v$heatmap_c14 <- pheatmap((data_c14_heat[rownames(data_c14_heat) %in% v$c14$GENE,]), 
                                     cluster_cols = FALSE, 
                                     color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                     border_color = NA, main = "c14", labels_row = v$c14$gene)
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$c14)){
     v$volcano_c14 <- ggplot(data = data_c14, aes(x=logFC, y=-log10(adj.P.Val), label=GENE)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$c14, aes(x=logFC, y=-log10(adj.P.Val)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$c14, aes(label=GENE)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 }) 
 #inflamedAT
 #heatmap
 observe({
   if (!is.null(v$inflamed)) {
     if (input$mode == "normalized") {
       if (nrow(v$inflamed) == 1) {
         v$heatmap_inflAT <- pheatmap(v$inflamed[,inflamed_select()], 
                                  scale = "row", cluster_cols = FALSE, cluster_rows= FALSE, 
                                  color = farben[[input$gradient_col]], border_color = NA, 
                                  main = "Inflamed AT", labels_row = v$inflamed$SYMBOL)
       } else {
         v$heatmap_inflAT <- pheatmap(v$inflamed[,inflamed_select()],
                                  scale = "row", cluster_cols = TRUE, 
                                  color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                  border_color = NA, main = "Inflamed AT", labels_row = v$inflamed$SYMBOL)
       }
     } else {
       if (nrow(v$inflamed) == 1) {
         v$heatmap_inflAT <- pheatmap(v$inflamed[,inflamed_select()], 
                                  cluster_cols = FALSE, cluster_rows= FALSE, 
                                  color = farben[[input$gradient_col]], border_color = NA, 
                                  main = "Inflamed AT", labels_row = v$inflamed$SYMBOL)
       } else {
         v$heatmap_inflAT <- pheatmap(v$inflamed[,inflamed_select()],
                                  cluster_cols = FALSE, 
                                  color = farben[[input$gradient_col]], clustering_method = input$cluster_method,
                                  border_color = NA, main = "Inflamed AT", labels_row = v$inflamed$SYMBOL)
       }
     }
   }
 })
 #Volcano
 observe({
   if(!is.null(v$inflamed)){
     v$volcano_inflAT <- ggplot(data = inflamed_select2(), aes(x=logFC, y=-log10(adj.P.Val))) +
       geom_point(color="#EEF2F2") + geom_point(data = inflamed_select2()[inflamed_select2()$Gene_symbol %in% v$inflamed$SYMBOL,], aes(x=logFC, y=-log10(adj.P.Val)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = inflamed_select2()[inflamed_select2()$Gene_symbol %in% v$inflamed$SYMBOL,], aes(label=Gene_symbol)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
#    
#Outputs
  # Define output for Hypoxia heatmap
  output$heatmap_hyp <- renderPlot({
    if (!is.null(v$heatmap_hyp)) v$heatmap_hyp
  })
  # Define output for Hypoxia Volcano  
  output$volcano_hyp <- renderPlot({
    if (!is.null(v$volcano_hyp)) v$volcano_hyp
  })
  # data table Hypoxia
  output$dt_hyp <- DT::renderDataTable({
    if (!is.null(v$hyp)) {
      datatable(v$hyp, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('P.Value', 'adj.P.Val'), digits = 3) %>%
        formatRound(columns=c('logFC', 'AveExpr',"t","B"), digits=3)
    }
  })
  
  # Define output for tnf heatmap
  output$heatmap_tnf <- renderPlot({
    if (!is.null(v$heatmap_tnf)) v$heatmap_tnf
  })
  # Define output for tnf Volcano  
  output$volcano_tnf <- renderPlot({
    if (!is.null(v$volcano_tnf)) v$volcano_tnf
  })
  # data table tnf
  output$dt_tnf <- DT::renderDataTable({
    if(input$cell_line == "hMADS"){
      if (!is.null(v$tnf)) {
        datatable(v$tnf, extensions = 'Buttons',
                  options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
          formatSignif(columns = c('pvalue', 'padj','bonf'), digits = 3) %>%
          formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
      }
    }
    if(input$cell_line == "SGBS"){
      if (!is.null(v$stim)) {
        datatable(data_stim$TNFa[data_stim$TNFa$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                  options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
          formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
          formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
      }
    }
  })
  
  # Define output for insulin heatmap
  output$heatmap_ins <- renderPlot({
    if (!is.null(v$heatmap_ins)) v$heatmap_ins
  })
  # Define output for insulin Volcano  
  output$volcano_ins <- renderPlot({
    if (!is.null(v$volcano_ins)) v$volcano_ins
  })
  # data table insulin
  output$dt_ins <- DT::renderDataTable({
    if(input$cell_line == "hMADS"){
      if (!is.null(v$ins)) {
        datatable(v$ins, extensions = 'Buttons',
                  options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
          formatSignif(columns = c('pvalue', 'padj'), digits = 3) %>%
          formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE','log2FoldChange_shrunken'), digits=3)
      }
    }
    if(input$cell_line == "SGBS"){
      if (!is.null(v$stim)) {
        datatable(data_stim$Insulin[data_stim$Insulin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                  options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
          formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
          formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
      }
    }
  })
  
  # Define output for gls heatmap
  output$heatmap_gls <- renderPlot({
    if (!is.null(v$heatmap_gls)) v$heatmap_gls
  })
  # Define output for gls Volcano  
  output$volcano_gls <- renderPlot({
    if (!is.null(v$volcano_gls)) v$volcano_gls
  })
  # data table gls
  output$dt_gls <- DT::renderDataTable({
    if (!is.null(v$gls)) {
      datatable(v$gls, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','bonf'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  })
  
  # Define output for aqp7 heatmap
  output$heatmap_aqp7 <- renderPlot({
    if (!is.null(v$heatmap_aqp7)) v$heatmap_aqp7
  })
  # Define output for aqp7 Volcano  
  output$volcano_aqp7 <- renderPlot({
    if (!is.null(v$volcano_aqp7)) v$volcano_aqp7
  })
  # data table aqp7
  output$dt_aqp7 <- DT::renderDataTable({
    if (!is.null(v$aqp7)) {
      datatable(v$aqp7, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('P.Value', 'adj.P.Val'), digits = 3) %>%
        formatRound(columns=c('AveExpr', 'logFC','t','B'), digits=3)
    }
  })
  
  # Define output for cd248 heatmap
  output$heatmap_cd248 <- renderPlot({
    if (!is.null(v$heatmap_cd248)) v$heatmap_cd248
  })
  # Define output for cd248 Volcano  
  output$volcano_cd248 <- renderPlot({
    if (!is.null(v$volcano_cd248)) v$volcano_cd248
  })
  # data table cd248
  output$dt_cd248 <- DT::renderDataTable({
    if (!is.null(v$cd248)) {
      datatable(v$cd248, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('P.Value', 'adj.P.Val'), digits = 3) %>%
        formatRound(columns=c('AveExpr', 'logFC','t','B'), digits=3)
    }
  })
  
  # Define output for c14 heatmap
  output$heatmap_c14 <- renderPlot({
    if (!is.null(v$heatmap_c14)) v$heatmap_c14
  })
  # Define output for c14 Volcano  
  output$volcano_c14 <- renderPlot({
    if (!is.null(v$volcano_c14)) v$volcano_c14
  })
  # data table c14
  output$dt_c14 <- DT::renderDataTable({
    if (!is.null(v$c14)) {
      datatable(v$c14, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('P.Value', 'adj.P.Val'), digits = 3) %>%
        formatRound(columns=c('AveExpr', 'logFC','t','B'), digits=3)
    }
  })
#il6
  # Define output for il6 heatmap
  output$heatmap_il6 <- renderPlot({
    if (!is.null(v$heatmap_il6)) v$heatmap_il6
  })
  # Define output for il6 Volcano  
  output$volcano_il6 <- renderPlot({
    if (!is.null(v$volcano_il6)) v$volcano_il6
  })
  # data table il6
  output$dt_il6 <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$IL6[data_stim$IL6$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  })  
  
  #tgf
  # Define output for tgf heatmap
  output$heatmap_tgf <- renderPlot({
    if (!is.null(v$heatmap_tgf)) v$heatmap_tgf
  })
  # Define output for tgf Volcano  
  output$volcano_tgf <- renderPlot({
    if (!is.null(v$volcano_tgf)) v$volcano_tgf
  })
  # data table tgf
  output$dt_tgf <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$TGFB1[data_stim$TGFB1$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  })  

  #rosi
  # Define output for rosi heatmap
  output$heatmap_rosi <- renderPlot({
    if (!is.null(v$heatmap_rosi)) v$heatmap_rosi
  })
  # Define output for rosi Volcano  
  output$volcano_rosi <- renderPlot({
    if (!is.null(v$volcano_rosi)) v$volcano_rosi
  })
  # data table rosi
  output$dt_rosi <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Rosiglitazone[data_stim$Rosiglitazone$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  })  
  
  #dexa
  # Define output for dexa heatmap
  output$heatmap_dexa <- renderPlot({
    if (!is.null(v$heatmap_dexa)) v$heatmap_dexa
  })
  # Define output for dexa Volcano  
  output$volcano_dexa <- renderPlot({
    if (!is.null(v$volcano_dexa)) v$volcano_dexa
  })
  # data table dexa
  output$dt_dexa <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Dexamethasone[data_stim$Dexamethasone$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #glucose
  # Define output for glucose heatmap
  output$heatmap_glucose <- renderPlot({
    if (!is.null(v$heatmap_glucose)) v$heatmap_glucose
  })
  # Define output for glucose Volcano  
  output$volcano_glucose <- renderPlot({
    if (!is.null(v$volcano_glucose)) v$volcano_glucose
  })
  # data table glucose
  output$dt_glucose <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Glucose[data_stim$Glucose$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #iso
  # Define output for iso heatmap
  output$heatmap_iso <- renderPlot({
    if (!is.null(v$heatmap_iso)) v$heatmap_iso
  })
  # Define output for iso Volcano  
  output$volcano_iso <- renderPlot({
    if (!is.null(v$volcano_iso)) v$volcano_iso
  })
  # data table iso
  output$dt_iso <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Isoprenalin[data_stim$Isoprenalin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #igf
  # Define output for igf heatmap
  output$heatmap_igf <- renderPlot({
    if (!is.null(v$heatmap_igf)) v$heatmap_igf
  })
  # Define output for igf Volcano  
  output$volcano_igf <- renderPlot({
    if (!is.null(v$volcano_igf)) v$volcano_igf
  })
  # data table igf
  output$dt_igf <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$IGF1[data_stim$IGF1$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #ibmx
  # Define output for ibmx heatmap
  output$heatmap_ibmx <- renderPlot({
    if (!is.null(v$heatmap_ibmx)) v$heatmap_ibmx
  })
  # Define output for ibmx Volcano  
  output$volcano_ibmx <- renderPlot({
    if (!is.null(v$volcano_ibmx)) v$volcano_ibmx
  })
  # data table ibmx
  output$dt_ibmx <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$IBMX[data_stim$IBMX$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #lauro
  # Define output for lauro heatmap
  output$heatmap_lauro <- renderPlot({
    if (!is.null(v$heatmap_lauro)) v$heatmap_lauro
  })
  # Define output for lauro Volcano  
  output$volcano_lauro <- renderPlot({
    if (!is.null(v$volcano_lauro)) v$volcano_lauro
  })
  # data table lauro
  output$dt_lauro <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Lauroylcarnitine[data_stim$Lauroylcarnitine$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #met
  # Define output for met heatmap
  output$heatmap_met <- renderPlot({
    if (!is.null(v$heatmap_met)) v$heatmap_met
  })
  # Define output for met Volcano  
  output$volcano_met <- renderPlot({
    if (!is.null(v$volcano_met)) v$volcano_met
  })
  # data table met
  output$dt_met <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Metformin[data_stim$Metformin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #atorva
  # Define output for atorva heatmap
  output$heatmap_atorva <- renderPlot({
    if (!is.null(v$heatmap_atorva)) v$heatmap_atorva
  })
  # Define output for atorva Volcano  
  output$volcano_atorva <- renderPlot({
    if (!is.null(v$volcano_atorva)) v$volcano_atorva
  })
  # data table atorva
  output$dt_atorva <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Atorvastatin[data_stim$Atorvastatin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #adipo
  # Define output for adipo heatmap
  output$heatmap_adipo <- renderPlot({
    if (!is.null(v$heatmap_adipo)) v$heatmap_adipo
  })
  # Define output for adipo Volcano  
  output$volcano_adipo <- renderPlot({
    if (!is.null(v$volcano_adipo)) v$volcano_adipo
  })
  # data table adipo
  output$dt_adipo <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Adiponectin[data_stim$Adiponectin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #lep
  # Define output for lep heatmap
  output$heatmap_lep <- renderPlot({
    if (!is.null(v$heatmap_lep)) v$heatmap_lep
  })
  # Define output for lep Volcano  
  output$volcano_lep <- renderPlot({
    if (!is.null(v$volcano_lep)) v$volcano_lep
  })
  # data table lep
  output$dt_lep <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Leptin[data_stim$Leptin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #sb
  # Define output for sb heatmap
  output$heatmap_sb <- renderPlot({
    if (!is.null(v$heatmap_sb)) v$heatmap_sb
  })
  # Define output for sb Volcano  
  output$volcano_sb <- renderPlot({
    if (!is.null(v$volcano_sb)) v$volcano_sb
  })
  # data table sb
  output$dt_sb <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$SB203580[data_stim$SB203580$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #sp
  # Define output for sp heatmap
  output$heatmap_sp <- renderPlot({
    if (!is.null(v$heatmap_sp)) v$heatmap_sp
  })
  # Define output for sp Volcano  
  output$volcano_sp <- renderPlot({
    if (!is.null(v$volcano_sp)) v$volcano_sp
  })
  # data table sp
  output$dt_sp <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$SP600125[data_stim$SP600125$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #uo
  # Define output for uo heatmap
  output$heatmap_uo <- renderPlot({
    if (!is.null(v$heatmap_uo)) v$heatmap_uo
  })
  # Define output for uo Volcano  
  output$volcano_uo <- renderPlot({
    if (!is.null(v$volcano_uo)) v$volcano_uo
  })
  # data table uo
  output$dt_uo <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$UO126[data_stim$UO126$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  
  #wort
  # Define output for wort heatmap
  output$heatmap_wort <- renderPlot({
    if (!is.null(v$heatmap_wort)) v$heatmap_wort
  })
  # Define output for wort Volcano  
  output$volcano_wort <- renderPlot({
    if (!is.null(v$volcano_wort)) v$volcano_wort
  })
  # data table wort
  output$dt_wort <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Wortmanin[data_stim$Wortmanin$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  }) 
  
  #ra
  # Define output for ra heatmap
  output$heatmap_ra <- renderPlot({
    if (!is.null(v$heatmap_ra)) v$heatmap_ra
  })
  # Define output for ra Volcano  
  output$volcano_ra <- renderPlot({
    if (!is.null(v$volcano_ra)) v$volcano_ra
  })
  # data table ra
  output$dt_ra <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$RetinoicAcid[data_stim$RetinoicAcid$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  })  
  
  #decano
  # Define output for decano heatmap
  output$heatmap_decano <- renderPlot({
    if (!is.null(v$heatmap_decano)) v$heatmap_decano
  })
  # Define output for decano Volcano  
  output$volcano_decano <- renderPlot({
    if (!is.null(v$volcano_decano)) v$volcano_decano
  })
  # data table decano
  output$dt_decano <- DT::renderDataTable({
    if (!is.null(v$stim)) {
      datatable(data_stim$Decanoyllcarnithine[data_stim$Decanoyllcarnithine$SYMBOL %in% v$stim$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','stat'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
    }
  })  
  
  #inflamedAT
  # Define output for decano heatmap
  output$heatmap_inflAT <- renderPlot({
    if (!is.null(v$heatmap_inflAT)) v$heatmap_inflAT
  })
  # Define output for decano Volcano  
  output$volcano_inflAT <- renderPlot({
    if (!is.null(v$volcano_inflAT)) v$volcano_inflAT
  })
  output$dt_inflAT <- DT::renderDataTable({
    if (!is.null(v$inflamed)) {
      datatable(inflamed_select2()[inflamed_select2()$Gene_symbol %in% v$inflamed$SYMBOL,], extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('P.Value', 'adj.P.Val'), digits = 3) %>%
        formatRound(columns=c('AveExpr', 'logFC','t','B'), digits=3)
    }
  })
      
#Download handlers for plots
  output$download_pdf_hyp1 <- downloadHandler(
    filename = function() { paste('Hypoxia_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_hyp, device = "pdf") }
  )

  output$download_pdf_hyp2 <- downloadHandler(
    filename = function() { paste('Hypoxia_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_hyp, device = "pdf") }
  )
  
# Download handlers for RDS objects
  output$download_gg_hyp1 <- downloadHandler(
      filename = function() { paste('Hypoxia_Heat_', Sys.Date(), '.rds', sep='') },
      content = function(file) { saveRDS(v$heatmap_hyp, file) }
    )

  output$download_gg_hyp2 <- downloadHandler(
      filename = function() { paste('Hypoxia_Volcano_', Sys.Date(), '.rds', sep='') },
      content = function(file) { saveRDS(v$volcano_hyp, file) }
    )

#Download handlers for plots
  output$download_pdf_hyp1 <- downloadHandler(
    filename = function() { paste('Hypoxia_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_hyp, device = "pdf") }
  )

  output$download_pdf_hyp2 <- downloadHandler(
    filename = function() { paste('Hypoxia_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_hyp, device = "pdf") }
  )
  
# Download handlers for RDS objects
  output$download_gg_hyp1 <- downloadHandler(
      filename = function() { paste('Hypoxia_Heat_', Sys.Date(), '.rds', sep='') },
      content = function(file) { saveRDS(v$heatmap_hyp, file) }
    )

  output$download_gg_hyp2 <- downloadHandler(
      filename = function() { paste('Hypoxia_Volcano_', Sys.Date(), '.rds', sep='') },
      content = function(file) { saveRDS(v$volcano_hyp, file) }
    )
  
  #Download handlers for plots
  output$download_pdf_tnf1 <- downloadHandler(
    filename = function() { paste('tnf_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_tnf, device = "pdf") }
  )
  
  output$download_pdf_tnf2 <- downloadHandler(
    filename = function() { paste('tnf_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_tnf, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_tnf1 <- downloadHandler(
    filename = function() { paste('tnf_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_tnf, file) }
  )
  
  output$download_gg_tnf2 <- downloadHandler(
    filename = function() { paste('tnf_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_tnf, file) }
  )
  
  #Download handlers for plots
  output$download_pdf_il61 <- downloadHandler(
    filename = function() { paste('il6_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_il6, device = "pdf") }
  )
  
  output$download_pdf_il62 <- downloadHandler(
    filename = function() { paste('il6_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_il6, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_il61 <- downloadHandler(
    filename = function() { paste('il6_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_il6, file) }
  )
  
  output$download_gg_il62 <- downloadHandler(
    filename = function() { paste('il6_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_il6, file) }
  ) 

  #Download handlers for plots
  output$download_pdf_tgf1 <- downloadHandler(
    filename = function() { paste('tgf_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_tgf, device = "pdf") }
  )
  
  output$download_pdf_tgf2 <- downloadHandler(
    filename = function() { paste('tgf_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_tgf, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_tgf1 <- downloadHandler(
    filename = function() { paste('tgf_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_tgf, file) }
  )
  
  output$download_gg_tgf2 <- downloadHandler(
    filename = function() { paste('tgf_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_tgf, file) }
  )  
  
  #Download handlers for plots
  output$download_pdf_inflAT1 <- downloadHandler(
    filename = function() { paste('inflAT_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_inflAT, device = "pdf") }
  )
  
  output$download_pdf_inflAT2 <- downloadHandler(
    filename = function() { paste('inflAT_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_inflAT, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_inflAT1 <- downloadHandler(
    filename = function() { paste('inflAT_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_inflAT, file) }
  )
  
  output$download_gg_inflAT2 <- downloadHandler(
    filename = function() { paste('inflAT_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_inflAT, file) }
  )
  
  #Download handlers for plots
  output$download_pdf_ins1 <- downloadHandler(
    filename = function() { paste('ins_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_ins, device = "pdf") }
  )
  
  output$download_pdf_ins2 <- downloadHandler(
    filename = function() { paste('ins_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_ins, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_ins1 <- downloadHandler(
    filename = function() { paste('ins_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_ins, file) }
  )
  
  output$download_gg_ins2 <- downloadHandler(
    filename = function() { paste('ins_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_ins, file) }
  )  
  
  #Download handlers for plots
  output$download_pdf_igf1 <- downloadHandler(
    filename = function() { paste('igf_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_igf, device = "pdf") }
  )
  
  output$download_pdf_igf2 <- downloadHandler(
    filename = function() { paste('igf_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_igf, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_igf1 <- downloadHandler(
    filename = function() { paste('igf_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_igf, file) }
  )
  
  output$download_gg_igf2 <- downloadHandler(
    filename = function() { paste('igf_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_igf, file) }
  )
  #Download handlers for plots
  output$download_pdf_adipo1 <- downloadHandler(
    filename = function() { paste('adipo_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_adipo, device = "pdf") }
  )
  
  output$download_pdf_adipo2 <- downloadHandler(
    filename = function() { paste('adipo_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_adipo, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_adipo1 <- downloadHandler(
    filename = function() { paste('adipo_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_adipo, file) }
  )
  
  output$download_gg_adipo2 <- downloadHandler(
    filename = function() { paste('adipo_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_adipo, file) }
  )
  #Download handlers for plots
  output$download_pdf_lep1 <- downloadHandler(
    filename = function() { paste('lep_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_lep, device = "pdf") }
  )
  
  output$download_pdf_lep2 <- downloadHandler(
    filename = function() { paste('lep_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_lep, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_lep1 <- downloadHandler(
    filename = function() { paste('lep_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_lep, file) }
  )
  
  output$download_gg_lep2 <- downloadHandler(
    filename = function() { paste('lep_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_lep, file) }
  )      
  #Download handlers for plots
  output$download_pdf_rosi1 <- downloadHandler(
    filename = function() { paste('rosi_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_rosi, device = "pdf") }
  )
  
  output$download_pdf_rosi2 <- downloadHandler(
    filename = function() { paste('rosi_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_rosi, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_rosi1 <- downloadHandler(
    filename = function() { paste('rosi_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_rosi, file) }
  )
  
  output$download_gg_rosi2 <- downloadHandler(
    filename = function() { paste('rosi_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_rosi, file) }
  )
  
  #Download handlers for plots
  output$download_pdf_dexa1 <- downloadHandler(
    filename = function() { paste('dexa_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_dexa, device = "pdf") }
  )
  
  output$download_pdf_dexa2 <- downloadHandler(
    filename = function() { paste('dexa_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_dexa, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_dexa1 <- downloadHandler(
    filename = function() { paste('dexa_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_dexa, file) }
  )
  
  output$download_gg_dexa2 <- downloadHandler(
    filename = function() { paste('dexa_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_dexa, file) }
  )  
  
  #Download handlers for plots
  output$download_pdf_iso1 <- downloadHandler(
    filename = function() { paste('iso_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_iso, device = "pdf") }
  )
  
  output$download_pdf_iso2 <- downloadHandler(
    filename = function() { paste('iso_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_iso, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_iso1 <- downloadHandler(
    filename = function() { paste('iso_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_iso, file) }
  )
  
  output$download_gg_iso2 <- downloadHandler(
    filename = function() { paste('iso_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_iso, file) }
  )  

  #Download handlers for plots
  output$download_pdf_ibmx1 <- downloadHandler(
    filename = function() { paste('ibmx_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_ibmx, device = "pdf") }
  )
  
  output$download_pdf_ibmx2 <- downloadHandler(
    filename = function() { paste('ibmx_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_ibmx, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_ibmx1 <- downloadHandler(
    filename = function() { paste('ibmx_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_ibmx, file) }
  )
  
  output$download_gg_ibmx2 <- downloadHandler(
    filename = function() { paste('ibmx_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_ibmx, file) }
  )   
  #Download handlers for plots
  output$download_pdf_met1 <- downloadHandler(
    filename = function() { paste('met_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_met, device = "pdf") }
  )
  
  output$download_pdf_met2 <- downloadHandler(
    filename = function() { paste('met_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_met, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_met1 <- downloadHandler(
    filename = function() { paste('met_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_met, file) }
  )
  
  output$download_gg_met2 <- downloadHandler(
    filename = function() { paste('met_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_met, file) }
  )

  #Download handlers for plots
  output$download_pdf_atorva1 <- downloadHandler(
    filename = function() { paste('atorva_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_atorva, device = "pdf") }
  )
  
  output$download_pdf_atorva2 <- downloadHandler(
    filename = function() { paste('atorva_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_atorva, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_atorva1 <- downloadHandler(
    filename = function() { paste('atorva_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_atorva, file) }
  )
  
  output$download_gg_atorva2 <- downloadHandler(
    filename = function() { paste('atorva_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_atorva, file) }
  )  
  #Download handlers for plots
  output$download_pdf_sb1 <- downloadHandler(
    filename = function() { paste('sb_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_sb, device = "pdf") }
  )
  
  output$download_pdf_sb2 <- downloadHandler(
    filename = function() { paste('sb_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_sb, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_sb1 <- downloadHandler(
    filename = function() { paste('sb_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_sb, file) }
  )
  
  output$download_gg_sb2 <- downloadHandler(
    filename = function() { paste('sb_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_sb, file) }
  )  
  #Download handlers for plots
  output$download_pdf_sp1 <- downloadHandler(
    filename = function() { paste('sp_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_sp, device = "pdf") }
  )
  
  output$download_pdf_sp2 <- downloadHandler(
    filename = function() { paste('sp_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_sp, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_sp1 <- downloadHandler(
    filename = function() { paste('sp_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_sp, file) }
  )
  
  output$download_gg_sp2 <- downloadHandler(
    filename = function() { paste('sp_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_sp, file) }
  ) 
  
  #Download handlers for plots
  output$download_pdf_uo1 <- downloadHandler(
    filename = function() { paste('uo_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_uo, device = "pdf") }
  )
  
  output$download_pdf_uo2 <- downloadHandler(
    filename = function() { paste('uo_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_uo, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_uo1 <- downloadHandler(
    filename = function() { paste('uo_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_uo, file) }
  )
  
  output$download_gg_uo2 <- downloadHandler(
    filename = function() { paste('uo_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_uo, file) }
  )  
  
  #Download handlers for plots
  output$download_pdf_glucose1 <- downloadHandler(
    filename = function() { paste('glucose_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_glucose, device = "pdf") }
  )
  
  output$download_pdf_glucose2 <- downloadHandler(
    filename = function() { paste('glucose_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_glucose, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_glucose1 <- downloadHandler(
    filename = function() { paste('glucose_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_glucose, file) }
  )
  
  output$download_gg_glucose2 <- downloadHandler(
    filename = function() { paste('glucose_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_glucose, file) }
  )  

  #Download handlers for plots
  output$download_pdf_ra1 <- downloadHandler(
    filename = function() { paste('ra_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_ra, device = "pdf") }
  )
  
  output$download_pdf_ra2 <- downloadHandler(
    filename = function() { paste('ra_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_ra, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_ra1 <- downloadHandler(
    filename = function() { paste('ra_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_ra, file) }
  )
  
  output$download_gg_ra2 <- downloadHandler(
    filename = function() { paste('ra_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_ra, file) }
  )  
  
  #Download handlers for plots
  output$download_pdf_lauro1 <- downloadHandler(
    filename = function() { paste('lauro_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_lauro, device = "pdf") }
  )
  
  output$download_pdf_lauro2 <- downloadHandler(
    filename = function() { paste('lauro_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_lauro, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_lauro1 <- downloadHandler(
    filename = function() { paste('lauro_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_lauro, file) }
  )
  
  output$download_gg_lauro2 <- downloadHandler(
    filename = function() { paste('lauro_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_lauro, file) }
  )  
  
  #Download handlers for plots
  output$download_pdf_decano1 <- downloadHandler(
    filename = function() { paste('decano_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_decano, device = "pdf") }
  )
  
  output$download_pdf_decano2 <- downloadHandler(
    filename = function() { paste('decano_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_decano, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_decano1 <- downloadHandler(
    filename = function() { paste('decano_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_decano, file) }
  )
  
  output$download_gg_decano2 <- downloadHandler(
    filename = function() { paste('decano_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_decano, file) }
  )
  
  #Download handlers for plots
  output$download_pdf_gls1 <- downloadHandler(
    filename = function() { paste('gls_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_gls, device = "pdf") }
  )
  
  output$download_pdf_gls2 <- downloadHandler(
    filename = function() { paste('gls_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_gls, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_gls1 <- downloadHandler(
    filename = function() { paste('gls_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_gls, file) }
  )
  
  output$download_gg_gls2 <- downloadHandler(
    filename = function() { paste('gls_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_gls, file) }
  )    
  
  #Download handlers for plots
  output$download_pdf_aqp71 <- downloadHandler(
    filename = function() { paste('aqp7_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_aqp7, device = "pdf") }
  )
  
  output$download_pdf_aqp72 <- downloadHandler(
    filename = function() { paste('aqp7_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_aqp7, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_aqp71 <- downloadHandler(
    filename = function() { paste('aqp7_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_aqp7, file) }
  )
  
  output$download_gg_aqp72 <- downloadHandler(
    filename = function() { paste('aqp7_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_aqp7, file) }
  )   
  
  #Download handlers for plots
  output$download_pdf_c141 <- downloadHandler(
    filename = function() { paste('c14_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_c14, device = "pdf") }
  )
  
  output$download_pdf_c142 <- downloadHandler(
    filename = function() { paste('c14_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_c14, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_c141 <- downloadHandler(
    filename = function() { paste('c14_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_c14, file) }
  )
  
  output$download_gg_c142 <- downloadHandler(
    filename = function() { paste('c14_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_c14, file) }
  ) 
  
  #Download handlers for plots
  output$download_pdf_cd2481 <- downloadHandler(
    filename = function() { paste('cd248_Heat_', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$heatmap_cd248, device = "pdf") }
  )
  
  output$download_pdf_cd2482 <- downloadHandler(
    filename = function() { paste('cd248_Volcano', Sys.Date(), '.pdf', sep='') },
    content = function(file) { ggsave(file, plot = v$volcano_cd248, device = "pdf") }
  )
  
  # Download handlers for RDS objects
  output$download_gg_cd2481 <- downloadHandler(
    filename = function() { paste('cd248_Heat_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$heatmap_cd248, file) }
  )
  
  output$download_gg_cd2482 <- downloadHandler(
    filename = function() { paste('cd248_Volcano_', Sys.Date(), '.rds', sep='') },
    content = function(file) { saveRDS(v$volcano_cd248, file) }
  )    
 })
