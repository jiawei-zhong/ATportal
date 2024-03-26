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
  v <- reactiveValues(hyp=NULL, tnf=NULL, ins=NULL, gls=NULL)
  
  # Reactive expression to load data based on user input
  observeEvent(input$start, {
    gene_split <- strsplit(input$gene_id, split = "\n")  
    v$hyp <- data_hyp[data_hyp$GENE %in% gene_split[[1]],]
    v$tnf <- data_tnf[data_tnf$gene %in% gene_split[[1]],]
    v$ins <- data_ins[data_ins$genes %in% gene_split[[1]],]
    v$gls <- data_gls[data_gls$gene %in% gene_split[[1]],]
    v$cd248 <- data_CD248[data_CD248$GENE %in% gene_split[[1]],]
    v$aqp7 <- data_AQP7[data_AQP7$GENE %in% gene_split[[1]],]
    v$c14 <- data_C14[data_C14$GENE %in% gene_split[[1]],]
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

 #TNF 
 #heatmap
 observe({
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
 })
 
 #Volcano
 observe({
   if(!is.null(v$tnf)){
     v$volcano_tnf <- ggplot(data = data_tnf, aes(x=log2FoldChange, y=-log10(padj), label=gene)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$tnf, aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$tnf, aes(label=gene)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })

 #insulin 
 #heatmap
 observe({
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
 })
 
 #Volcano
 observe({
   if(!is.null(v$ins)){
     v$volcano_ins <- ggplot(data = data_ins, aes(x=log2FoldChange, y=-log10(padj), label=genes)) +
       geom_point(color="#EEF2F2") + geom_point(data = v$ins, aes(x=log2FoldChange, y=-log10(padj)), color="#E2C744") +
       theme_classic() +
       geom_text_repel(data = v$ins, aes(label=genes)) + #not yet implemented in ggplotly
       geom_hline(yintercept=-log10(0.05), col="darkgrey")+ scale_y_continuous(expand = c(0,0))
   }
 })
 
 
 
 
#Define conditional UI for knock downs
 output$knockdowns <- renderUI({
   # if specific knock out line selected
   if(input$ko_gene == "GLS") {
     # return normal output with correct description etc. for GLS
     return(list(              includeHTML("htmls/description_novo.html"),
                               #div(actionButton("toggle_table2", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                               #div(id = "gene_table2", style = "display:none;",includeHTML("htmls/outputs_FACS.html")), 
                               #using class = "gene_table" for everything broke the code somehow
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
                                          includeHTML("htmls/details_novo.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_novo.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/29116032/", "Acosta et al.; PMID: 29116032"),  "(data) and Zhong et al. (portal).",
                               br(),
                               "Raw data for this data set can be downloaded from the GEO database:", tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", "GSE100795")
     )
     )
   }
   if(input$ko_gene == "AQP7") {
     # return normal output with correct description etc. for AQP7
     return(list(              includeHTML("htmls/description_novo.html"),
                               #div(actionButton("toggle_table2", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                               #div(id = "gene_table2", style = "display:none;",includeHTML("htmls/outputs_FACS.html")), 
                               #using class = "gene_table" for everything broke the code somehow
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_aqp7"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_AQP7", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_AQP7", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_aqp7"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_AQP72", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_AQP72", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_aqp7")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_novo.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_novo.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/29116032/", "Acosta et al.; PMID: 29116032"),  "(data) and Zhong et al. (portal).",
                               br(),
                               "Raw data for this data set can be downloaded from the GEO database:", tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", "GSE100795")
     )
     )
   }
      if(input$ko_gene == "CD248") {
     # return normal output with correct description etc. for CD248
     return(list(              includeHTML("htmls/description_novo.html"),
                               #div(actionButton("toggle_table2", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                               #div(id = "gene_table2", style = "display:none;",includeHTML("htmls/outputs_FACS.html")), 
                               #using class = "gene_table" for everything broke the code somehow
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_cd248"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_CD248", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_CD248", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_cd248"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_CD2482", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_CD2482", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_cd248")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_novo.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_novo.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/29116032/", "Acosta et al.; PMID: 29116032"),  "(data) and Zhong et al. (portal).",
                               br(),
                               "Raw data for this data set can be downloaded from the GEO database:", tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", "GSE100795")
     )
     )
   }
   if(input$ko_gene == "C14orf180") {
     # return normal output with correct description etc. for C14orf180
     return(list(              includeHTML("htmls/description_novo.html"),
                               #div(actionButton("toggle_table2", HTML("&dtrif;"), style="font-size:16px; padding: 0px;"), "Show available output options", class="highlight"),
                               #div(id = "gene_table2", style = "display:none;",includeHTML("htmls/outputs_FACS.html")), 
                               #using class = "gene_table" for everything broke the code somehow
                               br(),
                               tabsetPanel(
                                 tabPanel("Heatmap",
                                          br(),
                                          plotOutput(outputId = "heatmap_c14"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_C14", label = "Download PDF", class = "butt"), 
                                          downloadButton(outputId = "download_gg_C14", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))), 
                                 tabPanel("Volcano",
                                          br(),
                                          shinycustomloader::withLoader(plotOutput(outputId = "volcano_c14"), type="html", loader="dnaspin"),
                                          br(),
                                          downloadButton(outputId = "download_pdf_C142", label = "Download PDF", class = "butt"),
                                          downloadButton(outputId = "download_gg_C142", label = "Download ggplot2", class = "butt"),
                                          tags$head(tags$style(".butt {background-color:  #1a4659;} .butt{color: #E2C744;} .butt{border-color: #E2C744;}"))),
                                 tabPanel("Table",
                                          br(),
                                          DT::dataTableOutput("dt_c14")),
                                 tabPanel("Details",
                                          includeHTML("htmls/details_novo.html")
                                 ),
                                 tabPanel("Interpretation",
                                          includeHTML("htmls/interpretation_novo.html")
                                 )),
                               tags$br(),
                               tags$hr(),
                               "If you want to use this figure in your publication, please cite:", tags$a(href="https://pubmed.ncbi.nlm.nih.gov/29116032/", "Acosta et al.; PMID: 29116032"),  "(data) and Zhong et al. (portal).",
                               br(),
                               "Raw data for this data set can be downloaded from the GEO database:", tags$a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100795", "GSE100795")
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
    if (!is.null(v$tnf)) {
      datatable(v$tnf, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj','bonf'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE'), digits=3)
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
    if (!is.null(v$ins)) {
      datatable(v$ins, extensions = 'Buttons',
                options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')) )%>% 
        formatSignif(columns = c('pvalue', 'padj'), digits = 3) %>%
        formatRound(columns=c('baseMean', 'log2FoldChange','lfcSE','log2FoldChange_shrunken'), digits=3)
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
  
  # Download handlers for plots
#   output$download_pdf1 <- downloadHandler(
#     filename = function() { paste('timecourse', Sys.Date(), '.pdf', sep='') },
#     content = function(file) { ggsave(file, plot = v$plot, device = "pdf") }
#   )
#   
#   output$download_pdf2 <- downloadHandler(
#     filename = function() { paste('timecourse_heat', Sys.Date(), '.pdf', sep='') },
#     content = function(file) { ggsave(file, plot = v$heatmap, device = "pdf") }
#   )
#   
#   output$download_pdf3 <- downloadHandler(
#     filename = function() { paste('SVF_line_', Sys.Date(), '.pdf', sep='') },
#     content = function(file) { ggsave(file, plot = v$plot2, device = "pdf") }
#   )
#   
#   output$download_pdf4 <- downloadHandler(
#     filename = function() { paste('SVF_heat_', Sys.Date(), '.pdf', sep='') },
#     content = function(file) { ggsave(file, plot = v$heatmap2, device = "pdf") }
#   )
#   
#   output$download_pdf7 <- downloadHandler(
#     filename = function() { paste('WAT_fractionation_', Sys.Date(), '.pdf', sep='') },
#     content = function(file) { ggsave(file, plot = v$heatmap_novo, device = "pdf") }
#   )
#   
#   output$download_pdf8 <- downloadHandler(
#     filename = function() { paste('tissue_FANTOM_', Sys.Date(), '.pdf', sep='') },
#     content = function(file) { ggsave(file, plot = v$heatmap_fantom, device = "pdf") }
#   )
#   
#   # Download handlers for RDS objects
#   output$download_gg1 <- downloadHandler(
#     filename = function() { paste('timecourse', Sys.Date(), '.rds', sep='') },
#     content = function(file) { saveRDS(v$plot, file) }
#   )
#   
#   output$download_gg2 <- downloadHandler(
#     filename = function() { paste('timecourse_heat', Sys.Date(), '.rds', sep='') },
#     content = function(file) { saveRDS(v$heatmap, file) }
#   )
#   
#   output$download_gg3 <- downloadHandler(
#     filename = function() { paste('SVF_line_', Sys.Date(), '.rds', sep='') },
#     content = function(file) { saveRDS(v$plot2, file) }
#   )
#   
#   output$download_gg4 <- downloadHandler(
#     filename = function() { paste('SVF_heat_', Sys.Date(), '.rds', sep='') },
#     content = function(file) { saveRDS(v$heatmap2, file) }
#   )
#   
#   output$download_gg7 <- downloadHandler(
#     filename = function() { paste('WAT_fractionation_', Sys.Date(), '.rds', sep='') },
#     content = function(file) { saveRDS(v$heatmap_novo, file) }
#   )
#   
#   output$download_gg8 <- downloadHandler(
#     filename = function() { paste('tissue_FANTOM_', Sys.Date(), '.rds', sep='') },
#     content = function(file) { saveRDS(v$heatmap_fantom, file) }
#   )
 })
