  # Load packages ----
  suppressMessages(library(shiny))
  suppressMessages(library(shinythemes))
  suppressMessages(library(shinyhelper))
  suppressMessages(library(shinyWidgets))
  suppressMessages(library(shinydashboard))  # for box()
  suppressMessages(library(shinycustomloader))
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
  suppressMessages(library(shinycustomloader))
  
  # Depots Gene                     Feature_low violin_low          rv_depotgene
  # Depots                          DT_low                          rv_depot
  # Gene Subcluster                 Feature_high violin_high        rv_genesubcluster
  # Subcluster                      DT_high                         rv_subcluster
  
  
  shinyServer(function(input, output, session) {
    
    rv_depotgene <- reactiveValues(depot = "META_all", gene = "LEP")
    rv_depot <- reactiveValues(depot = "META_all")
    rv_genesubcluster <- reactiveValues(gene = "LEP", subcluster = "WAT_all_Lymphoid")
    rv_subcluster <- reactiveValues(gene = "LEP", subcluster = "WAT_all_Lymphoid")
    
    # rv_depotgene <- reactiveValues(depot = NULL, gene = NULL)
    # rv_depot <- reactiveValues(depot = NULL)
    # rv_genesubcluster <- reactiveValues(gene = NULL, subcluster = NULL)
    # rv_subcluster <- reactiveValues(gene = NULL, subcluster = NULL)
    
    updateSelectizeInput(session = session,
                         inputId = 'gene_for_featureplot',
                         choices = gene_all,
                         selected = "LEP",
                         server = TRUE)
    
    
    observeEvent(input$Depots, {
      
      
      rv_depotgene$depot <- input$Depots
      rv_depot$depot <- input$Depots
      
      
      # rv_depotgene$depot <- input$Depots
      # rv_depot$depot <- input$Depots
      
      
      
      if (input$Depots=="META_all") {
        depots <- "all" 
        rv_genesubcluster$subcluster <- 'WAT_all_Lymphoid'
        rv_subcluster$subcluster <- 'WAT_all_Lymphoid'
      }
      if (input$Depots=="META_all_sc") {
        depots <- "sc"
        rv_genesubcluster$subcluster <- 'WAT_sc_Lymphoid'
        rv_subcluster$subcluster <- 'WAT_sc_Lymphoid'
      }
      if (input$Depots=="META_all_om") {
        depots <- "om"
        rv_genesubcluster$subcluster <- 'WAT_om_Lymphoid'
        rv_subcluster$subcluster <- 'WAT_om_Lymphoid'
      }
      if (input$Depots=="META_all_pvat") {
        depots <- "pvat"
        rv_genesubcluster$subcluster <- 'WAT_pvat_Lymphoid'
        rv_subcluster$subcluster <- 'WAT_pvat_Lymphoid'
      }
      
      
      # if (input$Depots=="META_all") {
      #   depots <- "all" 
      #   rv_genesubcluster$subcluster <- 'WAT_all_Lymphoid'
      #   rv_subcluster$subcluster <- 'WAT_all_Lymphoid'
      # }
      # if (input$Depots=="META_all_sc") {
      #   depots <- "sc"
      #   rv_genesubcluster$subcluster <- 'WAT_sc_Lymphoid'
      #   rv_subcluster$subcluster <- 'WAT_sc_Lymphoid'
      # }
      # if (input$Depots=="META_all_om") {
      #   depots <- "om"
      #   rv_genesubcluster$subcluster <- 'WAT_om_Lymphoid'
      #   rv_subcluster$subcluster <- 'WAT_om_Lymphoid'
      # }
      # if (input$Depots=="META_all_pvat") {
      #   depots <- "pvat"
      #   rv_genesubcluster$subcluster <- 'WAT_pvat_Lymphoid'
      #   rv_subcluster$subcluster <- 'WAT_pvat_Lymphoid'
      # }
        
      
      
      
      

      if (input$Depots=="META_all") {
        choices_temp <- c("T, NK & NKT"="Lymphoid", "mono. & macro."="Myeloid", "vascular"="Vascular", "B" = "B")
      } else {
        if (input$Depots=="META_all_pvat") {
          choices_temp <- c("T, NK & NKT"="Lymphoid", "mono. & macro."="Myeloid", "FAPs"="FAPs", "vascular"="Vascular")
        } else {
          choices_temp <- c("T, NK & NKT"="Lymphoid", "mono. & macro."="Myeloid", "FAPs"="FAPs", "vascular"="Vascular", "B" = "B")
        }
      }
      updateSelectInput(session = session,
                        inputId = 'Subcluster',
                        label = 'Subcluster: ',
                        selected = 'Lymphoid',
                        choices = choices_temp)
      
    })
    
    
  
    
    
    observeEvent(input$gene_for_featureplot, {
      
      if (input$gene_for_featureplot != "") {
        rv_depotgene$gene <- input$gene_for_featureplot
        rv_genesubcluster$gene <- input$gene_for_featureplot
      }
      
    })
    
    
    
    observeEvent(input$Subcluster, {
      
      
      if (input$Depots=="META_all") {
        rv_genesubcluster$subcluster <- paste0('WAT_all_',input$Subcluster)
        rv_subcluster$subcluster <- paste0('WAT_all_',input$Subcluster)
      }
      if (input$Depots=="META_all_sc") {
        rv_genesubcluster$subcluster <- paste0('WAT_sc_',input$Subcluster)
        rv_subcluster$subcluster <- paste0('WAT_sc_',input$Subcluster)
      }
      if (input$Depots=="META_all_om") {
        rv_genesubcluster$subcluster <- paste0('WAT_om_',input$Subcluster)
        rv_subcluster$subcluster <- paste0('WAT_om_',input$Subcluster)
      }
      if (input$Depots=="META_all_pvat") {
        rv_genesubcluster$subcluster <- paste0('WAT_pvat_',input$Subcluster)
        rv_subcluster$subcluster <- paste0('WAT_pvat_',input$Subcluster)
      }
      
      
      # if (input$Depots=="META_all") {
      #   rv_genesubcluster$subcluster <- paste0('WAT_all_',input$Subcluster)
      #   rv_subcluster$subcluster <- paste0('WAT_all_',input$Subcluster)
      # }
      # if (input$Depots=="META_all_sc") {
      #   rv_genesubcluster$subcluster <- paste0('WAT_sc_',input$Subcluster)
      #   rv_subcluster$subcluster <- paste0('WAT_sc_',input$Subcluster)
      # }
      # if (input$Depots=="META_all_om") {
      #   rv_genesubcluster$subcluster <- paste0('WAT_om_',input$Subcluster)
      #   rv_subcluster$subcluster <- paste0('WAT_om_',input$Subcluster)
      # }
      # if (input$Depots=="META_all_pvat") {
      #   rv_genesubcluster$subcluster <- paste0('WAT_pvat_',input$Subcluster)
      #   rv_subcluster$subcluster <- paste0('WAT_pvat_',input$Subcluster)
      # }
      
      
    })
    
    
    # FeaturePlot low
    plotdata_FeaturePlot_low <- reactive({
      
      p1 <- DimPlot_scCustom(seurat_object = get(rv_depotgene$depot), repel = T)[[1]] + 
        theme(aspect.ratio=1, 
              text = element_text(family = "Red Hat Display")) + 
        scale_color_manual(values = discrete_color) +
        NoAxes() + 
        NoLegend()
      
      p2 <- FeaturePlot_scCustom(seurat_object = get(rv_depotgene$depot), features = rv_depotgene$gene, colors_use = viridis::viridis(n = 10, option = "D"))[[1]] + 
        labs(color=rv_depotgene$gene) +
        xlim(layer_scales(p1, 1)$x$get_limits()) +
        ylim(layer_scales(p1, 1)$y$get_limits()) +
        theme(aspect.ratio=1,
              text = element_text(family = "Red Hat Display"),
              plot.title = element_blank(),
              legend.title = element_text(family = "Red Hat Display"),
              legend.text = element_text(family = "Red Hat Display"),
              legend.position = c(0.85,0.85)) + NoAxes()
      
      legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(rv_depotgene$depot), repel = T)[[1]] + 
                                       theme(aspect.ratio=1,
                                             text = element_text(family = "Red Hat Display"),
                                             legend.position = "bottom",
                                             legend.text = element_text(size = 10, family = "Red Hat Display")) + 
                                       scale_color_manual(values = discrete_color) +
                                       NoAxes() + 
                                       guides(color = guide_legend(nrow = 1, override.aes = list(size=4)))))
      
      plot_grid(plotlist = list(plot_grid(plotlist = list(p1,p2),ncol = 2),legend),ncol = 1,rel_heights = c(1,0.2))
    })
    
    output$FeaturePlot_low <- renderPlot({
      plotdata_FeaturePlot_low()
    })
    
    output$pdf_FeaturePlot_low <- downloadHandler(
      filename = function() {
        paste0(Sys.time(), ".pdf")
      },
      content = function(fname) {
        ggsave(filename = fname, plot = plotdata_FeaturePlot_low(), height = 10, width = 10, units = "in", device = cairo_pdf)
      }
    )
    
    output$ui_FeaturePlot_low <- renderUI({
      list(plotOutput("FeaturePlot_low",height = "500px") %>% withLoader(type="html", loader="dnaspin"),
           br(),
           downloadButton(outputId = 'pdf_FeaturePlot_low', label = "Download pdf", style= "simple", color = "primary", size = "sm")
      )
    })
    
    
    
    # ViolinPlot low
    plotdata_ViolinPlot_low <- reactive({
      
      req(rv_depotgene$depot, rv_depotgene$gene)
      
      VlnPlot_scCustom(seurat_object = get(rv_depotgene$depot), features = rv_depotgene$gene, pt.size=0) + 
        theme(text = element_text(family = "Red Hat Display"),
              axis.title.x = element_blank()) + # 隐藏 X 轴标题)
        scale_fill_manual(values = discrete_color)
      
    })
    
    output$ViolinPlot_low <- renderPlot({
      plotdata_ViolinPlot_low()
    })
    
    output$pdf_ViolinPlot_low <- downloadHandler(
      filename = function() {
        paste0(Sys.time(), ".pdf")
      },
      content = function(fname) {
        ggsave(filename = fname, plot = plotdata_ViolinPlot_low(), height = 10, width = 10, units = "in", device = cairo_pdf)
      }
    )
    
    output$ui_ViolinPlot_low <- renderUI({
      list(plotOutput("ViolinPlot_low",height = "500px") %>% withLoader(type="html", loader="dnaspin"),
           br(),
           downloadButton(outputId = 'pdf_ViolinPlot_low', label = "Download pdf", style= "simple", color = "primary", size = "sm")
      )
    })
    
    
    # Marker low
    output$DT_low <- DT::renderDataTable(server = FALSE,{return(get(paste0(rv_depot$depot,'_marker')))},
                                         extensions = c('Buttons'),
                                         options = list(scrollX = TRUE,
                                                        pageLength = 10,
                                                        lengthMenu = c(10, 25, 50, 100),
                                                        dom = 'Blfrtip',
                                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                         )
    )
    
    output$ui_DT_low <- renderUI({
      DT::dataTableOutput("DT_low") %>% withLoader(type="html", loader="dnaspin")
    })
    
    
    
    
    # FeaturePlot high
    plotdata_FeaturePlot_high <- reactive({
      
      req(rv_genesubcluster$subcluster, rv_genesubcluster$gene)
      
      p1 <- DimPlot_scCustom(seurat_object = get(rv_genesubcluster$subcluster))[[1]] + 
        theme(aspect.ratio=1, 
              text = element_text(family = "Red Hat Display")) + 
        scale_color_manual(values = discrete_color_gradient(length(levels(get(rv_genesubcluster$subcluster))))) +
        NoAxes() + 
        NoLegend()
      
      p2 <- FeaturePlot_scCustom(seurat_object = get(rv_genesubcluster$subcluster), features = rv_genesubcluster$gene, colors_use = viridis::viridis(n = 10, option = "D"))[[1]] + 
        labs(color=rv_genesubcluster$gene) +
        xlim(layer_scales(p1, 1)$x$get_limits()) +
        ylim(layer_scales(p1, 1)$y$get_limits()) +
        theme(aspect.ratio=1,
              text = element_text(family = "Red Hat Display"),
              plot.title = element_blank(),
              legend.title = element_text(family = "Red Hat Display"),
              legend.text = element_text(family = "Red Hat Display"),
              legend.position = c(0.9,0.15)) + NoAxes()
      
      legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(rv_genesubcluster$subcluster), repel = T)[[1]] + 
                                       theme(aspect.ratio=1,
                                             text = element_text(family = "Red Hat Display"),
                                             legend.position = "bottom",
                                             legend.text = element_text(size = 10, family = "Red Hat Display")) + 
                                       scale_color_manual(values = discrete_color_gradient(length(levels(get(rv_genesubcluster$subcluster))))) +
                                       NoAxes() + 
                                       guides(color = guide_legend(nrow = 3, override.aes = list(size=4)))))
      
      plot_grid(plotlist = list(plot_grid(plotlist = list(p1,p2),ncol = 2),legend),ncol = 1,rel_heights = c(1,0.2))
    })
    
    output$FeaturePlot_high <- renderPlot({
      plotdata_FeaturePlot_high()
    })
    
    output$pdf_FeaturePlot_high <- downloadHandler(
      filename = function() {
        paste0(Sys.time(), ".pdf")
      },
      content = function(fname) {
        ggsave(filename = fname, plot = plotdata_FeaturePlot_high(), height = 10, width = 10, units = "in", device = cairo_pdf)
      }
    )
    
    output$ui_FeaturePlot_high <- renderUI({
      list(plotOutput("FeaturePlot_high",height = "500px") %>% withLoader(type="html", loader="dnaspin"),
           br(),
           downloadButton(outputId = 'pdf_FeaturePlot_high', label = "Download pdf", style= "simple", color = "primary", size = "sm")
      )
    })
    
    
    # ViolinPlot high
    plotdata_ViolinPlot_high <- reactive({
      req(rv_genesubcluster$subcluster, rv_genesubcluster$gene)
      
      VlnPlot_scCustom(seurat_object = get(rv_genesubcluster$subcluster), features = rv_genesubcluster$gene, pt.size=0) + 
        theme(text = element_text(family = "Red Hat Display"),
              axis.title.x = element_blank()) + # 隐藏 X 轴标题)
        scale_fill_manual(values = discrete_color_gradient(length(levels(get(rv_genesubcluster$subcluster)))))
    })
    
    output$ViolinPlot_high <- renderPlot({
      plotdata_ViolinPlot_high()
    })
    
    output$pdf_ViolinPlot_high <- downloadHandler(
      filename = function() {
        paste0(Sys.time(), ".pdf")
      },
      content = function(fname) {
        ggsave(filename = fname, plot = plotdata_ViolinPlot_high(), height = 10, width = 10, units = "in", device = cairo_pdf)
      }
    )
    
    output$ui_ViolinPlot_high <- renderUI({
      list(plotOutput("ViolinPlot_high",height = "500px") %>% withLoader(type="html", loader="dnaspin"),
           br(),
           downloadButton(outputId = 'pdf_ViolinPlot_high', label = "Download pdf", style= "simple", color = "primary", size = "sm")
      )
    })
    
    
    # Marker high
    
    output$DT_high <- DT::renderDataTable(server = FALSE,{return(get(paste0(rv_subcluster$subcluster,'_marker')))},
                                          extensions = c('Buttons'),
                                          options = list(scrollX = TRUE,
                                                         pageLength = 10,
                                                         lengthMenu = c(10, 25, 50, 100),
                                                         dom = 'Blfrtip',
                                                         buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                          )
    )
    
    
    output$ui_DT_high <- renderUI({
      DT::dataTableOutput("DT_high") %>% withLoader(type="html", loader="dnaspin")
    })
    
    
    
    
  })
  
  

