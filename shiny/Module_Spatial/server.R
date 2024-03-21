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
suppressMessages(library(extrafont))
suppressMessages(library(shinycustomloader))


# Dataset Gene Feature violon                  rv_datasetgene
# Dataset Table                                rv_dataset
# Dataset Slide SpatialDimplot                 rv_datasetslide
# Dataset Gene Slide SpatialFeatureplot        rv_datasetgeneslide


shinyServer(function(input, output, session) {
  
  rv_datasetgene <- reactiveValues(dataset = "Jesper_et_al_baseline", gene = "LEP")
  rv_dataset <- reactiveValues(dataset = "Jesper_et_al_baseline")
  rv_datasetslide <- reactiveValues(dataset = "Jesper_et_al_baseline", slide = "S42")
  rv_datasetgeneslide <- reactiveValues(dataset = "Jesper_et_al_baseline", gene = "LEP", slide = "S42")
  
  # rv_datasetgene <- reactiveValues()
  # rv_dataset <- reactiveValues()
  # rv_datasetslide <- reactiveValues()
  # rv_datasetgeneslide <- reactiveValues()
  
  updateSelectizeInput(session = session,
                       inputId = 'gene_for_featureplot_STx',
                       choices = gene_all,
                       selected = "LEP",
                       server = TRUE)
  
  observeEvent(input$Dataset_STx, {
    
    # rv_datasetgene <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx)
    # rv_dataset <- reactiveValues(dataset = input$Dataset_STx)
    # rv_datasetslide <- reactiveValues(dataset = input$Dataset_STx, slide = input$Slide_STx)
    # rv_datasetgeneslide <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx, slide = input$Slide_STx)
    rv_datasetgene$dataset <- input$Dataset_STx
    rv_dataset$dataset <- input$Dataset_STx
    rv_datasetslide$dataset <- input$Dataset_STx
    rv_datasetgeneslide$dataset <- input$Dataset_STx
    
    rv_datasetslide$slide <- ifelse(input$Dataset_STx == "Jesper_et_al_baseline", "S42", "S41")
    rv_datasetgeneslide$slide <- ifelse(input$Dataset_STx == "Jesper_et_al_baseline", "S42", "S41")

    updateSelectInput(session = session,
                      inputId = 'Slide_STx',
                      label = 'Slide: ',
                      selected = ifelse(input$Dataset_STx=="Jesper_et_al_baseline","S42","S41"),
                      choices = names(get(input$Dataset_STx)@images))

  })
  
  
  observeEvent(input$gene_for_featureplot_STx, {
    
    # rv_datasetgene <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx)
    # rv_datasetgeneslide <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx, slide = input$Slide_STx)
    
    if (input$gene_for_featureplot_STx != "") {
      rv_datasetgene$gene <- input$gene_for_featureplot_STx
      rv_datasetgeneslide$gene <- input$gene_for_featureplot_STx
    }
    
    
  })
  
  
  observeEvent(input$Slide_STx, {
    
    # rv_datasetslide <- reactiveValues(dataset = input$Dataset_STx, slide = input$Slide_STx)
    # rv_datasetgeneslide <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx, slide = input$Slide_STx)
    rv_datasetslide$slide <- input$Slide_STx
    rv_datasetgeneslide$slide <- input$Slide_STx

  })
  
  
  # FeaturePlot
  plotdata_FeaturePlot_STx <- reactive({
    p1 <- DimPlot_scCustom(seurat_object = get(rv_datasetgene$dataset))[[1]] + theme(aspect.ratio=1,text = element_text(family = "Arial")) + NoAxes() + NoLegend()
    
    
    p2 <- FeaturePlot_scCustom(seurat_object = get(rv_datasetgene$dataset), features = rv_datasetgene$gene, colors_use = viridis::viridis(n = 10, option = "D"))[[1]] + 
      labs(color=rv_datasetgene$gene) +
      xlim(layer_scales(p1, 1)$x$get_limits()) +
      ylim(layer_scales(p1, 1)$y$get_limits()) +
      theme(aspect.ratio=1,
            text = element_text(family = "Arial"),
            plot.title = element_blank(),
            legend.title = element_text(family = "Arial"),
            legend.text = element_text(family = "Arial"),
            legend.position = c(0.85,0.85)) + NoAxes()
    
    legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(rv_datasetgene$dataset))[[1]] + 
                                     theme(aspect.ratio=1,
                                           text = element_text(family = "Arial"),
                                           legend.position = "bottom",
                                           legend.text = element_text(size = 9, family = "Arial")) + 
                                     guides(fill = guide_legend(nrow = 5)) + NoAxes()))
    
    plot_grid(plotlist = list(plot_grid(plotlist = list(p1,p2),ncol = 2),legend),ncol = 1,rel_heights = c(1,0.2))
  })
  
  output$FeaturePlot_STx <- renderPlot({
    plotdata_FeaturePlot_STx()
  })
  
  output$pdf_FeaturePlot_STx <- downloadHandler(
    filename = function() {
      paste0(Sys.time(), ".pdf")
    },
    content = function(fname) {
      ggsave(filename = fname, plot = plotdata_FeaturePlot_STx(), height = 15, width = 10, units = "in", device = cairo_pdf)
    }
  )
  
  output$ui_FeaturePlot_STx <- renderUI({
    list(plotOutput("FeaturePlot_STx",height = "600px") %>% withLoader(type="html", loader="dnaspin"),
         br(),
         downloadButton(outputId = 'pdf_FeaturePlot_STx', label = "Download pdf", style= "simple", color = "primary", size = "sm")
    )
  })
  
  
  
  # ViolinPlot
  plotdata_ViolinPlot_STx <- reactive({
    
    VlnPlot_scCustom(seurat_object = get(rv_datasetgene$dataset), features = rv_datasetgene$gene) + 
      theme(text = element_text(family = "Arial"),
            axis.title.x = element_blank())  # 隐藏 X 轴标题)
    
  })
  
  output$ViolinPlot_STx <- renderPlot({
    plotdata_ViolinPlot_STx()
  })
  
  output$pdf_ViolinPlot_STx <- downloadHandler(
    filename = function() {
      paste0(Sys.time(), ".pdf")
    },
    content = function(fname) {
      ggsave(filename = fname, plot = plotdata_ViolinPlot_STx(), height = 15, width = 10, units = "in", device = cairo_pdf)
    }
  )
  
  output$ui_ViolinPlot_STx <- renderUI({
    list(plotOutput("ViolinPlot_STx",height = "600px") %>% withLoader(type="html", loader="dnaspin"),
         br(),
         downloadButton(outputId = 'pdf_ViolinPlot_STx', label = "Download pdf", style= "simple", color = "primary", size = "sm")
    )
  })
  
  
  # Marker
  output$DT_STx <- DT::renderDataTable(server = FALSE,{return(get(paste0(rv_dataset$dataset,'_marker')))},
                                       extensions = c('Buttons'),
                                       options = list(scrollX = TRUE,
                                                      pageLength = 10,
                                                      lengthMenu = c(10, 25, 50, 100),
                                                      dom = 'Blfrtip',
                                                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                       )
  )
  
  output$ui_DT_STx <- renderUI({
    DT::dataTableOutput("DT_STx") %>% withLoader(type="html", loader="dnaspin")
  })
  
  
  # SpatialDimPlot
  plotdata_SpatialDimPlot_STx <- reactive({
    legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(rv_datasetslide$dataset))[[1]] +
                                     theme(aspect.ratio=1,
                                           text = element_text(family = "Arial"),
                                           legend.position = "bottom",
                                           legend.text = element_text(size = 9, family = "Arial")) +
                                     guides(fill = guide_legend(nrow = 5)) + NoAxes()))
    
    p1 <- SpatialDimPlot(object = get(rv_datasetslide$dataset), images = rv_datasetslide$slide)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Arial"),legend.position = "none") + NoAxes() +
      scale_fill_manual(values = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
    
    background <- SpatialDimPlot(object = get(rv_datasetslide$dataset), images = rv_datasetslide$slide, alpha = 0)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Arial"),legend.position = "none") + NoAxes() +
      scale_fill_manual(values = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
    
    plot_grid(plotlist = list(plot_grid(plotlist = list(p1,background),ncol = 2),legend),ncol = 1,rel_heights = c(1,0.2))
  })
  
  output$SpatialDimPlot_STx <- renderPlot({
    plotdata_SpatialDimPlot_STx()
  })
  
  output$pdf_SpatialDimPlot_STx <- downloadHandler(
    filename = function() {
      paste0(Sys.time(), ".pdf")
    },
    content = function(fname) {
      ggsave(filename = fname, plot = plotdata_SpatialDimPlot_STx(), height = 15, width = 10, units = "in", device = cairo_pdf)
    }
  )
  
  output$ui_SpatialDimPlot_STx <- renderUI({
    list(plotOutput("SpatialDimPlot_STx",height = "600px") %>% withLoader(type="html", loader="dnaspin"),
         br(),
         downloadButton(outputId = 'pdf_SpatialDimPlot_STx', label = "Download pdf", style= "simple", color = "primary", size = "sm")
    )
  })
  
  
  # SpatialFeaturePlot
  plotdata_SpatialFeaturePlot_STx <- reactive({
    p1 <- SpatialFeaturePlot(object = get(rv_datasetgeneslide$dataset), features = rv_datasetgeneslide$gene, images = rv_datasetgeneslide$slide)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Arial"),legend.position = "none") + NoAxes()
    
    legend <- as_ggplot(get_legend(SpatialFeaturePlot(object = get(rv_datasetgeneslide$dataset), features = rv_datasetgeneslide$gene, images = rv_datasetgeneslide$slide)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Arial"),legend.position = "bottom") + NoAxes()))
    
    background <- SpatialFeaturePlot(object = get(rv_datasetgeneslide$dataset), features = rv_datasetgeneslide$gene, images = rv_datasetgeneslide$slide, alpha = 0)[[1]] + NoLegend() + NoAxes()
    
    plot_grid(plotlist = list(plot_grid(plotlist = list(p1,background),ncol = 2),legend),ncol = 1,rel_heights = c(1,0.2))
  })
  
  output$SpatialFeaturePlot_STx <- renderPlot({
    plotdata_SpatialFeaturePlot_STx()
  })
  
  output$pdf_SpatialFeaturePlot_STx <- downloadHandler(
    filename = function() {
      paste0(Sys.time(), ".pdf")
    },
    content = function(fname) {
      ggsave(filename = fname, plot = plotdata_SpatialFeaturePlot_STx(), height = 15, width = 10, units = "in", device = cairo_pdf)
    }
  )
  
  output$ui_SpatialFeaturePlot_STx <- renderUI({
    list(plotOutput("SpatialFeaturePlot_STx",height = "600px") %>% withLoader(type="html", loader="dnaspin"),
         br(),
         downloadButton(outputId = 'pdf_SpatialFeaturePlot_STx', label = "Download pdf", style= "simple", color = "primary", size = "sm")
    )
  })
  
})




