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


# Dataset Feature Feature violon               rv_datasetfeature
# Dataset Table                                rv_dataset
# Dataset Slide SpatialDimplot                 rv_datasetslide
# Dataset Feature Slide SpatialFeatureplot     rv_datasetfeatureslide


shinyServer(function(input, output, session) {
  
  rv_datasetfeature <- reactiveValues(dataset = "Jesper_et_al_baseline", feature = "LEP")
  rv_dataset <- reactiveValues(dataset = "Jesper_et_al_baseline")
  rv_datasetslide <- reactiveValues(dataset = "Jesper_et_al_baseline", slide = "S42")
  rv_datasetfeatureslide <- reactiveValues(dataset = "Jesper_et_al_baseline", feature = "LEP", slide = "S42")
  


  
  # dataset update
  observeEvent(input$Dataset_STx, {
    
    # rv_datasetfeature <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx)
    # rv_dataset <- reactiveValues(dataset = input$Dataset_STx)
    # rv_datasetslide <- reactiveValues(dataset = input$Dataset_STx, slide = input$Slide_STx)
    # rv_datasetfeatureslide <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx, slide = input$Slide_STx)
    rv_datasetfeature$dataset <- input$Dataset_STx
    rv_dataset$dataset <- input$Dataset_STx
    rv_datasetslide$dataset <- input$Dataset_STx
    rv_datasetfeatureslide$dataset <- input$Dataset_STx
    
    rv_datasetslide$slide <- ifelse(input$Dataset_STx == "Jesper_et_al_baseline", "S42", "S41")
    rv_datasetfeatureslide$slide <- ifelse(input$Dataset_STx == "Jesper_et_al_baseline", "S42", "S41")

    updateSelectInput(session = session,
                      inputId = 'Slide_STx',
                      label = 'Slide: ',
                      selected = ifelse(input$Dataset_STx=="Jesper_et_al_baseline","S42","S41"),
                      choices = names(get(input$Dataset_STx)@images))

  })
  
  
  # feature update
  updateSelectizeInput(session = session,
                       inputId = 'gene_for_featureplot_STx',
                       choices = c("None"="",gene_all),
                       selected = "LEP",
                       server = TRUE)

  updateSelectizeInput(session = session,
                       inputId = 'deconvolution_for_featureplot_STx',
                       choices = c("None"="",cell_type_all),
                       selected = "",
                       server = TRUE)

  
  observe({
    
    if (input$gene_for_featureplot_STx != "" && input$gene_for_featureplot_STx != "None") {
      rv_datasetfeature$feature <- input$gene_for_featureplot_STx
      rv_datasetfeatureslide$feature <- input$gene_for_featureplot_STx
      updateSelectizeInput(session = session, inputId = "deconvolution_for_featureplot_STx", selected = "")
      
    } else if (input$deconvolution_for_featureplot_STx != "" && input$deconvolution_for_featureplot_STx != "None") {
      rv_datasetfeature$feature <- input$deconvolution_for_featureplot_STx
      rv_datasetfeatureslide$feature <- input$deconvolution_for_featureplot_STx
      updateSelectizeInput(session = session, inputId = "gene_for_featureplot_STx", selected = "")
      
    } else {
      rv_datasetfeature$feature <- NULL
      rv_datasetfeatureslide$feature <- NULL
    }


  })
  
  
  # 监听选择输入，并根据需要更新其他输入
  observeEvent(input$gene_for_featureplot_STx, {
    if (input$gene_for_featureplot_STx != "") {
      updateSelectizeInput(session, "deconvolution_for_featureplot_STx", selected = "")
    }
  })
  
  observeEvent(input$deconvolution_for_featureplot_STx, {
    if (input$deconvolution_for_featureplot_STx != "") {
      updateSelectizeInput(session, "gene_for_featureplot_STx", selected = "")
    }
  })
  
  
  # slide update
  observeEvent(input$Slide_STx, {
    
    # rv_datasetslide <- reactiveValues(dataset = input$Dataset_STx, slide = input$Slide_STx)
    # rv_datasetfeatureslide <- reactiveValues(dataset = input$Dataset_STx, gene = input$gene_for_featureplot_STx, slide = input$Slide_STx)
    rv_datasetslide$slide <- input$Slide_STx
    rv_datasetfeatureslide$slide <- input$Slide_STx

  })
  
  
  # FeaturePlot
  plotdata_FeaturePlot_STx <- reactive({
    p1 <- DimPlot_scCustom(seurat_object = get(rv_datasetfeature$dataset))[[1]] + 
    theme(aspect.ratio=1,text = element_text(family = "Red Hat Display")) + 
    scale_color_manual(values = discrete_color_gradient(length(levels(get(rv_datasetfeature$dataset))))) +
    NoAxes() + 
    NoLegend()
    
    
    p2 <- FeaturePlot_scCustom(seurat_object = get(rv_datasetfeature$dataset), features = rv_datasetfeature$feature, colors_use = viridis::viridis(n = 10, option = "D"))[[1]] + 
      labs(color=rv_datasetfeature$feature) +
      xlim(layer_scales(p1, 1)$x$get_limits()) +
      ylim(layer_scales(p1, 1)$y$get_limits()) +
      theme(aspect.ratio=1,
            text = element_text(family = "Red Hat Display"),
            plot.title = element_blank(),
            legend.title = element_text(family = "Red Hat Display"),
            legend.text = element_text(family = "Red Hat Display"),
            legend.position = c(0.85,0.85)) + NoAxes()
    
    legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(rv_datasetfeature$dataset))[[1]] + 
                                     theme(aspect.ratio=1,
                                           text = element_text(family = "Red Hat Display"),
                                           legend.position = "bottom",
                                           legend.text = element_text(size = 9, family = "Red Hat Display")) + 
                                     scale_color_manual(values = discrete_color_gradient(length(levels(get(rv_datasetfeature$dataset))))) +
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
    
    VlnPlot_scCustom(seurat_object = get(rv_datasetfeature$dataset), features = rv_datasetfeature$feature, pt.size=0) + 
      scale_fill_manual(values = discrete_color_gradient(length(levels(get(rv_datasetfeature$dataset))))) +
      theme(text = element_text(family = "Red Hat Display"),
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
                                           text = element_text(family = "Red Hat Display"),
                                           legend.position = "bottom",
                                           legend.text = element_text(size = 9, family = "Red Hat Display")) +
                                     scale_color_manual(values = discrete_color_gradient(length(levels(get(rv_datasetslide$dataset))))) +
                                     guides(fill = guide_legend(nrow = 5)) + NoAxes()))
    
    p1 <- SpatialDimPlot(object = get(rv_datasetslide$dataset), images = rv_datasetslide$slide)[[1]] + 
      theme(aspect.ratio=1,
          text = element_text(family = "Red Hat Display"),
          legend.position = "none") + 
      scale_fill_manual(values = discrete_color_gradient(length(levels(get(rv_datasetslide$dataset))))) +
      NoAxes() 
      
    
    background <- SpatialDimPlot(object = get(rv_datasetslide$dataset), images = rv_datasetslide$slide, alpha = 0)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Red Hat Display"),legend.position = "none") + NoAxes() +
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
    p1 <- SpatialFeaturePlot(object = get(rv_datasetfeatureslide$dataset), features = rv_datasetfeatureslide$feature, images = rv_datasetfeatureslide$slide)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Red Hat Display"),legend.position = "none") + NoAxes()
    
    legend <- as_ggplot(get_legend(SpatialFeaturePlot(object = get(rv_datasetfeatureslide$dataset), features = rv_datasetfeatureslide$feature, images = rv_datasetfeatureslide$slide)[[1]] + theme(aspect.ratio=1,text = element_text(family = "Red Hat Display"),legend.position = "bottom") + NoAxes()))
    
    background <- SpatialFeaturePlot(object = get(rv_datasetfeatureslide$dataset), features = rv_datasetfeatureslide$feature, images = rv_datasetfeatureslide$slide, alpha = 0)[[1]] + NoLegend() + NoAxes()
    
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




