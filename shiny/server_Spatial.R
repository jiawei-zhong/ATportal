observeEvent(input$Dataset_STx, {
  updateSelectizeInput(session = session,
                       inputId = 'gene_for_featureplot_STx',
                       choices = sort(rownames(get(input$Dataset_STx))),
                       selected = NULL,
                       server = TRUE)
  

  updatePickerInput(session = session,
                    inputId = 'Slide_STx',
                    label = 'Slide: ',
                    # selected = ,
                    choices = names(get(input$Dataset_STx)@images))
  

})




output$DimPlot_STx <- renderPlot({
  DimPlot_scCustom(seurat_object = get(input$Dataset_STx),)[[1]] + theme(aspect.ratio=1) + NoAxes()
  # DimPlot(get(input$Dataset_STx), label = T) + theme(aspect.ratio=1) + NoAxes()
})
output$ui_DimPlot_STx <- renderUI({
  shinycustomloader::withLoader(plotOutput("DimPlot_STx"), type="html", loader="dnaspin")
})

output$FeaturePlot_STx <- renderPlot({
  # FeaturePlot(object = get(input$Dataset_STx),features = input$gene_for_featureplot_STx)[[1]] + theme(aspect.ratio=1) + NoAxes()
  Plot_Density_Custom(seurat_object = get(input$Dataset_STx), features = input$gene_for_featureplot_STx) + theme(aspect.ratio=1) + NoAxes()
})
output$ui_FeaturePlot_STx <- renderUI({
  shinycustomloader::withLoader(plotOutput("FeaturePlot_STx"), type="html", loader="dnaspin")
})

output$ViolinPlot_STx <- renderPlotly({
  ggplotly(Stacked_VlnPlot(seurat_object = get(input$Dataset_STx), features = input$gene_for_featureplot_STx, x_lab_rotate = TRUE)[[1]])
})
output$ui_ViolinPlot_STx <- renderUI({
  shinycustomloader::withLoader(plotlyOutput("ViolinPlot_STx"), type="html", loader="dnaspin")
})

output$DT_STx <- DT::renderDataTable(server = FALSE,{return(get(paste0(input$Dataset_STx,'_marker')))},
                                     extensions = c('Buttons'),
                                     options = list(scrollX = TRUE,
                                                    pageLength = 10,
                                                    lengthMenu = c(10, 25, 50, 100),
                                                    dom = 'Blfrtip',
                                                    buttons = list(
                                                      list(extend = "csv", text = "Download Current Page", filename = "page",
                                                           exportOptions = list(modifier = list(page = "current"))),
                                                      list(extend = "csv", text = "Download Full matrix", filename = "data",
                                                           exportOptions = list(modifier = list(page = "all"))))))
output$ui_DT_STx <- renderUI({
  shinycustomloader::withLoader(DT::dataTableOutput("DT_STx"), type="html", loader="dnaspin")
})


output$SpatialDimPlot_STx <- renderPlot({
  SpatialDimPlot(object = get(input$Dataset_STx), images = input$Slide_STx)[[1]] + theme(aspect.ratio=1) + NoAxes() +
    scale_fill_manual(values = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))
})
output$ui_SpatialDimPlot_STx <- renderUI({
  shinycustomloader::withLoader(plotOutput("SpatialDimPlot_STx"), type="html", loader="dnaspin")
})

output$SpatialFeaturePlot_STx <- renderPlot({
  SpatialFeaturePlot(object = get(input$Dataset_STx), features = input$gene_for_featureplot_STx, images = input$Slide_STx)[[1]] + theme(aspect.ratio=1) + NoAxes() + theme(legend.position = "right")
})
output$ui_SpatialFeaturePlot_STx <- renderUI({
  shinycustomloader::withLoader(plotOutput("SpatialFeaturePlot_STx"), type="html", loader="dnaspin")
})


