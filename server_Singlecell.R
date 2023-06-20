observeEvent(input$Depots, {
  updateSelectizeInput(session = session,
                       inputId = 'gene_for_featureplot',
                       choices = sort(rownames(get(input$Depots))),
                       selected = NULL,
                       server = TRUE)
  
  # updatePickerInput(session = session,
  #                   inputId = 'gene_for_featureplot',
  #                   selected = NULL,
  #                   label = 'Gene: ',
  #                   choices = sort(rownames(get(input$Depots))))
  
  # updateSelectInput(session, "gene_for_featureplot", choices = sort(rownames(get(input$Depots))), selected = NULL)
  if (input$Depots=="META_all") {
    choices_temp <- c("Lymphoid"="Lymphoid", "Myeloid"="Myeloid", "Vascular"="Vascular", "B cell" = "B")
  } else {
    if (input$Depots=="META_all_pvat") {
      choices_temp <- c("Lymphoid"="Lymphoid", "Myeloid"="Myeloid", "FAPs"="FAPs", "Vascular"="Vascular")
    } else {
      choices_temp <- c("Lymphoid"="Lymphoid", "Myeloid"="Myeloid", "FAPs"="FAPs", "Vascular"="Vascular", "B cell" = "B")
    }
  }
  updatePickerInput(session = session,
                    inputId = 'Subcluster',
                    label = 'Subcluster: ',
                    selected = 'Lymphoid',
                    choices = choices_temp)
})



# updateSelectizeInput(session = session, inputId = 'gene_for_featureplot', choices = gene, server = TRUE)

# output$DimPlot <- renderPlotly({
#   a <- data.frame(get(input$Depots)@reductions$umap@cell.embeddings)
#   a <- cbind(a,data.frame(ident=Idents(get(input$Depots))))
#   a <- ggplot(a, aes(x=UMAP_1, y=UMAP_2, color=ident)) +
#     geom_point(size=0.5) +
#     theme_classic() +
#     # theme(aspect.ratio=1) +
#     theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
#     NoAxes()
#   ggplotly(a,width = 600,height = 500)
# })
# 
# output$FeaturePlot <- renderPlotly({
#   a <- FeaturePlot(object = get(input$Depots),features = input$gene_for_featureplot)[[1]] + 
#     NoAxes()
#   ggplotly(a,width = 600,height = 500)
# })

subcluster_seurat <- reactive({
  if (input$Depots=="META_all") {
    depots <- "all"
  }
  if (input$Depots=="META_all_sc") {
    depots <- "sc"
  }
  if (input$Depots=="META_all_om") {
    depots <- "om"
  }
  if (input$Depots=="META_all_pvat") {
    depots <- "pvat"
  }
  
  paste0('WAT_',depots,"_",input$Subcluster)
  
})



output$DimPlot_low <- renderPlot({
  DimPlot_scCustom(seurat_object = get(input$Depots), repel = T)[[1]] + theme(aspect.ratio=1) + NoAxes()
  # DimPlot(get(input$Depots), label = T) + theme(aspect.ratio=1) + NoAxes()
})
output$ui_DimPlot_low <- renderUI({
  shinycustomloader::withLoader(plotOutput("DimPlot_low"), type="html", loader="dnaspin")
})

output$FeaturePlot_low <- renderPlot({
  # FeaturePlot(object = get(input$Depots),features = input$gene_for_featureplot)[[1]] + theme(aspect.ratio=1) + NoAxes()
  Plot_Density_Custom(seurat_object = get(input$Depots), features = input$gene_for_featureplot) + theme(aspect.ratio=1) + NoAxes()
})
output$ui_FeaturePlot_low <- renderUI({
  shinycustomloader::withLoader(plotOutput("FeaturePlot_low"), type="html", loader="dnaspin")
})

output$ViolinPlot_low <- renderPlotly({
  ggplotly(Stacked_VlnPlot(seurat_object = get(input$Depots), features = input$gene_for_featureplot, x_lab_rotate = TRUE)[[1]])
})
output$ui_ViolinPlot_low <- renderUI({
  shinycustomloader::withLoader(plotlyOutput("ViolinPlot_low"), type="html", loader="dnaspin")
})

output$DT_low <- DT::renderDataTable(server = FALSE,{return(get(paste0(input$Depots,'_marker')))},
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
output$ui_DT_low <- renderUI({
  shinycustomloader::withLoader(DT::dataTableOutput("DT_low"), type="html", loader="dnaspin")
})


output$DimPlot_high <- renderPlot({
  DimPlot_scCustom(seurat_object = get(subcluster_seurat()))[[1]] + theme(aspect.ratio=1) + NoAxes()
  # DimPlot(get(subcluster_seurat()), label = T) + theme(aspect.ratio=1) + NoAxes()
})
output$ui_DimPlot_high <- renderUI({
  shinycustomloader::withLoader(plotOutput("DimPlot_high"), type="html", loader="dnaspin")
})

output$FeaturePlot_high <- renderPlot({
  # FeaturePlot(object = get(subcluster_seurat()),features = input$gene_for_featureplot)[[1]] + theme(aspect.ratio=1) + NoAxes()
  Plot_Density_Custom(seurat_object = get(subcluster_seurat()), features = input$gene_for_featureplot) + theme(aspect.ratio=1) + NoAxes()
})
output$ui_FeaturePlot_high <- renderUI({
  shinycustomloader::withLoader(plotOutput("FeaturePlot_high"), type="html", loader="dnaspin")
})

output$ViolinPlot_high <- renderPlotly({
  ggplotly(Stacked_VlnPlot(seurat_object = get(subcluster_seurat()), features = input$gene_for_featureplot, x_lab_rotate = TRUE)[[1]])
})
output$ui_ViolinPlot_high <- renderUI({
  shinycustomloader::withLoader(plotlyOutput("ViolinPlot_high"), type="html", loader="dnaspin")
})

output$DT_high <- DT::renderDataTable(server = FALSE,{return(get(paste0(subcluster_seurat(),'_marker')))},
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
output$ui_DT_high <- renderUI({
  shinycustomloader::withLoader(DT::dataTableOutput("DT_high"), type="html", loader="dnaspin")
})



# output$ui_FeaturePlot <- renderUI({
#   if (!input$gene_for_featureplot %in% c("Not selected","",NULL)) {
#     output$FeaturePlot <- renderPlot({
#       # FeaturePlot(object = get(input$Depots),features = input$gene_for_featureplot)[[1]] + theme(aspect.ratio=1) + NoAxes()
#       Plot_Density_Custom(seurat_object = get(input$Depots), features = input$gene_for_featureplot) + theme(aspect.ratio=1) + NoAxes()
#     })
#     shinycustomloader::withLoader(plotOutput("FeaturePlot"), type="html", loader="dnaspin")
#   } else {
#     code("Please selece a gene")
#   }
# })
# 
# if (input$gene_for_featureplot=="Not selected") {
#   output$FeaturePlot <- renderPlot({
#     # FeaturePlot(object = get(input$Depots),features = input$gene_for_featureplot)[[1]] + theme(aspect.ratio=1) + NoAxes()
#     Plot_Density_Custom(seurat_object = get(input$Depots), features = input$gene_for_featureplot) + theme(aspect.ratio=1) + NoAxes()
#   })
#   output$ui_FeaturePlot <- renderUI({
#     shinycustomloader::withLoader(plotOutput("FeaturePlot"), type="html", loader="dnaspin")
#   })
# } else {
#   output$ui_FeaturePlot <- renderUI({
#     code("Please selece a gene")
#   })
# }