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
suppressMessages(library(egg))


shinyServer(function(input, output, session) {

  gene_reactive <- reactiveVal("")

  observe({
    if (!is.null(input$data_from_html) && substr(input$data_from_html,1,5)=="gene:") {
      gene_reactive(gsub(pattern = "gene:",replacement = "",x = input$data_from_html))
    }
  })

  output$ui <- renderUI(
    if (gene_reactive()=="") {
      h1("undefined",align = "center")
    } else {
      if (gene_reactive() %in% rownames(META_all)) {

        output$plot <- renderPlot({


          p2 <- FeaturePlot_scCustom(seurat_object = META_all, features = gene_reactive(), colors_use = viridis::viridis(n = 10, option = "D"))[[1]] +
            labs(color=gene_reactive()) +
            xlim(layer_scales(p1, 1)$x$get_limits()) +
            ylim(layer_scales(p1, 1)$y$get_limits()) +
            theme(aspect.ratio=1,
                  text = element_text(family = "Red Hat Display"),
                  plot.title = element_blank(),
                  legend.title = element_text(family = "Red Hat Display"),
                  legend.text = element_text(family = "Red Hat Display"),
                  legend.position = c(0.85,0.85))

          p3 <- VlnPlot_scCustom(seurat_object = META_all, features = gene_reactive(), pt.size=0) +
            theme(aspect.ratio=1,
                  plot.title = element_blank(),
                  axis.title.x = element_blank(),
                  text = element_text(family = "Red Hat Display")) +
            NoLegend() +
            # scale_fill_manual(values = discrete_color)
            scale_fill_manual(values = discrete_color_gradient(length(levels(META_all))))

          egg::ggarrange(p1, p2, p3, nrow = 1)

        })

        plotOutput("plot",height="400px")

      } else {
        h1(paste0(gene_reactive(),' is not in this module'),align = "center")
      }
    }
  )
})


