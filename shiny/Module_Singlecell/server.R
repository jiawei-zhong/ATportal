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


shinyServer(function(input, output, session) {
  
  
  
  # rv_cohortgene <- reactiveValues(cohort = NULL, gene = NULL)
  # rv_cohort <- reactiveValues(cohort = NULL)
  # rv_genesubcluster <- reactiveValues(gene = NULL, subcluster = NULL)
  # rv_subcluster <- reactiveValues(gene = NULL, subcluster = NULL)
  
  updateSelectizeInput(session = session,
                       inputId = 'Gene_bl',
                       choices = gene_all,
                       selected = "LEP",
                       server = TRUE)
  
  observeEvent(input$Cohort_bl, {
    if (input$Cohort_bl=="Massier_et_al") {
      choices_temp <- c("all" = "all","T, NK & NKT" = "Lymphoid", "mono. & macro." = "Myeloid", "vascular" = "Vascular", "B" = "B", "FAPs subcutaneous"="FAPs_sc", "FAPs omental"="FAPs_om", "FAPs perivascular"="FAPs_pvat")
    } 
    if (input$Cohort_bl=="Hinte_et_al") {
      choices_temp <- c("omAT_LTSS" = "omAT_LTSS", "omAT_MTSS" = "omAT_MTSS", "scAT_LTSS" = "scAT_LTSS", "scAT_NEFA" = "scAT_NEFA")
    }
    if (input$Cohort_bl=="Reinisch_et_al") {
      choices_temp <- c("sc" = "sc", "sc_Adipo" = "sc_Adipo", "sc_APCs" = "sc_APCs", "sc_Immune" = "sc_Immune", "vis" = "vis", "vis_Adipo" = "vis_Adipo", "vis_APCs" = "vis_APCs", "vis_Immune" = "vis_Immune", "vis_Meso" = "vis_Meso")
    }
    # if (input$cohort=="Massier_et_al") {
    #     choices_temp <- c("T, NK & NKT"="Lymphoid", "mono. & macro."="Myeloid", "vascular"="Vascular", "B" = "B")
    # } 
    
    updateSelectInput(session = session,
                      inputId = 'Subcluster_bl',
                      label = 'Subcluster',
                      choices = choices_temp)
  })
  
  
  
  
  # FeaturePlot
  plotdata_FeaturePlot <- eventReactive(input$SearchButton_bl, {
    p1 <- DimPlot_scCustom(seurat_object = get(paste0(input$Cohort_bl,"_",input$Subcluster_bl)), repel = T)[[1]] + 
      theme(aspect.ratio=1, 
            text = element_text(family = "Red Hat Display")) + 
      scale_color_manual(values = discrete_color_gradient(length(levels(get(paste0(input$Cohort_bl,"_",input$Subcluster_bl)))))) +
      NoAxes() + 
      NoLegend()
    
    p2 <- FeaturePlot_scCustom(seurat_object = get(paste0(input$Cohort_bl,"_",input$Subcluster_bl)), features = input$Gene_bl, colors_use = viridis::viridis(n = 10, option = "D"))[[1]] + 
      labs(color=input$Gene_bl) +
      xlim(layer_scales(p1, 1)$x$get_limits()) +
      ylim(layer_scales(p1, 1)$y$get_limits()) +
      theme(aspect.ratio=1,
            text = element_text(family = "Red Hat Display"),
            plot.title = element_blank(),
            legend.title = element_text(family = "Red Hat Display"),
            legend.text = element_text(family = "Red Hat Display"),
            legend.position = c(0.85,0.85)) + 
      NoAxes()
    
    # legend <- as_ggplot(get_legend(DimPlot_scCustom(seurat_object = get(paste0(input$Cohort_bl,"_",input$Subcluster_bl)), repel = T)[[1]] + 
    #                                  theme(aspect.ratio=1,
    #                                        text = element_text(family = "Red Hat Display"),
    #                                        legend.position = "bottom",
    #                                        legend.text = element_text(size = 10, family = "Red Hat Display")) + 
    #                                  scale_color_manual(values = discrete_color_gradient(length(levels(get(paste0(input$Cohort_bl,"_",input$Subcluster_bl)))))) +
    #                                  NoAxes() + 
    #                                  guides(color = guide_legend(nrow = 3, override.aes = list(size=4)))))
    
    plot_grid(plotlist = list(plot_grid(plotlist = list(p1,p2),ncol = 2),get(paste0("legend_",input$Cohort_bl,"_",input$Subcluster_bl))),ncol = 1,rel_heights = c(1,0.2))
  })
  
  
  
  
  
  output$ui_FeaturePlot <- renderUI({
    if (input$SearchButton_bl) {

      output$FeaturePlot <- renderPlot({
        plotdata_FeaturePlot()
      })

      output$pdf_FeaturePlot <- downloadHandler(
        filename = function() {
          paste0(Sys.time(), ".pdf")
        },
        content = function(fname) {
          ggsave(filename = fname, plot = plotdata_FeaturePlot(), height = 10, width = 10, units = "in", device = cairo_pdf)
        }
      )

      return(list(plotOutput("FeaturePlot",height = "500px")
                  # %>% withLoader(type="html", loader="dnaspin")
                  ,
                  br(),
                  downloadButton(outputId = 'pdf_FeaturePlot', label = "Download pdf", style= "simple", color = "primary", size = "sm")
      ))
    } else {
      return(NULL)
    }
  })
  
  
  
  # ViolinPlot
  plotdata_ViolinPlot <- eventReactive(input$SearchButton_bl, {
    
    VlnPlot_scCustom(seurat_object = get(paste0(input$Cohort_bl,"_",input$Subcluster_bl)), features = input$Gene_bl, pt.size=0) + 
      theme(text = element_text(family = "Red Hat Display"),
            axis.title.x = element_blank(), legend.position = "none") + # 隐藏 X 轴标题)
      scale_fill_manual(values = discrete_color_gradient(length(levels(get(paste0(input$Cohort_bl,"_",input$Subcluster_bl))))))
    
  })
  
  
  
  output$ui_ViolinPlot <- renderUI({
    if (input$SearchButton_bl) {
      
      output$ViolinPlot <- renderPlot({
        plotdata_ViolinPlot()
      })
      
      output$pdf_ViolinPlot <- downloadHandler(
        filename = function() {
          paste0(Sys.time(), ".pdf")
        },
        content = function(fname) {
          ggsave(filename = fname, plot = plotdata_ViolinPlot(), height = 10, width = 10, units = "in", device = cairo_pdf)
        }
      )
      
      return(list(plotOutput("ViolinPlot",height = "500px") 
                  # %>% withLoader(type="html", loader="dnaspin")
                  ,
                  br(),
                  downloadButton(outputId = 'pdf_ViolinPlot', label = "Download pdf", style= "simple", color = "primary", size = "sm")
      ))
    } else {
      return(NULL)
    }
  })
  
  
  output$ui_citation <- renderUI({if (input$SearchButton_bl) {
    if (!is.null(plotdata_FeaturePlot())) {
      citation_temp <- citation[citation$Cohort_ID %in% input$Cohort_bl,3:9,drop=F] %>% unique()
      
      citation_list <- list()
      for (i in 1:nrow(citation_temp)) {
        
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
  
  
  # # Marker
  # output$DT <- DT::renderDataTable(server = FALSE,{return(get(paste0(input$Cohort_bl,"_",input$Subcluster_bl,'_marker')))},
  #                                      extensions = c('Buttons'),
  #                                      options = list(scrollX = TRUE,
  #                                                     pageLength = 10,
  #                                                     lengthMenu = c(10, 25, 50, 100),
  #                                                     dom = 'Blfrtip',
  #                                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
  #                                      )
  # )
  # 
  # output$ui_DT <- renderUI({
  #     DT::dataTableOutput("DT") %>% withLoader(type="html", loader="dnaspin")
  # })
  
  
})



