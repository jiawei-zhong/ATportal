# Load packages ----

suppressMessages(library(shiny))
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

fluidPage(
  tags$script("
        console.log('Shiny App loaded.');

        // 当Shiny与服务器建立连接时发送消息
        $(document).on('shiny:connected', function(event) {
            window.parent.postMessage('Shiny is ready', '*');
        });

        window.addEventListener('message', function(event) {
            console.log('Received message:', event.data);
            Shiny.setInputValue('data_from_html', event.data);
        });
    "),
  shinycustomloader::withLoader(uiOutput("ui", height = "250px"), type="html", loader="dnaspin")
)

