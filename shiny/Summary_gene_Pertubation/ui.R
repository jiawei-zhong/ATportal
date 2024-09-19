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
suppressMessages(library(fresh))
suppressMessages(library(scCustomize))

fluidPage(
  use_googlefont("Red Hat Display"),
  use_theme(create_theme(
    theme = "default",
    bs_vars_font(
      family_sans_serif = "'Red Hat Display', cursive"
    )
  )),
  # tags$head(
  #   tags$style('
  #     .container-fluid { padding-left: 0px; padding-right: 0px; }
  #   ')
  # ),
  tags$head(
    tags$script(src="https://cdnjs.cloudflare.com/ajax/libs/iframe-resizer/4.3.1/iframeResizer.contentWindow.min.js"),
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
        ")
  ),
  # tags$script("
  #       console.log('Shiny App loaded.');
  # 
  #       // 当Shiny与服务器建立连接时发送消息
  #       $(document).on('shiny:connected', function(event) {
  #           window.parent.postMessage('Shiny is ready', '*');
  #       });
  # 
  #       window.addEventListener('message', function(event) {
  #           console.log('Received message:', event.data);
  #           Shiny.setInputValue('data_from_html', event.data);
  #       });
  #   "),
  uiOutput("ui", height = "auto")
)

