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

#### Set width and margin in CSS ####
body.css <- "
#body_div {
width: 1200px; 
margin: 15px auto; 
}
h6 { /* Header 4 */
    font-size: 18px;
    line-height: 1.5;
  color: black;
}
h1 { /* Header 4 */
    /font-family: 'Lobster', cursive;
  font-weight: 250;
  line-height: 1.1;
  color:midnightblue;
}
#PRO_tabPanel_div {
height: 1000px;
}
#ML_tabPanel_div {
height: 3000px;
}
#DE_tabPanel_div {
height: 2000px;
}
#CORR_tabPanel_div {
height: 1000px;
}
.landing-page-box { width:19%; height:32vh; background-color:white;margin-right:1%; 
                    border: 3px solid #DDDDDD; margin-bottom: 2px; float: left; 
                    transition: 0.5s ease; position: relative; object-fit: scale-down;}
                    
.landing-page-box:hover, .landing-page-box-about:hover {-webkit-transform: scale(1.05); 
                          -ms-transform: scale(1.05); transform: scale(1.05); }
.landing-page-icon {width:100%; height:26vh; background-color: white;
                     border: 0 ; position: absolute; object-fit: scale-down; }
.landing-page-box-title {font-size: 20px; text-align:center; color: darkblue;
                          font-weight: bold; background-color: none; width:100%; 
                          max-height: 40px; margin-top: 1vh;}
.landing-page-button { position: absolute; top:0; width: 100%; 
                        height: 100%; opacity: 0;}
"


# timeoutSeconds <- 5
# inactivity <- sprintf("function idleTimer() {
# var t = setTimeout(logout, %s);
# window.onmousemove = resetTimer; // catches mouse movements
# window.onmousedown = resetTimer; // catches mouse movements
# window.onclick = resetTimer;     // catches mouse clicks
# window.onscroll = resetTimer;    // catches scrolling
# window.onkeypress = resetTimer;  //catches keyboard actions
# 
# function logout() {
# Shiny.setInputValue('timeOut', '%ss')
# }
# 
# function resetTimer() {
# clearTimeout(t);
# t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
# }
# }
# idleTimer();", timeoutSeconds*1000, timeoutSeconds, timeoutSeconds*1000)

# js_code <- "$(document).on('shiny:value', function(event) {
#   
#   if (event.target.id === 'PRO.num.of.expressed.lipid') {
#     window.scrollBy(0,400);
#   }
#   if (event.target.id === 'DE.pair.detect') {
#     window.scrollBy(0,600);
#   }
#   if (event.target.id === 'ML.ROC.all') {
#     window.scrollBy(0,400);
#   }
# });"


shinyUI(fluidPage(
  #includeScript("www/scrolldown.js"),
  #tags$script(js_code),
  #### CSS/shinythemes/shinyjs contorl ####
  tags$style(body.css),
  #tags$script(inactivity),
  # tags$head(tags$style(
  #   type="text/css",
  #   "#DE.species.enrichment.kegg.pathview img {max-width: 100%; width: 100%; height: auto}"
  # )),
  theme = shinythemes::shinytheme("flatly"),
  # shinyjs::useShinyjs(),
  
  #### Navbar page ####
  # div(id = 'body_div',
  ## Application title
  # titlePanel("LipidSig: a web-based tool for lipidomic data analysis"),
  
  # NavbarPage
  navbarPage(title = 'WAT portal',
             id = 'narbarpage',
             
             
             ####################################
             #####                          #####
             #####           Home           #####
             #####                          #####
             ####################################
             
             source("ui_Home.R", local = TRUE)$value,  #tabPanel #Homepage
             
             ####################################
             #####                          #####
             #####         Clinical         #####
             #####                          #####
             ####################################
             
             # source("ui_Clinical.R", local = TRUE)$value  #tabPanel #Clinical module
             
             navbarMenu(title = "Clinical",
                        source("ui_Clinical_Heatmap.R", local = TRUE)$value,  #tabPanel #Clinical module # Heatmap
                        source("ui_Clinical_Forest.R", local = TRUE)$value,  #tabPanel #Clinical module # Forest
                        source("ui_Clinical_Scatter.R", local = TRUE)$value,  #tabPanel #Clinical module # Scatter plot
                        source("ui_Clinical_Weightloss.R", local = TRUE)$value,  #tabPanel #Clinical module # Weight Loss
                        source("ui_Clinical_BMI.R", local = TRUE)$value,  #tabPanel #Clinical module # BMI
                        source("ui_Clinical_HOMA.R", local = TRUE)$value,  #tabPanel #Clinical module # HOMA
                        source("ui_Clinical_Sex.R", local = TRUE)$value  #tabPanel #Clinical module # Sex difference
             ),
             
             ####################################
             #####                          #####
             #####          Depots          #####
             #####                          #####
             ####################################
             
             navbarMenu(title = "Depots",
                        source("ui_Depots_Heatmap.R", local = TRUE)$value,
                        source("ui_Depots_Forest.R", local = TRUE)$value,  #tabPanel #Depots module)
                        source("ui_Depots_Boxplot.R", local = TRUE)$value,
                        source("ui_Depots_Volcano.R", local = TRUE)$value,
                        source("ui_Depots_GO.R", local = TRUE)$value
             ),             
             ####################################
             #####                          #####
             #####     Characterization     #####
             #####                          #####
             ####################################
             
             source("ui_Characterization.R", local = TRUE)$value,  #tabPanel #Characterization module
             
             ####################################
             #####                          #####
             #####       Single-cell        #####
             #####                          #####
             ####################################
             
             source("ui_Singlecell.R", local = TRUE)$value,  #tabPanel #Single cell module
             
             ####################################
             #####                          #####
             #####         Spatial          #####
             #####                          #####
             ####################################
             
             source("ui_Spatial.R", local = TRUE)$value,  #tabPanel #Spatial module
             
             ####################################
             #####                          #####
             #####       Perturbation       #####
             #####                          #####
             ####################################
             
             source("ui_Perturbation.R", local = TRUE)$value,  #tabPanel #Perturbation module
             
             ####################################
             #####                          #####
             #####        Ethnicity         #####
             #####                          #####
             ####################################
             
             source("ui_Ethnicity.R", local = TRUE)$value,  #tabPanel #Ethnicity module
             
             ####################################
             #####                          #####
             #####          About           #####
             #####                          #####
             ####################################
             
             source("ui_About.R", local = TRUE)$value  #tabPanel #About page
             
             
             
  ) #navbarPage #narbarpage
  # ) #div #body_div
))








