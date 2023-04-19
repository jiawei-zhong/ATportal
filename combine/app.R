# Load packages ----
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Biobase))
suppressMessages(library(Hmisc))  
suppressMessages(library(grid))
suppressMessages(library(cowplot))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinydashboard))  # for box()
suppressMessages(library(ggplotify))
suppressMessages(library(callr))





r_bg(function(){
  shiny::runApp(
    appDir = "../clinical_heatmap",
    port = 3939,
    launch.browser = F,
    host = "127.0.0.1"
  )}, supervise = TRUE)

r_bg(function(){
  shiny::runApp(
    appDir = "../clinical_forestplot",
    port = 4040,
    launch.browser = F,
    host = "127.0.0.1"
  )}, supervise = TRUE)

r_bg(function(){
  shiny::runApp(
    appDir = "../app1",
    port = 4141,
    launch.browser = F,
    host = "127.0.0.1"
  )}, supervise = TRUE)

r_bg(function(){
  shiny::runApp(
    appDir = "../app2",
    port = 4242,
    launch.browser = F,
    host = "127.0.0.1"
  )}, supervise = TRUE)

ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
                navbarPage(title = "WAT portal",
                           tabPanel(title = "Home", "This part is under progress"),
                           navbarMenu(title = "Meta analysis",
                                      tabPanel(title = "Heatmap",
                                               # titlePanel("Heatmap for gene expression and clinical measures"),
                                               tags$iframe(src="http://127.0.0.1:3939/", height="1200vh", width="100%", frameborder="0", scrolling="yes")),
                                      tabPanel(title = "Forest plot",
                                               # titlePanel("Forest plot for gene expression and clinical measure"),
                                               tags$iframe(src="http://127.0.0.1:4040/", height="1200vh", width="100%", frameborder="0", scrolling="yes"))),
                           tabPanel(title = "Module 2", tags$iframe(src="https://www.google.com/", height="1200vh", width="100%", frameborder="0", scrolling="yes")),
                           tabPanel(title = "Module 3", "Module 3 is under progress"),
                           tabPanel(title = "Module 4", "Module 4 is under progress"),
                           tabPanel(title = "Module 5", "Module 5 is under progress"),
                           tabPanel("Shinyapp1", tags$iframe(src="http://127.0.0.1:4141/", height="1200vh", width="100%", frameborder="0", scrolling="yes")),
                           tabPanel("Shinyapp2", tags$iframe(src="http://127.0.0.1:4242/", height="1200vh", width="100%", frameborder="0", scrolling="yes")),
                           tabPanel(title = "About", "This part is under progress")))

server <- function(input, output, session) {
  
}


# Run the app ----
shinyApp(ui = ui, server = server)