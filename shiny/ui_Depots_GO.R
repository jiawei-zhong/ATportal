tabPanel(title = 'GESA for biological pathways',
  
  tabsetPanel(
    
    tabPanel('Enrichment', 
             
             sidebarLayout(
               position = 'left',
               
               sidebarPanel(
                 
                 uiOutput('cohort_go'),
                 uiOutput('ref_depot_go'),
                 uiOutput('qry_depot_go'),
                 selectInput(inputId = 'p.adjust_go',
                             label = 'P-value adjustment method.',
                             choices = list('FDR'='fdr','Holm'='holm','Hochberg'='hochberg',
                                            'Hommel'='hommel','Bonferroni'='bonferroni',
                                            'BH'='BH','BY'='BY','None'='none'),
                             selected = 'FDR',
                             multiple = F),
                 pickerInput(inputId = 'cor_method_go',
                             label = 'Correlation method.',
                             selected = 'pearson',
                             choices = list('Pearson' = 'pearson', 'Spearman' = 'spearman'),
                             inline = T,
                             options = list(size = 10, `live-search` = T, `actions-box` = T)),
                 actionButton(inputId = 'SearchButton_go',
                              label = 'Search',
                              style = 'simple',
                              color = 'primary',
                              size = 'sm')
                 
               ),
               
               mainPanel(
                 
                 withLoader(plotOutput(outputId = 'graph_go'), loader='dnaspin'),
                 downloadButton(outputId = 'download_GO_enrch'),
                 dataTableOutput(outputId = 'DT_go')
                 
               )
             )
    ),
    tabPanel('Data Table', 
             sidebarLayout(
               sidebarPanel(width = 3,
                            checkboxInput(inputId='show_DE_go',label='Show only significant pathways.',value=T),
                            sliderInput(inputId = 'P_threshold_go',label = 'P-value cut-off.',
                                        min = 0, max = 1, value = 0.05, step = 0.01)),
               mainPanel(downloadButton(outputId = 'download_DT_go_csv','Download as .csv'), 
                         downloadButton(outputId = 'download_DT_go_xlsx','Download as .xlsx'),
                         dataTableOutput(outputId = 'tab_DT_go'))
             ),
    )
  )
)
