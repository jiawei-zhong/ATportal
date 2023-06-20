tabPanel(title = 'Volcano-plot',
  tabsetPanel(
    
    tabPanel(title = 'Volcano',
             
             sidebarLayout(
               position = 'left',
               
               sidebarPanel(
                 uiOutput('cohort_vc'),
                 uiOutput('ref_depot_vc'),
                 uiOutput('qry_depot_vc'),
                 
                 uiOutput('FC_threshold_vc'),
                 uiOutput('P_value_threshold_vc'),
                 uiOutput('selected_genes_vc'),
                 selectInput(inputId = 'p.adjust_vc',
                             label = 'P-value adjustment method.',
                             choices = list('FDR'='fdr','Holm'='holm','Hochberg'='hochberg',
                                            'Hommel'='hommel','Bonferroni'='bonferroni',
                                            'BH'='BH','BY'='BY','None'='none'),
                             selected = 'FDR',
                             multiple = F),
                 selectInput(inputId = 'disease_status_vc',
                             label = 'Disease statuses',
                             choices = list('Healthy'='Healthy','Obese'='Obese','Cancer'='Cancer'),
                             selected = list('Healthy'='Healthy','Obese'='Obese','Cancer'='Cancer'),
                             multiple = T),
                 actionButton(inputId = 'SearchButton_vc',
                              label = 'Search',
                              style = 'simple',
                              color = 'primary',
                              size = 'sm')
               ),
               
               mainPanel(
                 withLoader(plotlyOutput(outputId = 'graph_vc'),loader='dnaspin'),
                 br(),
                 br(),
                 downloadButton('download_volcano_pdf','Download as .pdf')
                 #downloadButton('download_volcano_rdata','Download as .RData')
               )
             )
    ), 
    tabPanel('Data Table',
             sidebarLayout(
               sidebarPanel(width=3,checkboxInput(inputId = 'show_DE_vc',label = 'Only show DE genes.',value = T)),
               mainPanel(downloadButton(outputId = 'download_DT_vc_csv','Download as .csv'), 
                         downloadButton(outputId = 'download_DT_vc_xlsx','Download as .xlsx'),
                         dataTableOutput(outputId = 'tab_DT_vc')))
    )
  )
)