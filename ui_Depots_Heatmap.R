tabPanel(title = 'Heatmap',
  
  sidebarLayout(
    position = 'left',
    sidebarPanel(
      
      selectInput(inputId = 'cohort_ht',
                  label = 'Cohorts to include.',
                  choices = list('EMIF'='EMIF','Krieg et al.'='Krieg','Schleinitz et al.'='Schleinitz',
                                 'Keller et al'='Keller','Mazaki et al'='Mazaki',
                                 'Hoggard et al'='Hoggard'),
                  selected = list('EMIF'='EMIF','Krieg et al.'='Krieg','Schleinitz et al.'='Schleinitz',
                                  'Keller et al'='Keller','Mazaki et al'='Mazaki',
                                  'Hoggard et al'='Hoggard'),
                  multiple = T),
      uiOutput('ref_depot_ht'),
      uiOutput('qry_depot_ht'),
      
      selectInput(inputId = 'disease_status_ht',
                  label = 'Disease statuses',
                  choices = list('Healthy'='Healthy','Obese'='Obese','Cancer'='Cancer'),
                  selected = list('Healthy'='Healthy','Obese'='Obese','Cancer'='Cancer'),
                  multiple = T),
      textAreaInput(inputId = 'gene_lst_ht',
                    label = 'Input list of genes',
                    value = '',
                    width = '300px',
                    height = '300px',
                    placeholder = 'Minimum 2 genes.'),
      pickerInput(inputId = 'cor_method_ht',
                  label = 'Correlation method.',
                  selected = 'pearson',
                  choices = list('Pearson' = 'pearson', 'Spearman' = 'spearman'),
                  inline = T,
                  options = list(size = 10, `live-search` = T, `actions-box` = T)),
      actionButton(inputId = 'SearchButton_ht',
                   label = 'Search',
                   style = 'simple',
                   color = 'primary',
                   size = 'sm')
      
    ),
    mainPanel(
      fluidRow(
        splitLayout(
          cellWidths=c('50%','50%'),
          withLoader(plotOutput(outputId='H1_heatmap_ht'),loader='dnaspin'),
          withLoader(plotOutput(outputId='H2_heatmap_ht'),loader='dnaspin')
        ),
        downloadButton(outputId = 'download_H1_pdf', label = 'Download Corr Mat as .pdf'),
        downloadButton(outputId = 'download_H2_pdf', label = 'Download P Mat as .pdf')
      )
    )
  )
)