tabPanel(title='Boxplot',
  sidebarLayout(
    position = 'left',
    
    sidebarPanel(
      
      selectInput(inputId = 'cohort_bx',
                  label = 'Cohorts to include.',
                  choices = list('EMIF'='EMIF','Krieg et al.'='Krieg','Schleinitz et al.'='Schleinitz',
                                 'Keller et al'='Keller','Mazaki et al'='Mazaki',
                                 'Hoggard et al'='Hoggard'),
                  selected = list('EMIF'='EMIF','Krieg et al.'='Krieg','Schleinitz et al.'='Schleinitz',
                                  'Keller et al'='Keller','Mazaki et al'='Mazaki',
                                  'Hoggard et al'='Hoggard'),
                  multiple = T),
      selectInput(inputId = 'ref_depot_bx',
                  label = 'Reference depot.',
                  choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                  selected = 'SAT Abdomen',
                  multiple = F),
      selectInput(inputId = 'qry_depot_bx',
                  label = 'Query depot.',
                  choices = list('SAT Abdomen'='SAT Abdomen','VAT Omentum'='VAT Omentum'),
                  selected = 'VAT Omentum',
                  multiple = F),
      selectInput(inputId = 'disease_status_bx',
                  label = 'Disease statuses',
                  choices = list('Healthy'='Healthy','Obese'='Obese','Cancer'='Cancer'),
                  selected = list('Healthy'='Healthy','Obese'='Obese','Cancer'='Cancer'),
                  multiple = T),
      textAreaInput(inputId = 'gene_bx',
                    label = 'Input gene name.',
                    value = 'ITLN1',
                    width = '300px',
                    height = '50px',
                    placeholder = 'Only one gene'),
      actionButton(inputId = 'SearchButton_bx',
                   label = 'Search',
                   style = 'simple',
                   color = 'primary',
                   size = 'sm')
    ),
    mainPanel(
      withLoader(plotlyOutput(outputId = 'graph_bx'),loader='dnaspin'),
      downloadButton(outputId = 'download_boxplot_pdf')
    )
  )
)
