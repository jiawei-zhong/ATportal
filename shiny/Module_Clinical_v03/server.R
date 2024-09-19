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
suppressMessages(library(tidyr))
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


# meta package source code: forest.R 3985 has changed from "italic(x)[df]^2, hq," to "chi[df]^2, hq," to compatible Red Hat Display


shinyServer(function(input, output, session) {

  ## setting of color
  output$ui_discrete_col_CB2 <- renderUI({
    if (input$discrete_col=="CB2") {
      # 返回PickerInput
      return(list(selectInput(inputId= "discrete_col_CB2", 
                              label="Select a ColorBrewer 2 colorscheme", 
                              choices = c("Accent","Dark2","Paired","Pastel1","Pastel2","Set1","Set2","Set3","Spectral",
                                          "RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG"), 
                              selected = "Dark2"))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_correlation_pc消失
      return(NULL)
    }
  })

  
  ## phenotype association ----
  
  candidate_traits_pc_RNA <- eventReactive(input$cohort_pc_RNA, {
    # 获取所选队列的所有特征
    all_traits <- lapply(input$cohort_pc_RNA, function(cohort) {
      colnames(get(paste0('cor_list_', cohort, "_pearson"))$r)
    })
    
    # 返回所选队列的特征的并集
    trait_vector[trait_vector %in% unique(unlist(all_traits))]
  })
  
  
  # 创建一个反应式值来跟踪是否显示警告信息
  showWarning_pc <- reactiveVal(FALSE)
  
  
  observeEvent(input$cohort_pc_RNA, {
    
    # 检查是否选择了META_sc或META_om，并且选择了多于一个选项
    if (any(input$cohort_pc_RNA %in% c("META_sc", "META_om")) && length(input$cohort_pc_RNA) > 1) {
      showWarning_pc(TRUE) # 设置显示警告
      # 初始化选择器，但保留警告信息显示
      updatePickerInput(session = session,
                        inputId = 'cohort_pc_RNA',
                        label = NULL,
                        selected = NULL,
                        choices = list(META= c("META sc"="META_sc",
                                               "META om"="META_om"),
                                       subcutaneous=c("Kerr, A. (2020) Baseline"="DEOSHeset_Baseline",               #  sc  phenotype
                                                      "Petrus, P. (2018) Baseline"="POeset_Baseline",                #  sc  phenotype
                                                      "Arner, P. (2018)"="SOWOTeset",                                #  sc  phenotype
                                                      "Keller, M. (2017) sc"="Keller_et_al_sc",                      #  sc  phenotype  sex
                                                      "Arner, E. (2012)"="RIKENeset",                                #  sc  phenotype
                                                      "Stančáková, A. (2012)"="METSIM_204eset",                      #  sc  phenotype
                                                      "Raulerson, C. (2019)"="METSIM_434eset",                       #  sc  phenotype
                                                      "Civelek, M. (2017)"="METSIM_770eset",                         #  sc  phenotype
                                                      "Krieg, L. (2021) sc"="Krieg_et_al_sc",                        #  sc  phenotype  sex
                                                      "Arner, P. (2016) sc"="EMIFeset_sc",                           #  sc  phenotype
                                                      "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",   #  sc  phenotype  sex
                                                      "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",  #  sc  phenotype  sex
                                                      "Winnier, DA. (2015)"="GSE64567eset",                          #  sc  phenotype  sex
                                                      "Nono Nankam, PA. (2020)"="Nankam_et_al",                      #  sc  phenotype
                                                      "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",           #  sc  phenotype  sex
                                                      "Salcedo-Tacuma, D. (2022) sc"="GSE188799eset_sc",             #  sc  phenotype
                                                      "MacLaren, RE. (2010) sc"="GSE15524eset_sc",                   #  sc  phenotype  sex
                                                      "Matualatupauw, JC. (2017)"="GSE87382eset",                    #  sc  phenotype  sex
                                                      "Defour, M. (2020)"="GSE154610eset",                           #  sc  phenotype  sex
                                                      "Johansson, LE. (2012) Baseline"="GSE35411eset_Baseline",      #  sc  phenotype  sex
                                                      "Du Plessis, J. (2015) sc"="GSE58979eset_sc",                  #  sc  phenotype  sex
                                                      "Van Bussel, IPG. (2017)"="GSE84046eset_Baseline",             #  sc  phenotype  sex
                                                      "Grundberg, E. (2012)"="E_TABM_1140",                          #  sc  phenotype
                                                      "Heinonen, S. (2017)"="GSE92405eset",                          #  sc  phenotype  sex
                                                      "Sharma, NK. (2016)"="GSE95674eset"                            #  sc  phenotype  sex
                                                      ),
                                       omental=c("Keller, M. (2017) om"="Keller_et_al_om",                           #  om  phenotype  sex
                                                 "Krieg, L. (2021) om"="Krieg_et_al_om",                             #  om  phenotype  sex
                                                 "Arner, P. (2016) om"="EMIFeset_om",                                #  om  phenotype
                                                 "Barberio, MD. (2019)"="GSE88837eset",                              #  om  phenotype       ethnicity
                                                 "Salcedo-Tacuma, D. (2022) om"="GSE188799eset_om",                  #  om  phenotype
                                                 "MacLaren, RE. (2010) om"="GSE15524eset_om",                         #  om  phenotype  sex
                                                 "Du Plessis, J. (2015) om"="GSE58979eset_om",                       #  om  phenotype  sex
                                                 "Aguilera, CM. (2015)"="GSE9624eset"                                #  om  phenotype
                                                 )
                        ),
                        choicesOpt = list(subtext = c("sc","om",rep("sc",22),rep("om",8))),
                        options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE))
    } else if (length(input$cohort_pc_RNA) == 0) {
      # 如果没有选择任何选项，则根据之前的交互可能需要保持警告显示
      # 不改变showWarning_pc的值，保持当前状态
    } else {
      showWarning_pc(FALSE) # 如果选择符合规则，则不显示警告
    }
    
    updatePickerInput(session = session,
                      inputId = 'trait_pc_RNA',
                      label = NULL,
                      # selected = '',
                      choices = candidate_traits_pc_RNA(),
                      # width = "300px",
                      options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = ifelse(test = length(input$cohort_pc_RNA)>1,yes = 1,no = 100), liveSearch = TRUE),
                      # multiple = FALSE
    )
    
    updateTextAreaInput(session = session,
                        inputId = "gene_pc",
                        label = NULL,
                        value = "",
                        placeholder = ifelse(test = length(input$cohort_pc_RNA)>1, yes = "Only 1 gene", no = "One per line, eg:\nPPARG\nADIPOQ"))
    
  }) # observeEvent
  
  
  # 根据showWarning_pc的值动态渲染警告信息
  output$ui_warning_pc <- renderUI({
    if (showWarning_pc()) {
      tags$div(style = "color: red; margin-top: 5px;",
               "META cannot be selected with other cohorts.")
    }
  })
  
  
  output$ui_correction_pc_RNA <- renderUI({
    if (!all(input$cohort_pc_RNA %in% c("META_sc", "META_om"))) {
      # 返回PickerInput
      return(list(pickerInput(inputId = "adjustment_pc_RNA",
                              label = "Adjustment",
                              choices = NULL,
                              # width = "300px",
                              options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
                              multiple = TRUE))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_correlation_pc消失
      return(NULL)
    }
  })
  
  
  observeEvent(input$trait_pc_RNA, {
    if (!all(input$cohort_pc_RNA %in% c("META_sc", "META_om"))) {
      
      temp <- c("Keller, M. (2017) sc"="Keller_et_al_sc",                            #  sc  phenotype  sex
                "Krieg, L. (2021) sc"="Krieg_et_al_sc",                              #  sc  phenotype  sex
                "Imbert, A. (2022) Baseline"="GSE141221_diogenes1_Baseline",         #  sc  phenotype  sex
                "Armenise, C. (2017) Baseline"="GSE95640_diogenes2_Baseline",        #  sc  phenotype  sex
                "Winnier, DA. (2015)"="GSE64567eset",                                #  sc  phenotype  sex
                "Vink, RG. (2017) Baseline"="GSE77962eset_Baseline",                 #  sc  phenotype  sex
                "Keller, M. (2017) om"="Keller_et_al_om",                            #  om  phenotype  sex
                "Krieg, L. (2021) om"="Krieg_et_al_om"                               #  om  phenotype  sex
      )
      
      if (any(input$cohort_pc_RNA %in% temp)) {
        updata_choices <- c("sex"="Gender", trait_vector[trait_vector %in% setdiff(candidate_traits_pc_RNA(),input$trait_pc_RNA)])
      } else {
        updata_choices <- trait_vector[trait_vector %in% setdiff(candidate_traits_pc_RNA(),input$trait_pc_RNA)]
      }
      
      updatePickerInput(session = session,
                        inputId = 'adjustment_pc_RNA',
                        label = NULL,
                        # selected = '',
                        choices = updata_choices,
                        # width = "300px",
                        options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
                        # multiple = TRUE
      )
    }
  })
  
  
  candidate_traits_pc_protein <- reactive({
    trait_vector[trait_vector %in% colnames(pData(Proteomicseset))[colSums(is.na(pData(Proteomicseset)))<nrow(pData(Proteomicseset))]]
  })
  
  
  # output$ui_trait_pc_RNA <- renderUI({
  #   if (input$proteomics_pc) {
  #     # 返回PickerInput
  #     return(list(hr(),
  #                 div("Trait (Proteomics)", class="highlight"),
  #                 br(),
  #                 pickerInput(inputId = "trait_pc_RNA",
  #                             label = NULL,
  #                             choices = candidate_traits_pc_protein(), 
  #                             # width = "300px",
  #                             options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
  #                             multiple = TRUE))
  #     )
  #   } else {
  #     # 返回NULL使ui_correction_pc_protein消失
  #     return(NULL)
  #   }
  # })
  
  
  # output$ui_correction_pc_protein <- renderUI({
  #   if (input$proteomics_pc) {
  #     # 返回PickerInput
  #     return(list(hr(),
  #                 div("Adjustment (Proteomics, optional)", class="highlight"),
  #                 br(),
  #                 pickerInput(inputId = "adjustment_pc_RNA",
  #                             label = NULL,
  #                             choices = NULL,
  #                             # width = "300px",
  #                             options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
  #                             multiple = TRUE))
  #     )
  #   } else {
  #     # 返回NULL使ui_correction_pc_protein消失
  #     return(NULL)
  #   }
  # })
  
  
  # observeEvent(input$trait_pc_RNA, {
  #   if (input$proteomics_pc) {
  #     updatePickerInput(session = session,
  #                       inputId = 'adjustment_pc_RNA',
  #                       label = NULL,
  #                       # selected = '',
  #                       choices = c("sex"="Gender",trait_vector[trait_vector %in% setdiff(candidate_traits_pc_protein(), input$trait_pc_RNA)]),
  #                       # width = "300px",
  #                       options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
  #                       # multiple = TRUE
  #     )
  #   }
  # })
  
  
  plotdata_pc <- eventReactive(input$SearchButton_pc, {
    
    # input$cohort_pc_RNA
    # input$trait_pc_RNA
    gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_pc),split='\n|\t|,')) %>% unique()
    # print(paste0('gene is ', gene))
    
    # check if gene is in cohorts
    temp <- c()
    for (i in input$cohort_pc_RNA) {
      temp <- c(temp,rownames(get(paste0('cor_list_', i, '_', input$correlation_method_pc))$r))
    }
    temp <- unique(temp)

    gene <- intersect(gene, temp)
    
    out <- list(has_plot_RNA = FALSE, plot_RNA = NULL, plottype_RNA = NULL, df_RNA = NULL, cohort_RNA = NULL, has_plot_protein = FALSE, plot_protein = NULL, plottype_protein = NULL, df_protein = NULL) # 初始化返回列表
    
    if (length(gene) > 0 & length(input$cohort_pc_RNA) > 0) {
      if (length(input$cohort_pc_RNA)==1) {
        # gene <- intersect(gene,rownames(get(paste0('cor_list_', input$cohort_pc_RNA, '_', input$correlation_method_pc))$r))
        
        if (!all(input$cohort_pc_RNA %in% c("META_sc", "META_om")) & length(input$adjustment_pc_RNA) > 0) {
          
          cor_df <- list()
          cor_df$r <- matrix(NA, nrow = length(gene), ncol = length(input$trait_pc_RNA), dimnames = list(gene, input$trait_pc_RNA))
          cor_df$p <- cor_df$r
          cor_df$n <- cor_df$r
          
          # 使用expand.grid生成所有基因和特征的组合
          combinations <- expand.grid(gene = gene, trait = input$trait_pc_RNA)
          
          # 对每一对组合应用函数
          results <- apply(combinations, 1, function(x) {
            genetemp <- x['gene']
            trait <- x['trait']
            temp <- t(rbind(Biobase::exprs(get(input$cohort_pc_RNA))[genetemp,,drop=F],t(pData(get(input$cohort_pc_RNA))[,c(trait,input$adjustment_pc_RNA)]))) %>% na.omit()
            
            # 如果存在Gender列，将其转换为数值型
            if ("Gender" %in% colnames(temp)) {
              temp[, "Gender"] <- ifelse(temp[, "Gender"] == "m", 1, 0)
            }
            
            # 保存行名
            rownames_temp <- rownames(temp)
            
            # 将temp矩阵的每一列转换为数值型
            temp <- apply(temp, 2, as.numeric)
            
            # 重新设置行名
            rownames(temp) <- rownames_temp
            
            cor_df_temp_temp <- pcor.test(x = temp[, genetemp], y = temp[, trait], z = temp[, input$adjustment_pc_RNA], method = input$correlation_method_pc)
            c(r = cor_df_temp_temp$estimate, p = cor_df_temp_temp$p.value, n = cor_df_temp_temp$n)
          })
          
          cor_df$r[] <- results['r',]
          cor_df$p[] <- results['p',]
          cor_df$n[] <- results['n',]
          
        } else {
          cor_df <- list()
          cor_df$r <- get(paste0('cor_list_', input$cohort_pc_RNA, '_', input$correlation_method_pc))$r[gene,input$trait_pc_RNA,drop=F]
          cor_df$p <- get(paste0('cor_list_', input$cohort_pc_RNA, '_', input$correlation_method_pc))$p[gene,input$trait_pc_RNA,drop=F]
          cor_df$n <- get(paste0('cor_list_', input$cohort_pc_RNA, '_', input$correlation_method_pc))$n[gene,input$trait_pc_RNA,drop=F]
        }
        
        
        out_df_RNA <- data.frame(cbind(cor_df$r,cor_df$p,cor_df$n))
        colnames(out_df_RNA) <- c(paste0(input$trait_pc_RNA,".r"), paste0(input$trait_pc_RNA,".p"), paste0(input$trait_pc_RNA,".n"))
        
        cor_df$p <- (-log2(cor_df$p))
        
        if (nrow(cor_df$r)==1 & ncol(cor_df$r)>=1) {
          temp <- data.frame(trait=colnames(cor_df$r),
                             trait_full=names(trait_vector)[match(colnames(cor_df$r), trait_vector)],
                             r=cor_df$r[1,],
                             p=cor_df$p[1,],
                             p_stars=ifelse(cor_df$p[1,] < -log2(0.05), "ns", ifelse(cor_df$p[1,] < -log2(0.01), "*", ifelse(cor_df$p[1,] < -log2(0.001), "**", ifelse(cor_df$p[1,] < -log2(0.0001), "***", "****")))))
          
          temp <- temp[order(temp$r,decreasing = T),] %>% na.omit()
          
          # P>=0.05  NS 
          # 0.01<P<=0.05  *
          # 0.001<P<=0.01  **
          # 0.0001<P<=0.001  ***
          # P<=0.0001
          
          
          out_plot_RNA <- ggbarplot(temp, x = "trait_full", y = "r", fill = "p", color = NA) +
            geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
            geom_text(aes(label = p_stars, y = ifelse(r > 0, r + (range(temp$r)[2]-range(temp$r)[1])*0.01, r - (range(temp$r)[2]-range(temp$r)[1])*0.01)), 
                      vjust = ifelse(temp$r > 0, 0, 1), color = "black") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  #axis.text.x = element_blank(),  # 隐藏 X 轴文本
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                  axis.title.x = element_blank(),  # 隐藏 X 轴标题
                  # axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,colour = "black")
                  axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) +  
            labs(y = paste0(stringr::str_to_title(input$correlation_method_pc), " correlation"), fill = "- Log10 P-value") +
            # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51]) +
            scale_fill_gradientn(colours = get(input$gradient_col)[ifelse(input$gradient_col == "virid",yes = 1, no = 26):51]) 
          
          out_plottype_RNA <- "barplot"
          
          
        } else if (nrow(cor_df$r)>1 & ncol(cor_df$r)==1) {
          temp <- data.frame(gene=rownames(cor_df$r),
                             r=cor_df$r[,1],
                             p=cor_df$p[,1],
                             p_stars=ifelse(cor_df$p[,1] < -log2(0.05), "ns", ifelse(cor_df$p[,1] < -log2(0.01), "*", ifelse(cor_df$p[,1] < -log2(0.001), "**", ifelse(cor_df$p[,1] < -log2(0.0001), "***", "****")))))
          
          temp <- temp[order(temp$r,decreasing = T),]
          
          # out_plot_RNA <- ggscatter(temp, x = "gene", y = "r", color = "p",size = 3)+
          #   theme(legend.position = "right",
          #         axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) +
          #   labs(x = "Genes", y = paste0(stringr::str_to_title(input$correlation_method_pc)," correlation"), color = "- Log10 P-value") +
          #   scale_color_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51])
          
          out_plot_RNA <- ggbarplot(temp, x = "gene", y = "r", fill = "p", color = NA) +
            geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
            geom_text(aes(label = p_stars, y = ifelse(r > 0, r + (range(temp$r)[2]-range(temp$r)[1])*0.01, r - (range(temp$r)[2]-range(temp$r)[1])*0.01)), 
                      vjust = ifelse(temp$r > 0, 0, 1), color = "black") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  #axis.text.x = element_blank(),  # 隐藏 X 轴文本
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                  axis.title.x = element_blank(),  # 隐藏 X 轴标题
                  # axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,colour = "black")
                  axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) + 
            labs(x = "Genes", y = paste0(stringr::str_to_title(input$correlation_method_pc)," correlation"), fill = "- Log10 P-value") +
            # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51]) +
            scale_fill_gradientn(colours = get(input$gradient_col)[ifelse(input$gradient_col == "virid",yes = 1, no = 26):51])
          
          out_plottype_RNA <- "barplot"
          
          
        } else if (nrow(cor_df$r)>1 & nrow(cor_df$r)>1) {
          row_order <- hclust(d = dist(cor_df$r))$order
          column_order <- hclust(d = dist(t(cor_df$r)))$order
          cor_df$r <- cor_df$r[row_order,column_order]
          cor_df$p <- cor_df$p[row_order,column_order]
          cor_df$r_full <- cor_df$r
          colnames(cor_df$r_full) <- names(trait_vector)[match(colnames(cor_df$r_full), trait_vector)]
          
          out_plot_RNA <- pheatmap::pheatmap(mat = cor_df$r_full,
                                             cluster_cols = T,
                                             cluster_rows = T,
                                             display_numbers = ifelse(test = cor_df$p > (-log2(0.001)),
                                                                      yes =  "***",
                                                                      no =  ifelse(test = cor_df$p > (-log2(0.01)),
                                                                                   yes =  "**",
                                                                                   no =  ifelse(test = cor_df$p > (-log2(0.05)),
                                                                                                yes =  "*",
                                                                                                no =  ""))),
                                             cellwidth = 20,
                                             cellheight = 20,
                                             border_color = NA,
                                             treeheight_col = 30,
                                             treeheight_row = 30,
                                             fontsize = 15,
                                             angle_col = 90,
                                             main = "Heatmap of correlation",
                                             silent = T,
                                             show_colnames = T,
                                             fontfamily= "Red Hat Display",
                                             clustering_method = input$cluster_method,
                                             # color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51),
                                             color = get(input$gradient_col)
          )
          
          out_plot_RNA <- as.ggplot(out_plot_RNA) + theme(text = element_text(family = "Red Hat Display"))
          
          out_plottype_RNA <- "heatmap"
          
        }
        
        if (input$cohort_pc_RNA %in% "META_sc") {
          out_cohort_RNA <- unique(meta_sc_combination$cohort[meta_sc_combination$gene %in% rownames(cor_df$r) & meta_sc_combination$trait %in% colnames(cor_df$r)])
          out_cohort_RNA <- cohort_name[cohort_name %in% out_cohort_RNA]
        } else if (input$cohort_pc_RNA %in% "META_om") {
          out_cohort_RNA <- unique(meta_om_combination$cohort[meta_om_combination$gene %in% rownames(cor_df$r) & meta_om_combination$trait %in% colnames(cor_df$r)])
          out_cohort_RNA <- cohort_name[cohort_name %in% out_cohort_RNA]
        } else {
          out_cohort_RNA <- cohort_name[cohort_name %in% input$cohort_pc_RNA]
        }
        
        
        
      } else if (length(input$cohort_pc_RNA) > 1 && (!all(input$cohort_pc_RNA %in% c("META_sc", "META_om")))) {
        
        if (length(input$adjustment_pc_RNA) > 0) {
          
          cor <- list(r=c(),p=c(),n=c(),ID=c())
          
          for (cohort in input$cohort_pc_RNA) {
            eset <- get(cohort)
            
            pData(eset) <- pData(eset)[,colSums(is.na(pData(eset))) < nrow(pData(eset))]
            
            if ((gene %in% rownames(eset)) & (input$trait_pc_RNA %in% colnames(pData(eset)))) {
              
              temp <- t(rbind(Biobase::exprs(eset)[gene,,drop=F],t(pData(eset)[,colnames(pData(eset)) %in% c(input$trait_pc_RNA,input$adjustment_pc_RNA),drop=F])))  %>% na.omit()
              
              # 如果存在Gender列，将其转换为数值型
              if ("Gender" %in% colnames(temp)) {
                temp[, "Gender"] <- ifelse(temp[, "Gender"] == "m", 1, 0)
              }
              
              # 保存行名
              rownames_temp <- rownames(temp)
              
              # 将temp矩阵的每一列转换为数值型
              temp <- apply(temp, 2, as.numeric)
              
              # 重新设置行名
              rownames(temp) <- rownames_temp
              
              if (any(input$adjustment_pc_RNA %in% colnames(pData(eset)))) {
                cor_df <- pcor.test(x = temp[, gene], y = temp[, input$trait_pc_RNA], z = temp[, colnames(temp) %in% input$adjustment_pc_RNA], method = input$correlation_method_pc)
                cor$r <- c(cor$r, cor_df$estimate)
                cor$p <- c(cor$p, cor_df$p.value)
                cor$n <- c(cor$n, cor_df$n)
                cor$ID <- c(cor$ID, names(cohort_name[cohort_name == cohort]))
                
              } else {
                cor$r <- c(cor$r, get(paste0('cor_list_', cohort, '_', input$correlation_method_pc))$r[gene,input$trait_pc_RNA])
                cor$p <- c(cor$p, get(paste0('cor_list_', cohort, '_', input$correlation_method_pc))$p[gene,input$trait_pc_RNA])
                cor$n <- c(cor$n, get(paste0('cor_list_', cohort, '_', input$correlation_method_pc))$n[gene,input$trait_pc_RNA])
                cor$ID <- c(cor$ID, names(cohort_name[cohort_name == cohort]))
              }
            }
          }
          
        } else {
          
          results <- lapply(input$cohort_pc_RNA, function(cohort) {
            data_name <- paste0('cor_list_', cohort, '_', input$correlation_method_pc)
            data <- get(data_name)
            
            # Check if the gene exists in the cohort
            if (gene[1] %in% rownames(data$r) & input$trait_pc_RNA %in% colnames(data$r)) {
              return(list(
                r = data$r[gene[1], input$trait_pc_RNA],
                p = data$p[gene[1], input$trait_pc_RNA],
                n = data$n[gene[1], input$trait_pc_RNA],
                ID = names(cohort_name[cohort_name == cohort])
              ))
            } else {
              return(NULL)  # Return NULL for cohorts that don't have the gene
            }
          })
          
          # Filter out NULL results
          results <- Filter(Negate(is.null), results)
          
          cor <- list(
            r = unlist(lapply(results, `[[`, "r")),
            p = unlist(lapply(results, `[[`, "p")),
            n = unlist(lapply(results, `[[`, "n")),
            ID = unlist(lapply(results, `[[`, "ID"))
          )
          
        }
        
        out_df_RNA <- data.frame(r=cor$r, p=cor$p, n=cor$n, cohort=cor$ID)
        
        out_plot_RNA <- meta::metacor(cor = cor$r, n = cor$n, studlab = cor$ID)
        
        out_plottype_RNA <- "metacor"
        
        out_cohort_RNA <- cohort_name[cor$ID]
        
      }
      
      out$has_plot_RNA <- TRUE
      out$plot_RNA <- out_plot_RNA
      out$plottype_RNA <- out_plottype_RNA
      out$df_RNA <- out_df_RNA
      out$cohort_RNA <- out_cohort_RNA
    }
    
    
    
    
    if (input$proteomics_pc==T) {
      
      gene_protein <- intersect(unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_pc),split='\n|\t|,')) %>% unique(), rownames(Proteomicseset))
      
      trait_protein <- intersect(input$trait_pc_RNA, candidate_traits_pc_protein())
      
      adjustment_protein <- intersect(input$adjustment_pc_RNA, candidate_traits_pc_protein())
      
      
      if (length(gene_protein) > 0 & length(trait_protein) > 0) {
        
        if (length(adjustment_protein) > 0) {
          
          cor_df <- list()
          cor_df$r <- matrix(NA, nrow = length(gene_protein), ncol = length(trait_protein), dimnames = list(gene_protein, trait_protein))
          cor_df$p <- cor_df$r
          cor_df$n <- cor_df$r
          
          # 使用expand.grid生成所有基因和特征的组合
          combinations <- expand.grid(gene = gene_protein, trait = trait_protein)
          
          # 对每一对组合应用函数
          results <- apply(combinations, 1, function(x) {
            genetemp <- x['gene']
            trait <- x['trait']
            temp <- t(rbind(Biobase::exprs(Proteomicseset)[genetemp,,drop=F],t(pData(Proteomicseset)[,c(trait,adjustment_protein)]))) %>% na.omit()
            
            # 如果存在Gender列，将其转换为数值型
            if ("Gender" %in% colnames(temp)) {
              temp[, "Gender"] <- ifelse(temp[, "Gender"] == "m", 1, 0)
            }
            
            # 保存行名
            rownames_temp <- rownames(temp)
            
            # 将temp矩阵的每一列转换为数值型
            temp <- apply(temp, 2, as.numeric)
            
            # 重新设置行名
            rownames(temp) <- rownames_temp
            
            cor_df_temp_temp <- pcor.test(x = temp[, genetemp], y = temp[, trait], z = temp[, adjustment_protein], method = input$correlation_method_pc)
            c(r = cor_df_temp_temp$estimate, p = cor_df_temp_temp$p.value, n = cor_df_temp_temp$n)
          })
          
          cor_df$r[] <- results['r',]
          cor_df$p[] <- results['p',]
          cor_df$n[] <- results['n',]
          
        } else {
          cor_df <- list()
          cor_df$r <- get(paste0('cor_list_Proteomicseset_', input$correlation_method_pc))$r[gene_protein,trait_protein,drop=F]
          cor_df$p <- get(paste0('cor_list_Proteomicseset_', input$correlation_method_pc))$p[gene_protein,trait_protein,drop=F]
          cor_df$n <- get(paste0('cor_list_Proteomicseset_', input$correlation_method_pc))$n[gene_protein,trait_protein,drop=F]
        }
        
        out_df_protein <- data.frame(cbind(cor_df$r,cor_df$p,cor_df$n))
        colnames(out_df_protein) <- c(paste0(trait_protein,".r"), paste0(trait_protein,".p"), paste0(trait_protein,".n"))
        
        cor_df$p <- (-log2(cor_df$p))
        
        if (nrow(cor_df$r)==1 & ncol(cor_df$r)>=1) {
          temp <- data.frame(trait=colnames(cor_df$r),
                             trait_full=names(trait_vector)[match(colnames(cor_df$r), trait_vector)],
                             r=cor_df$r[1,],
                             p=cor_df$p[1,],
                             p_stars=ifelse(cor_df$p[1,] < -log2(0.05), "ns", ifelse(cor_df$p[1,] < -log2(0.01), "*", ifelse(cor_df$p[1,] < -log2(0.001), "**", ifelse(cor_df$p[1,] < -log2(0.0001), "***", "****")))))
          
          temp <- temp[order(temp$r,decreasing = T),] %>% na.omit()
          
          # P>=0.05  NS 
          # 0.01<P<=0.05  *
          # 0.001<P<=0.01  **
          # 0.0001<P<=0.001  ***
          # P<=0.0001
          
          
          out_plot_protein <- ggbarplot(temp, x = "trait_full", y = "r", fill = "p", color = NA) +
            geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
            geom_text(aes(label = p_stars, y = ifelse(r > 0, r + (range(temp$r)[2]-range(temp$r)[1])*0.01, r - (range(temp$r)[2]-range(temp$r)[1])*0.01)), 
                      vjust = ifelse(temp$r > 0, 0, 1), color = "black") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  #axis.text.x = element_blank(),  # 隐藏 X 轴文本
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                  axis.title.x = element_blank(),  # 隐藏 X 轴标题
                  # axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,colour = "black")
                  axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) +  
            labs(y = paste0(stringr::str_to_title(input$correlation_method_pc), " correlation"), fill = "- Log10 P-value") +
            # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51]) +
            scale_fill_gradientn(colours = get(input$gradient_col)[ifelse(input$gradient_col == "virid",yes = 1, no = 26):51]) 
          
          out_plottype_protein <- "barplot"
          
          
        } else if (nrow(cor_df$r)>1 & ncol(cor_df$r)==1) {
          temp <- data.frame(gene=rownames(cor_df$r),
                             r=cor_df$r[,1],
                             p=cor_df$p[,1],
                             p_stars=ifelse(cor_df$p[,1] < -log2(0.05), "ns", ifelse(cor_df$p[,1] < -log2(0.01), "*", ifelse(cor_df$p[,1] < -log2(0.001), "**", ifelse(cor_df$p[,1] < -log2(0.0001), "***", "****")))))
          
          temp <- temp[order(temp$r,decreasing = T),]
          
          # out_plot_protein <- ggscatter(temp, x = "gene", y = "r", color = "p",size = 3)+
          #   theme(legend.position = "right",
          #         axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) +
          #   labs(x = "Genes", y = paste0(stringr::str_to_title(input$correlation_method_pc)," correlation"), color = "- Log10 P-value") +
          #   scale_color_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51])
          
          out_plot_protein <- ggbarplot(temp, x = "gene", y = "r", fill = "p", color = NA) +
            geom_hline(yintercept = 0, color = "black") +  # 添加 y=0 的线
            geom_text(aes(label = p_stars, y = ifelse(r > 0, r + (range(temp$r)[2]-range(temp$r)[1])*0.01, r - (range(temp$r)[2]-range(temp$r)[1])*0.01)), 
                      vjust = ifelse(temp$r > 0, 0, 1), color = "black") +
            theme(legend.position = "right",
                  text = element_text(family = "Red Hat Display"),
                  #axis.text.x = element_blank(),  # 隐藏 X 轴文本
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank(),  # 隐藏 X 轴刻度
                  axis.title.x = element_blank(),  # 隐藏 X 轴标题
                  # axis.text.x = element_text(angle=45,hjust = 1,vjust = 1,colour = "black")
                  axis.text.x = element_text(angle=90,hjust = 1,vjust = 0.5,colour = "black")) + 
            labs(x = "Genes", y = paste0(stringr::str_to_title(input$correlation_method_pc)," correlation"), fill = "- Log10 P-value") +
            # scale_fill_gradientn(colours = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)[35:51]) +
            scale_fill_gradientn(colours = get(input$gradient_col)[ifelse(input$gradient_col == "virid",yes = 1, no = 26):51])
          
          out_plottype_protein <- "barplot"
          
          
        } else if (nrow(cor_df$r)>1 & nrow(cor_df$r)>1) {
          row_order <- hclust(d = dist(cor_df$r))$order
          column_order <- hclust(d = dist(t(cor_df$r)))$order
          cor_df$r <- cor_df$r[row_order,column_order]
          cor_df$p <- cor_df$p[row_order,column_order]
          cor_df$r_full <- cor_df$r
          colnames(cor_df$r_full) <- names(trait_vector)[match(colnames(cor_df$r_full), trait_vector)]
          
          out_plot_protein <- pheatmap::pheatmap(mat = cor_df$r_full,
                                                 cluster_cols = T,
                                                 cluster_rows = T,
                                                 display_numbers = ifelse(test = cor_df$p > (-log2(0.001)),
                                                                          yes =  "***",
                                                                          no =  ifelse(test = cor_df$p > (-log2(0.01)),
                                                                                       yes =  "**",
                                                                                       no =  ifelse(test = cor_df$p > (-log2(0.05)),
                                                                                                    yes =  "*",
                                                                                                    no =  ""))),
                                                 cellwidth = 20,
                                                 cellheight = 20,
                                                 border_color = NA,
                                                 treeheight_col = 30,
                                                 treeheight_row = 30,
                                                 fontsize = 15,
                                                 angle_col = 90,
                                                 main = "Heatmap of correlation",
                                                 silent = T,
                                                 show_colnames = T,
                                                 fontfamily= "Red Hat Display",
                                                 clustering_method = input$cluster_method,
                                                 # color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51),
                                                 color = get(input$gradient_col)
          )
          
          out_plot_protein <- as.ggplot(out_plot_protein) + theme(text = element_text(family = "Red Hat Display"))
          
          out_plottype_protein <- "heatmap"
          
        }
        
        out$has_plot_protein <- TRUE
        out$plot_protein <- out_plot_protein
        out$plottype_protein <- out_plottype_protein
        out$df_protein <- out_df_protein
        
      }
    }
    
    return(out)
    
  }) # eventReactive
  
  
  
  
  ## Sex ----
  
  output$ui_statistics_sex <- renderUI({
    # 如果只选择了一个cohort
    if(length(input$cohort_sex) == 1) {
      # 返回PickerInput
      return(list(prettyRadioButtons(inputId = "statistics_sex",
                                     label = "Statistical test",
                                     selected = "wilcox.test",
                                     choices = c("t-Test" = "t.test", "Wilcoxon Test" = "wilcox.test"),
                                     inline = TRUE))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_group_sex消失
      return(NULL)
    }
  })
  
  
  output$ui_effect_sex <- renderUI({
    # 如果选择了多个cohort
    if(length(input$cohort_sex) > 1) {
      # 返回PickerInput
      return(list(prettyRadioButtons(inputId = "effect_sex",
                                     label = "Effect measure (Transcriptomics, optional, available for one gene)",
                                     selected = "standardized mean difference",
                                     choices = c("log fold change" = "log fold change", "standardized mean difference" = "standardized mean difference"),
                                     inline = TRUE))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_group_sex消失
      return(NULL)
    }
  })
  
  
  output$ui_statistics_sex_protein <- renderUI({
    # 如果只选择了一个cohort
    if(length(input$cohort_sex) > 1 & input$proteomics_sex==T) {
      # 返回PickerInput
      return(list(prettyRadioButtons(inputId = "statistics_sex_protein",
                                     label = "Statistical test (Proteomics)",
                                     selected = "wilcox.test",
                                     choices = c("t-Test" = "t.test", "Wilcoxon Test" = "wilcox.test"),
                                     inline = TRUE))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_group_sex消失
      return(NULL)
    }
  })
  
  
  plotdata_sex <- eventReactive(input$SearchButton_sex, {
    
    gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_sex),split='\n|\t|,')) %>% unique()
    
    # 使用lapply遍历所有的cohort并获取每个cohort中存在的基因  Use lapply to loop through all cohorts and get the genes present in each cohort
    all_genes_in_cohorts <- lapply(input$cohort_sex, function(cohort) {
      rownames(Biobase::exprs(get(cohort)))
    }) %>% unlist() %>% unique()
    
    gene <- intersect(gene, all_genes_in_cohorts)
    
    out <- list(has_plot_RNA = FALSE, plot_RNA = NULL, plottype_RNA = NULL, df_RNA = NULL, cohort_RNA = NULL, has_plot_protein = FALSE, plot_protein = NULL, plottype_protein = NULL, df_protein = NULL) # 初始化返回列表
    
    if (length(gene) > 0) {
      if (length(input$cohort_sex)==1) {
        eset <- get(input$cohort_sex)
        
        plot_list <- lapply(gene, function(i) {
          df <- rbind(
            data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$Gender == "m"]], sex = "Male"),
            data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$Gender == "f"]], sex = "Female")
          )
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2[1:2]
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
          }
          ggboxplot(data = df, x = "sex", y = "expression", fill = "sex", xlab = "", ylab = paste0(i, " expression")) +
            theme(aspect.ratio = 1, 
                  axis.title.x = element_blank(), 
                  text = element_text(family = "Red Hat Display")) + 
            NoLegend() +
            stat_compare_means(comparisons = list(c("Male", "Female")), method = input$statistics_sex) +
            scale_fill_manual(values = discrete_color)
        })
        
        out_plot_RNA <- plot_grid(plotlist = plot_list,align = "hv",nrow = ceiling(length(gene)/5))
        
        out_plottype_RNA <- "boxplot"
        
        out$row_RNA <- ceiling(length(gene)/5)
        
        # Create a long-format data frame
        out_df_RNA <- data.frame(
          gene = rep(rownames(Biobase::exprs(eset)[gene, , drop = FALSE]), ncol(Biobase::exprs(eset)[gene, , drop = FALSE])),
          expression = as.vector(Biobase::exprs(eset)[gene, , drop = FALSE]),
          sex = rep(ifelse(eset$Gender == "m", "Male", "Female"), each = nrow(Biobase::exprs(eset)[gene, , drop = FALSE]))
        )
        
        # Calculate the mean expression for each gene by sex
        out_df_RNA <- out_df_RNA %>%
          group_by(gene, sex) %>%
          summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = sex, values_from = mean_expression) %>% data.frame()
        
        # Set the rownames and remove the gene column
        rownames(out_df_RNA) <- out_df_RNA$gene
        out_df_RNA$gene <- NULL
        
        colnames(out_df_RNA) <- paste0("Mean_",colnames(out_df_RNA))
        
        out_cohort_RNA <- cohort_name[cohort_name %in% input$cohort_sex]
        
      } else if (length(gene) > 1 & length(input$cohort_sex) > 1) {
        df <- matrix(data = NA, nrow = length(gene), ncol = length(input$cohort_sex))
        rownames(df) <- gene
        colnames(df) <- input$cohort_sex
        
        for (j in input$cohort_sex) {
          eset <- get(j)
          for (i in gene) {
            if (i %in% rownames(Biobase::exprs(eset))) {
              df[i,j] <- log2(mean(Biobase::exprs(eset)[i,colnames(eset)[eset$Gender=="m"]])/mean(Biobase::exprs(eset)[i,colnames(eset)[eset$Gender=="f"]]))
            }
          }
          colnames(df)[colnames(df)==j] <- names(cohort_name[cohort_name %in% j])
        }
        
        out_plot_RNA <- pheatmap::pheatmap(mat = df,
                                           cluster_cols = ncol(df)>1,
                                           cluster_rows = ncol(df)>1,
                                           cellwidth = 20,
                                           cellheight = 20,
                                           border_color = NA,
                                           treeheight_col = 30,
                                           treeheight_row = 30,
                                           fontsize = 15,
                                           main = "Heatmap of log fold change",
                                           silent = T,
                                           show_colnames = T,
                                           fontfamily= "Red Hat Display",
                                           clustering_method = input$cluster_method,
                                           # color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
                                           color = get(input$gradient_col)
        )
        
        out_plot_RNA <- as.ggplot(out_plot_RNA) + theme(text = element_text(family = "Red Hat Display"))
        
        out_plottype_RNA <- "heatmap"
        
        out_df_RNA <- as.data.frame(df)
        
        out_cohort_RNA <- cohort_name[colnames(df)[colSums(is.na(df))<nrow(df)]]
        
      } else if (length(gene) == 1 & length(input$cohort_sex) > 1) {
        
        
        male_list <- list()
        female_list <- list()
        
        for (i in input$cohort_sex) {
          eset <- get(i)
          
          if (gene %in% rownames(Biobase::exprs(eset))) {
            
            male_list$n <- c(male_list$n, length(Biobase::exprs(eset)[gene,colnames(eset)[eset$Gender=="m"]]))
            male_list$mean <- c(male_list$mean, mean(Biobase::exprs(eset)[gene,colnames(eset)[eset$Gender=="m"]]))
            male_list$sd <- c(male_list$sd, sd(Biobase::exprs(eset)[gene,colnames(eset)[eset$Gender=="m"]]))
            male_list$ID <- c(male_list$ID, names(cohort_name[cohort_name %in% i]))
            
            female_list$n <- c(female_list$n, length(Biobase::exprs(eset)[gene,colnames(eset)[eset$Gender=="f"]]))
            female_list$mean <- c(female_list$mean, mean(Biobase::exprs(eset)[gene,colnames(eset)[eset$Gender=="f"]]))
            female_list$sd <- c(female_list$sd, sd(Biobase::exprs(eset)[gene,colnames(eset)[eset$Gender=="f"]]))
            female_list$ID <- c(female_list$ID, names(cohort_name[cohort_name %in% i]))
            
          }
          
        }
        
        
        out_plottype_RNA <- "metacont"
        
        out_df_RNA <- data.frame(Male_n=male_list$n,
                                 Male_mean=male_list$mean,
                                 Male_sd=male_list$sd,
                                 Female_n=female_list$n,
                                 Female_mean=female_list$mean,
                                 Female_sd=female_list$sd,
                                 Trait=female_list$ID
        )
        
        out$gene_RNA <- gene
        
        out$missing_gene_cohort_RNA <- names(cohort_name[cohort_name %in% setdiff(input$cohort_sex,cohort_name[male_list$ID])])
        
        out_cohort_RNA <- cohort_name[male_list$ID]
        
        
        if (input$effect_sex=="log fold change") {
          
          out_plot_RNA <- meta::metacont(n.e = male_list$n,
                                         mean.e = male_list$mean,
                                         sd.e = male_list$sd,
                                         n.c = female_list$n,
                                         mean.c = female_list$mean,
                                         sd.c = female_list$sd,
                                         studlab = male_list$ID,
                                         label.e = "Male",
                                         label.c = "Female",
                                         label.right = "Male",
                                         label.left = "Female",
                                         sm = "ROM",
                                         backtransf=FALSE)
          
          out_plot_RNA$sm <- "log2 FC"
          out_plot_RNA$TE <- log2(exp(out_plot_RNA$TE))
          out_plot_RNA$TE.common <- log2(exp(out_plot_RNA$TE.common))
          out_plot_RNA$TE.random <- log2(exp(out_plot_RNA$TE.random))
          out_plot_RNA$TE.fixed <- log2(exp(out_plot_RNA$TE.fixed))
          out_plot_RNA$lower <- log2(exp(out_plot_RNA$lower))
          out_plot_RNA$lower.common <- log2(exp(out_plot_RNA$lower.common))
          out_plot_RNA$lower.random <- log2(exp(out_plot_RNA$lower.random))
          out_plot_RNA$upper <- log2(exp(out_plot_RNA$upper))
          out_plot_RNA$upper.common <- log2(exp(out_plot_RNA$upper.common))
          out_plot_RNA$upper.random <- log2(exp(out_plot_RNA$upper.random))
          
          
        } else if (input$effect_sex=="standardized mean difference") {
          
          out_plot_RNA <- meta::metacont(n.e = male_list$n,
                                         mean.e = male_list$mean,
                                         sd.e = male_list$sd,
                                         n.c = female_list$n,
                                         mean.c = female_list$mean,
                                         sd.c = female_list$sd,
                                         studlab = male_list$ID,
                                         label.e = "Male",
                                         label.c = "Female",
                                         label.right = "Male",
                                         label.left = "Female",
                                         sm = "SMD",
                                         backtransf=FALSE,
                                         method.smd = "Hedges")
          
          
        }
        
      }
      
      out$has_plot_RNA <- TRUE
      out$plot_RNA <- out_plot_RNA
      out$plottype_RNA <- out_plottype_RNA
      out$df_RNA <- out_df_RNA
      out$cohort_RNA <- out_cohort_RNA
    }
    
    
    
    
    if (input$proteomics_sex==T) {
      
      gene_protein <- intersect(unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_sex),split='\n|\t|,')) %>% unique(), rownames(Proteomicseset))
      
      if (length(gene_protein) > 0) {
        
        eset <- Proteomicseset
        
        plot_list <- lapply(gene_protein, function(i) {
          df <- rbind(
            data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$Gender == "m"]], sex = "Male"),
            data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$Gender == "f"]], sex = "Female")
          ) %>% na.omit()
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2[1:2]
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
          }
          ggboxplot(data = df, x = "sex", y = "expression", fill = "sex", xlab = "", ylab = paste0(i, " expression")) +
            theme(aspect.ratio = 1, 
                  axis.title.x = element_blank(), 
                  text = element_text(family = "Red Hat Display"),
                  legend.position = "none") + 
            stat_compare_means(comparisons = list(c("Male", "Female")), method = ifelse(test = length(input$cohort_sex) == 1,yes = input$statistics_sex,no = input$statistics_sex_protein)) +
            scale_fill_manual(values = discrete_color)
        })
        
        out_plot_protein <- plot_grid(plotlist = plot_list,align = "hv",nrow = ceiling(length(gene_protein)/5))
        
        out_plottype_protein <- "boxplot"
        
        out$row_protein <- ceiling(length(gene_protein)/5)
        
        # Create a long-format data frame
        out_df_protein <- data.frame(
          gene = rep(rownames(Biobase::exprs(eset)[gene_protein, , drop = FALSE]), ncol(Biobase::exprs(eset)[gene_protein, , drop = FALSE])),
          expression = as.vector(Biobase::exprs(eset)[gene_protein, , drop = FALSE]),
          sex = rep(ifelse(eset$Gender == "m", "Male", "Female"), each = nrow(Biobase::exprs(eset)[gene_protein, , drop = FALSE]))
        )
        
        # Calculate the mean expression for each gene by sex
        out_df_protein <- out_df_protein %>%
          group_by(gene, sex) %>%
          summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = sex, values_from = mean_expression) %>% data.frame()
        
        # Set the rownames and remove the gene column
        rownames(out_df_protein) <- out_df_protein$gene
        out_df_protein$gene <- NULL
        
        colnames(out_df_protein) <- paste0("Mean_",colnames(out_df_protein))
        
        out$has_plot_protein <- TRUE
        out$plot_protein <- out_plot_protein
        out$plottype_protein <- out_plottype_protein
        out$df_protein <- out_df_protein
        
      }
      
    }
    
    return(out)
    
  }) # eventReactive
  
  
  
  
  ## BMI cross-sectional----
  
  output$ui_group_bmicross_RNA <- renderUI({
    
    # 如果只选择了一个cohort
    if(length(input$cohort_bmicross) == 1) {
      
      # 获取该cohort的BMI分组
      bmi_groups <- cohort_BMI_list[[input$cohort_bmicross]]
      
      if (c("Underweight","Healthy Weight", "Overweight") %in% bmi_groups %>% any()) {
        bmi_groups <- c(bmi_groups,"Non Obese")
      }
      
      if (c("Obese I","Obese II","Obese III") %in% bmi_groups %>% any()) {
        bmi_groups <- c(bmi_groups,"Obese")
      }
      
      name <- bmi_groups
      name[name=="Underweight"] <- "BMI < 18.5"
      name[name=="Healthy Weight"] <- "BMI: 18.5 ~ 25"
      name[name=="Overweight"] <- "BMI: 25 ~ 30"
      name[name=="Obese I"] <- "BMI: 30 ~ 35"
      name[name=="Obese II"] <- "BMI: 35 ~ 40"
      name[name=="Obese III"] <- "BMI > 40"
      name[name=="Non Obese"] <- "BMI < 30"
      name[name=="Obese"] <- "BMI > 30"
      # names(bmi_groups) <- name
      
      
      # 返回PickerInput
      
      return(list(pickerInput(inputId = "group_bmicross_RNA",
                              label = "BMI group",
                              choices = bmi_groups, 
                              selected = if (all(c("Obese","Non Obese") %in% bmi_groups)) {c("Obese","Non Obese")} else {NULL},
                              multiple = TRUE,
                              # width = "300px",
                              options = pickerOptions(actionsBox = T, size = 10, liveSearch = T),
                              choicesOpt = list(subtext = name))
      )
      )
      
      
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_group_bmi消失
      return(NULL)
    }
  })
  
  
  output$ui_group_bmicross_protein <- renderUI({
    
    # 如果选择了proteomics
    if(input$proteomics_bmicross == T & length(input$cohort_bmicross) > 1) {
      
      # 获取该cohort的BMI分组
      bmi_groups <- c("Healthy Weight", "Overweight", "Obese I", "Obese II", "Obese III", "Non Obese", "Obese")
      
      name <- bmi_groups
      name[name=="Underweight"] <- "BMI < 18.5"
      name[name=="Healthy Weight"] <- "BMI: 18.5 ~ 25"
      name[name=="Overweight"] <- "BMI: 25 ~ 30"
      name[name=="Obese I"] <- "BMI: 30 ~ 35"
      name[name=="Obese II"] <- "BMI: 35 ~ 40"
      name[name=="Obese III"] <- "BMI > 40"
      name[name=="Non Obese"] <- "BMI < 30"
      name[name=="Obese"] <- "BMI > 30"
      # names(bmi_groups) <- name
      
      
      # 返回PickerInput
      
      return(list(pickerInput(inputId = "group_bmicross_protein",
                              label = "BMI group (Proteomics)",
                              choices = bmi_groups,
                              selected = if (all(c("Obese","Non Obese") %in% bmi_groups)) {c("Obese","Non Obese")} else {NULL},
                              multiple = TRUE,
                              # width = "300px",
                              options = pickerOptions(actionsBox = T, size = 10, liveSearch = T),
                              choicesOpt = list(subtext = name))
      )
      )
      
      
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_group_bmi消失
      return(NULL)
    }
  })
  
  
  output$ui_statistics_bmicross <- renderUI({
    # 如果只选择了一个cohort
    if(length(input$cohort_bmicross) == 1) {
      # 返回PickerInput
      
      return(list(prettyRadioButtons(inputId = "statistics_bmicross",
                                     label = "Statistical test",
                                     selected = "wilcox.test",
                                     choices = c("t-Test" = "t.test", "Wilcoxon Test" = "wilcox.test"),
                                     inline = TRUE)
      )
      )
      
      
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_group_bmi消失
      return(NULL)
    }
  })
  
  
  output$ui_statistics_bmicross_protein <- renderUI({
    # 如果只选择了一个cohort
    if(length(input$cohort_bmicross) > 1 & input$proteomics_bmicross==T) {
      # 返回PickerInput
      return(list(prettyRadioButtons(inputId = "statistics_bmicross_protein",
                                     label = "Statistical test (Proteomics)",
                                     selected = "wilcox.test",
                                     choices = c("t-Test" = "t.test", "Wilcoxon Test" = "wilcox.test"),
                                     inline = TRUE))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_statistics_bmicross_protein消失
      return(NULL)
    }
  })
  
  
  output$ui_effect_bmicross <- renderUI({
    # 如果选择了多个cohort
    if(length(input$cohort_bmicross) > 1) {
      # 返回PickerInput
      return(list(prettyRadioButtons(inputId = "effect_bmicross",
                                     label = "Effect measure (Transcriptomics, optional, available for one gene)",
                                     selected = "standardized mean difference",
                                     choices = c("log fold change" = "log fold change", "standardized mean difference" = "standardized mean difference"),
                                     inline = TRUE))
             
      )
    } else {
      # 如果选择了多个cohort或没有选择，返回NULL使ui_effect_bmicross消失
      return(NULL)
    }
  })
  
  
  plotdata_bmicross <- eventReactive(input$SearchButton_bmicross, {
    
    gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_bmicross),split='\n|\t|,')) %>% unique()
    
    # 使用lapply遍历所有的cohort并获取每个cohort中存在的基因  Use lapply to loop through all cohorts and get the genes present in each cohort
    all_genes_in_cohorts <- lapply(input$cohort_bmicross, function(cohort) {
      rownames(Biobase::exprs(get(cohort)))
    }) %>% unlist() %>% unique()
    
    gene <- intersect(gene, all_genes_in_cohorts)
    
    out <- list(has_plot_RNA = FALSE, plot_RNA = NULL, plottype_RNA = NULL, df_RNA = NULL, cohort_RNA = NULL, has_plot_protein = FALSE, plot_protein = NULL, plottype_protein = NULL, df_protein = NULL) # 初始化返回列表
    
    if (length(gene) > 0) {
      if (length(input$cohort_bmicross)==1) {
        eset <- get(input$cohort_bmicross)
        eset <- eset[,!is.na(eset$BMI_catelogy)]
        
        compare_list <- Filter(function(combo) {
          !(combo[1] == "Obese" & combo[2] %in% c("Obese I", "Obese II", "Obese III")) &
            !(combo[2] == "Obese" & combo[1] %in% c("Obese I", "Obese II", "Obese III")) &
            !(combo[1] == "Non Obese" & combo[2] %in% c("Underweight", "Healthy Weight", "Overweight")) &
            !(combo[2] == "Non Obese" & combo[1] %in% c("Underweight", "Healthy Weight", "Overweight"))
        }, combn(input$group_bmicross_RNA, 2, simplify = FALSE))
        
        
        plot_list <- lapply(gene, function(i) {
          
          # 根据input$group_bmi的选择来确定哪些BMI类别应该被包括
          df_list <- lapply(input$group_bmicross_RNA, function(bmi_choice) {
            switch(bmi_choice,
                   "Underweight" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Underweight"]], BMI_catelogy = "Underweight"),
                   "Healthy Weight" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Healthy Weight"]], BMI_catelogy = "Healthy Weight"),
                   "Overweight" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Overweight"]], BMI_catelogy = "Overweight"),
                   "Obese I" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Obese I"]], BMI_catelogy = "Obese I"),
                   "Obese II" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Obese II"]], BMI_catelogy = "Obese II"),
                   "Obese III" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Obese III"]], BMI_catelogy = "Obese III"),
                   "Non Obese" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")]], BMI_catelogy = "Non Obese"),
                   "Obese" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]], BMI_catelogy = "Obese")
            )
          })
          
          # 合并所有的数据框
          df <- do.call(rbind, df_list)
          
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(length(unique(df$BMI_catelogy)))
          }
          # 为每个基因生成盒形图
          p <- ggboxplot(data = df, x = "BMI_catelogy", y = "expression", fill = "BMI_catelogy", xlab = "", ylab = paste0(i, " expression")) +
            theme(aspect.ratio = 1, 
                  axis.title.x = element_blank(), 
                  axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
                  text = element_text(family = "Red Hat Display")) +
            NoLegend() +
            stat_compare_means(comparisons = compare_list, method = input$statistics_bmicross) +
            scale_fill_manual(values = discrete_color)
          
          return(p)
        })
        
        
        # Create a function to get the average expression for a given gene and BMI category
        get_avg_expression <- function(gene, bmi_category) {
          switch(bmi_category,
                 "Underweight" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Underweight"], na.rm = TRUE),
                 "Healthy Weight" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Healthy Weight"], na.rm = TRUE),
                 "Overweight" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Overweight"], na.rm = TRUE),
                 "Obese I" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Obese I"], na.rm = TRUE),
                 "Obese II" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Obese II"], na.rm = TRUE),
                 "Obese III" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Obese III"], na.rm = TRUE),
                 "Non Obese" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")], na.rm = TRUE),
                 "Obese" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")], na.rm = TRUE)
          )
        }
        
        # Create the output data frame
        out_df_RNA <- as.data.frame(matrix(0, nrow = length(gene), ncol = length(input$group_bmicross_RNA)))
        rownames(out_df_RNA) <- gene
        colnames(out_df_RNA) <- input$group_bmicross_RNA
        
        # Fill the data frame with average expression values
        for (g in gene) {
          for (bmi_cat in input$group_bmicross_RNA) {
            out_df_RNA[g, bmi_cat] <- get_avg_expression(g, bmi_cat)
          }
        }
        
        out_plot_RNA <- plot_grid(plotlist = plot_list,align = "hv", nrow = ceiling(length(gene)/5))
        
        out_plottype_RNA <- "boxplot"
        
        out$row_RNA <- ceiling(length(gene)/5)
        
        out_cohort_RNA <- cohort_name[cohort_name %in% input$cohort_bmicross]
        
      } else if (length(gene) > 1 & length(input$cohort_bmicross) > 1) {
        
        # 去除没用Obese和Non obese的cohort
        cohort_bmicross <- Filter(function(i) {
          bmi_categories <- get(i)$BMI_catelogy
          any(bmi_categories %in% c("Obese I", "Obese II", "Obese III")) && 
            any(bmi_categories %in% c("Underweight", "Healthy Weight", "Overweight"))
        }, input$cohort_bmicross)
        
        df <- matrix(data = NA,nrow = length(gene),ncol = length(cohort_bmicross))
        rownames(df) <- gene
        colnames(df) <- cohort_bmicross
        
        for (j in cohort_bmicross) {
          eset <- get(j)
          for (i in gene) {
            if (i %in% rownames(Biobase::exprs(eset))) {
              df[i,j] <- log2(mean(Biobase::exprs(eset)[i,colnames(eset)[eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]])/mean(Biobase::exprs(eset)[i,colnames(eset)[eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")]]))
            }
          }
          colnames(df)[colnames(df)==j] <- names(cohort_name[cohort_name %in% j])
        }
        
        out_plot_RNA <- pheatmap::pheatmap(mat = df,
                                           cluster_cols = ncol(df)>1,
                                           cluster_rows = ncol(df)>1,
                                           cellwidth = 20,
                                           cellheight = 20,
                                           border_color = NA,
                                           treeheight_col = 30,
                                           treeheight_row = 30,
                                           fontsize = 15,
                                           main = "Heatmap of log fold change",
                                           silent = T,
                                           show_colnames = T,
                                           fontfamily= "Red Hat Display",
                                           clustering_method = input$cluster_method,
                                           # color = colorRampPalette(brewer.pal(n = 7, name ="PRGn"))(51)
                                           color = get(input$gradient_col)
        )
        
        out_plottype_RNA <- "heatmap"
        
        out_plot_RNA <- as.ggplot(out_plot_RNA) + theme(text = element_text(family = "Red Hat Display"))
        
        out_df_RNA <- as.data.frame(df)
        
        out_cohort_RNA <- cohort_name[colnames(df)[colSums(is.na(df))<nrow(df)]]
        
      } else if (length(gene) == 1 & length(input$cohort_bmicross) > 1) {
        
        cohort_bmicross <- Filter(function(i) {
          bmi_categories <- get(i)$BMI_catelogy
          any(bmi_categories %in% c("Obese I", "Obese II", "Obese III")) && 
            any(bmi_categories %in% c("Underweight", "Healthy Weight", "Overweight"))
        }, input$cohort_bmicross)
        
        
        obese_list <- list()
        nonobese_list <- list()
        # missing_gene_cohort <- c()
        
        for (i in cohort_bmicross) {
          eset <- get(i)
          
          if (gene %in% rownames(Biobase::exprs(eset))) {
            obese_list$n <- c(obese_list$n, length(Biobase::exprs(eset)[gene,colnames(eset)[eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]]))
            obese_list$mean <- c(obese_list$mean, mean(Biobase::exprs(eset)[gene,colnames(eset)[eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]]))
            obese_list$sd <- c(obese_list$sd, sd(Biobase::exprs(eset)[gene,colnames(eset)[eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]]))
            obese_list$ID <- c(obese_list$ID, names(cohort_name[cohort_name %in% i]))
            
            nonobese_list$n <- c(nonobese_list$n, length(Biobase::exprs(eset)[gene,colnames(eset)[eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")]]))
            nonobese_list$mean <- c(nonobese_list$mean, mean(Biobase::exprs(eset)[gene,colnames(eset)[eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")]]))
            nonobese_list$sd <- c(nonobese_list$sd, sd(Biobase::exprs(eset)[gene,colnames(eset)[eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")]]))
            nonobese_list$ID <- c(nonobese_list$ID, names(cohort_name[cohort_name %in% i]))
          }
          
        }
        
        out_plottype_RNA <- "metacont"
        
        out_df_RNA <- data.frame(Obese_n=obese_list$n,
                                 Obese_mean=obese_list$mean,
                                 Obese_sd=obese_list$sd,
                                 Nonobese_n=nonobese_list$n,
                                 Nonobese_mean=nonobese_list$mean,
                                 Nonobese_sd=nonobese_list$sd,
                                 Trait=nonobese_list$ID
        )
        
        out$gene_RNA <- gene
        
        out$missing_gene_cohort_RNA <- names(cohort_name[cohort_name %in% setdiff(input$cohort_bmicross,cohort_name[obese_list$ID])])
        
        out_cohort_RNA <- cohort_name[obese_list$ID]
        
        
        if (input$effect_bmicross=="log fold change") {
          
          out_plot_RNA <- meta::metacont(n.e = obese_list$n,
                                         mean.e = obese_list$mean,
                                         sd.e = obese_list$sd,
                                         n.c = nonobese_list$n,
                                         mean.c = nonobese_list$mean,
                                         sd.c = nonobese_list$sd,
                                         studlab = obese_list$ID,
                                         label.e = "obese",
                                         label.c = "nonobese",
                                         label.right = "Obese",
                                         label.left = "Non obese",
                                         sm = "ROM",
                                         backtransf = FALSE)
          
          out_plot_RNA$sm <- "log2 FC"
          out_plot_RNA$TE <- log2(exp(out_plot_RNA$TE))
          out_plot_RNA$TE.common <- log2(exp(out_plot_RNA$TE.common))
          out_plot_RNA$TE.random <- log2(exp(out_plot_RNA$TE.random))
          out_plot_RNA$TE.fixed <- log2(exp(out_plot_RNA$TE.fixed))
          out_plot_RNA$lower <- log2(exp(out_plot_RNA$lower))
          out_plot_RNA$lower.common <- log2(exp(out_plot_RNA$lower.common))
          out_plot_RNA$lower.random <- log2(exp(out_plot_RNA$lower.random))
          out_plot_RNA$upper <- log2(exp(out_plot_RNA$upper))
          out_plot_RNA$upper.common <- log2(exp(out_plot_RNA$upper.common))
          out_plot_RNA$upper.random <- log2(exp(out_plot_RNA$upper.random))
          
        } else if (input$effect_bmicross=="standardized mean difference") {
          
          out_plot_RNA <- meta::metacont(n.e = obese_list$n,
                                         mean.e = obese_list$mean,
                                         sd.e = obese_list$sd,
                                         n.c = nonobese_list$n,
                                         mean.c = nonobese_list$mean,
                                         sd.c = nonobese_list$sd,
                                         studlab = obese_list$ID,
                                         label.e = "obese",
                                         label.c = "nonobese",
                                         label.right = "Obese",
                                         label.left = "Non obese",
                                         sm = "SMD",
                                         method.smd = "Hedges")
          
        }
        
      }
      
      out$has_plot_RNA <- TRUE
      out$plot_RNA <- out_plot_RNA
      out$plottype_RNA <- out_plottype_RNA
      out$df_RNA <- out_df_RNA
      out$cohort_RNA <- out_cohort_RNA
    }
      
    if (input$proteomics_bmicross==T) {
      
      gene_protein <- intersect(unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_bmicross),split='\n|\t|,')) %>% unique(), rownames(Proteomicseset))
      
      if (length(gene_protein)>0) {
        
        if (length(input$cohort_bmicross) > 1) {
          group_bmicross_protein <- input$group_bmicross_protein
        } else {
          group_bmicross_protein <- setdiff(input$group_bmicross_RNA,"Underweight")
        }
        eset <- Proteomicseset
        eset <- eset[,!is.na(eset$BMI_catelogy)]
        
        compare_list <- Filter(function(combo) {
          !(combo[1] == "Obese" & combo[2] %in% c("Obese I", "Obese II", "Obese III")) &
            !(combo[2] == "Obese" & combo[1] %in% c("Obese I", "Obese II", "Obese III")) &
            !(combo[1] == "Non Obese" & combo[2] %in% c("Underweight", "Healthy Weight", "Overweight")) &
            !(combo[2] == "Non Obese" & combo[1] %in% c("Underweight", "Healthy Weight", "Overweight"))
        }, combn(group_bmicross_protein, 2, simplify = FALSE))
        
        
        plot_list <- lapply(gene_protein, function(i) {
          
          # 根据input$group_bmi的选择来确定哪些BMI类别应该被包括
          df_list <- lapply(group_bmicross_protein, function(bmi_choice) {
            switch(bmi_choice,
                   "Underweight" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Underweight"]], BMI_catelogy = "Underweight"),
                   "Healthy Weight" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Healthy Weight"]], BMI_catelogy = "Healthy Weight"),
                   "Overweight" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Overweight"]], BMI_catelogy = "Overweight"),
                   "Obese I" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Obese I"]], BMI_catelogy = "Obese I"),
                   "Obese II" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Obese II"]], BMI_catelogy = "Obese II"),
                   "Obese III" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy == "Obese III"]], BMI_catelogy = "Obese III"),
                   "Non Obese" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")]], BMI_catelogy = "Non Obese"),
                   "Obese" = data.frame(expression = Biobase::exprs(eset)[i, colnames(eset)[eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")]], BMI_catelogy = "Obese")
            )
          })
          
          # 合并所有的数据框
          df <- do.call(rbind, df_list) %>% na.omit()
          
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(length(unique(df$BMI_catelogy)))
          }
          # 为每个基因生成盒形图
          p <- ggboxplot(data = df, x = "BMI_catelogy", y = "expression", fill = "BMI_catelogy", xlab = "", ylab = paste0(i, " expression")) +
            theme(aspect.ratio = 1, 
                  axis.title.x = element_blank(), 
                  axis.text.x = element_text(angle=45, hjust = 1, vjust = 1),
                  text = element_text(family = "Red Hat Display"),
                  legend.position = "none") +
            stat_compare_means(comparisons = compare_list, method = ifelse(test = length(input$cohort_bmicross) == 1,yes = input$statistics_bmicross,no = input$statistics_bmicross_protein)) +
            scale_fill_manual(values = discrete_color)
          
          return(p)
        })
        
        
        # Create a function to get the average expression for a given gene and BMI category
        get_avg_expression <- function(gene, bmi_category) {
          switch(bmi_category,
                 "Underweight" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Underweight"], na.rm = TRUE),
                 "Healthy Weight" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Healthy Weight"], na.rm = TRUE),
                 "Overweight" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Overweight"], na.rm = TRUE),
                 "Obese I" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Obese I"], na.rm = TRUE),
                 "Obese II" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Obese II"], na.rm = TRUE),
                 "Obese III" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy == "Obese III"], na.rm = TRUE),
                 "Non Obese" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy %in% c("Underweight", "Healthy Weight", "Overweight")], na.rm = TRUE),
                 "Obese" = mean(Biobase::exprs(eset)[gene, eset$BMI_catelogy %in% c("Obese I", "Obese II", "Obese III")], na.rm = TRUE)
          )
        }
        
        # Create the output data frame
        out_df_protein <- as.data.frame(matrix(0, nrow = length(gene_protein), ncol = length(group_bmicross_protein)))
        rownames(out_df_protein) <- gene_protein
        colnames(out_df_protein) <- group_bmicross_protein
        
        # Fill the data frame with average expression values
        for (g in gene_protein) {
          for (bmi_cat in group_bmicross_protein) {
            out_df_protein[g, bmi_cat] <- get_avg_expression(g, bmi_cat)
          }
        }
        
        out_plot_protein <- plot_grid(plotlist = plot_list,align = "hv", nrow = ceiling(length(gene_protein)/5))
        
        out_plottype_protein <- "boxplot"
        
        out$row_protein <- ceiling(length(gene_protein)/5)
        
        out$has_plot_protein <- TRUE
        out$plot_protein <- out_plot_protein
        out$plottype_protein <- out_plottype_protein
        out$df_protein <- out_df_protein
        
      }
      
    }
    
    return(out)
    
  }) # eventReactive
  
  
  
  
  ## BMI longitudinal----
  
  candidate_traits_bmilong <- eventReactive(input$cohort_bmilong, {
    # 获取所选队列的所有特征
    all_traits <- lapply(input$cohort_bmilong, function(cohort) {
      colnames(pData(get(cohort)))
    })
    
    # 返回所选队列的特征的并集
    unique(unlist(all_traits))
  })
  
  
  observeEvent(input$cohort_bmilong, {
    updatePickerInput(session = session,
                      inputId = 'trait',
                      label = 'Traits',
                      # selected = '',
                      choices = candidate_traits_bmilong(),
                      # width = "300px",
                      options = pickerOptions(actionsBox = TRUE, size = 10, maxOptions = 100, liveSearch = TRUE),
                      # multiple = FALSE
    )
  }) # observeEvent
  
  
  plotdata_bmilong <- eventReactive(input$SearchButton_bmilong, {
    
    gene <- unlist(strsplit(gsub(pattern = " ",replacement = "",x = input$gene_bmilong),split='\n|\t|,')) %>% unique()
    
    # 使用lapply遍历所有的cohort并获取每个cohort中存在的基因  Use lapply to loop through all cohorts and get the genes present in each cohort
    all_genes_in_cohorts <- lapply(input$cohort_bmilong, function(cohort) {
      rownames(Biobase::exprs(get(cohort)))
    }) %>% unlist() %>% unique()
    
    gene <- intersect(gene, all_genes_in_cohorts)
    
    out <- list(has_plot_RNA = FALSE, plot_RNA = NULL, plottype_RNA = NULL, df_RNA = NULL, cohort_RNA = NULL) # 初始化返回列表
    
    if (length(gene) > 0 | length(input$trait_bmilong) > 0) {
      if (length(gene)>0) {
        plot_list <- lapply(gene, function(current_gene) {
          gene_plots <- lapply(input$cohort_bmilong, function(current_cohort) {
            if (current_gene %in% rownames(Biobase::exprs(get(current_cohort)))) {
              
              temp_pdat <- pData(get(current_cohort))
              temp_pdat$expression <- Biobase::exprs(get(current_cohort))[current_gene,]
              temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
              temp_pdat <- temp_pdat[order(temp_pdat$subject),]
              
              if (length(unique(temp_pdat$group)) > 1) {
                facet_by <- "group"
              } else {
                facet_by <- NULL
              }
              
              if (current_cohort %in% c("GSE141221_diogenes1_WL","GSE95640_diogenes2_WL", "GSE77962eset_WL", "GSE35411eset_WL")) {
                xlab <- "weeks after diet"
              } else {
                xlab <- "years after bariatric surgery"
              }
              
              if (input$discrete_col == "default") {
                discrete_color <- portalcol2
              } else {
                discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(length(unique(temp_pdat$time_point)))
              }
              p <- ggboxplot(data = temp_pdat,x = "time_point",y = "expression",fill = "time_point", xlab = xlab, ylab = paste0(current_gene, " expression"), title = names(cohort_name[cohort_name %in% current_cohort]) %>% gsub(pattern = " weight loss",replacement = "")) +
                # facet_wrap(facets = facet_by, scales = "free_x") +
                geom_line(aes(group=subject), linetype = ifelse(nrow(temp_pdat)/length(unique(temp_pdat$group))/length(unique(temp_pdat$time_point))<50,yes = "dashed",no = NA)) +
                stat_compare_means(comparisons = combn(unique(as.character(temp_pdat$time_point)), 2, simplify = FALSE), paired = TRUE, method = input$statistics_bmilong) +
                theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
                theme(panel.border = element_blank()) +
                theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
                      axis.line.y=element_line(linetype=1,color="black",size=0.9),
                      aspect.ratio = length(unique(temp_pdat$group))* 1.0,
                      text = element_text(family = "Red Hat Display"),
                      strip.background = element_rect(fill = "white", colour = "white"),
                      strip.text = element_text(size = 12),
                      legend.position = "none") +
                scale_fill_manual(values = discrete_color) +
                scale_x_discrete(drop = TRUE)
              
              p <- facet(p, facet.by = facet_by,scales = "free_x")
              
              if (length(unique(temp_pdat$group)) > 1) {
                for (group_num in 1:length(unique(temp_pdat$group))) {
                  p <- p + stat_compare_means(comparisons = list(as.character(unique(temp_pdat$time_point[temp_pdat$group==unique(temp_pdat$group)[group_num]]))), paired = TRUE, method = input$statistics_bmilong)
                }
              }
              
              p
              
            } else {
              ggplot() + theme_void() + theme(aspect.ratio=1)
            }
          })
          
          return(gene_plots)
        })
        
        # 结果是一个嵌套列表，其中每个外部元素对应一个基因，每个内部元素对应一个eset的ggboxplot
        
        # temp <- unlist(plot_list, recursive = FALSE)
        # 
        # 
        # for (i in 1:length(temp)) {
        #   if (i==1) {
        #     a <- paste0('p1 <- egg::ggarrange(','temp[[',i,']]')
        #   } else {
        #     a <- paste0(a,',','temp[[',i,']]')
        #   }
        # }
        # a <- paste0(a,', ncol=length(input$cohort_bmilong), draw = F, newpage = F)')
        # 
        # eval(parse(text=a))
        
        # p1 <- plot_grid(plotlist = unlist(plot_list, recursive = FALSE),align = "hv",ncol = length(input$cohort_bmilong))
        p1 <- unlist(plot_list, recursive = FALSE)
        
        df_1 <- do.call(rbind, lapply(input$cohort_bmilong, function(current_cohort) {
          eset <- get(current_cohort)
          do.call(rbind, lapply(gene, function(current_gene) {
            if (current_gene %in% rownames(Biobase::exprs(eset))) {
              data.frame(
                cohort = names(cohort_name[cohort_name == current_cohort]),
                time_point = unique(pData(eset)$time_point),
                gene_or_trait = current_gene,
                mean_value = sapply(unique(pData(eset)$time_point), function(tp) {
                  mean(Biobase::exprs(eset)[current_gene, pData(eset)$time_point == tp], na.rm = TRUE)
                })
              )
            } else {
              NULL
            }
          }))
        }))
        
      }
      
      if (length(input$trait_bmilong)>0) {
        plot_list <- lapply(input$trait_bmilong, function(current_trait) {
          trait_plots <- lapply(input$cohort_bmilong, function(current_cohort) {
            if (current_trait %in% colnames(pData(get(current_cohort)))) {
              
              temp_pdat <- pData(get(current_cohort))
              temp_pdat <- temp_pdat[order(temp_pdat$time_point),]
              temp_pdat <- temp_pdat[order(temp_pdat$subject),]
              
              if (length(unique(temp_pdat$group)) > 1) {
                facet_by <- "group"
              } else {
                facet_by <- NULL
              }
              
              if (current_cohort %in% c("GSE141221_diogenes1_WL","GSE95640_diogenes2_WL", "GSE77962eset_WL", "GSE35411eset_WL")) {
                xlab <- "weeks after diet"
              } else {
                xlab <- "years after bariatric surgery"
              }
              
              if (input$discrete_col == "default") {
                discrete_color <- portalcol2
              } else {
                discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(length(unique(temp_pdat$time_point)))
              }
              p <- ggboxplot(data = temp_pdat,x = "time_point",y = current_trait,fill = "time_point",xlab = xlab,ylab = current_trait, title = names(cohort_name[cohort_name %in% current_cohort]) %>% gsub(pattern = " weight loss",replacement = "")) +
                # facet_wrap(facets = facet_by, scales = "free_x") +
                geom_line(aes(group=subject), linetype = ifelse(nrow(temp_pdat)/length(unique(temp_pdat$group))/length(unique(temp_pdat$time_point))<50,yes = "dashed",no = NA)) +
                stat_compare_means(comparisons = combn(unique(as.character(temp_pdat$time_point)), 2, simplify = FALSE), paired = TRUE, method = input$statistics_bmilong) +
                theme(panel.grid.major=element_line(colour=NA),panel.grid.minor=element_line(colour=NA)) +
                theme(panel.border = element_blank()) +
                theme(axis.line.x=element_line(linetype=1,color="black",size=0.9),
                      axis.line.y=element_line(linetype=1,color="black",size=0.9),
                      aspect.ratio = length(unique(temp_pdat$group))* 1.0,
                      text = element_text(family = "Red Hat Display"),
                      strip.background = element_rect(fill = "white", colour = "white"),
                      strip.text = element_text(size = 12),
                      legend.position = "none") +
                scale_fill_manual(values = discrete_color) +
                scale_x_discrete(drop = TRUE)
              
              p <- facet(p, facet.by = facet_by,scales = "free_x")
              
              if (length(unique(temp_pdat$group)) > 1) {
                for (group_num in 1:length(unique(temp_pdat$group))) {
                  p <- p + stat_compare_means(comparisons = list(as.character(unique(temp_pdat$time_point[temp_pdat$group==unique(temp_pdat$group)[group_num]]))), paired = TRUE, method = input$statistics_bmilong)
                }
              }
              
              p
              
            } else {
              ggplot() + theme_void() + theme(aspect.ratio=1)
            }
          })
          
          return(trait_plots)
        })
        
        # 结果是一个嵌套列表，其中每个外部元素对应一个基因，每个内部元素对应一个eset的ggboxplot
        
        # p2 <- plot_grid(plotlist = unlist(plot_list, recursive = FALSE),align = "hv",ncol = length(input$cohort_bmilong))
        p2 <- unlist(plot_list, recursive = FALSE)
        
        df_2 <- do.call(rbind, lapply(input$cohort_bmilong, function(current_cohort) {
          eset <- get(current_cohort)
          do.call(rbind, lapply(input$trait_bmilong, function(current_trait) {
            if (current_trait %in% colnames(pData(eset))) {
              data.frame(
                cohort = current_cohort,
                time_point = unique(pData(eset)$time_point),
                gene_or_trait = current_trait,
                mean_value = sapply(unique(pData(eset)$time_point), function(tp) {
                  mean(pData(eset)[pData(eset)$time_point == tp, current_trait], na.rm = TRUE)
                })
              )
            } else {
              NULL
            }
          }))
        }))
        
      }
      
      if (length(gene)>0 & length(input$trait_bmilong)>0) {
        # out$plot_RNA <- plot_grid(plotlist = list(p1,p2),align = "hv",ncol = 1, rel_heights = c(1/length(input$trait_bmilong),1/length(gene)))
        out$plot_RNA <- c(p1, p2)
        out$df_RNA <- rbind(df_1,df_2)
      } else if (length(gene)>0 & length(input$trait_bmilong)==0) {
        out$plot_RNA <- p1
        out$df_RNA <- df_1
      } else if (length(gene)==0 & length(input$trait_bmilong)>0) {
        out$plot_RNA <- p2
        out$df_RNA <- df_2
      }
      
      out$has_plot_RNA <- TRUE

      out$row_RNA <- length(gene) + length(input$trait_bmilong)
      
      out$col_RNA <- length(input$cohort_bmilong)
      
      out$cohort_RNA <- cohort_name[cohort_name %in% input$cohort_bmilong]    
    }
    return(out)
    
  }) # eventReactive
  
  
  
  ## output phenotype association ----
  
  output$ui_plot_pc_RNA <- renderUI({if (input$SearchButton_pc) {
    if (plotdata_pc()$has_plot_RNA) {
      output$plot_pc_RNA <- renderPlot({
        if (plotdata_pc()$plottype_RNA %in% c("heatmap", "barplot")) {
          plotdata_pc()$plot_RNA
        } else if (plotdata_pc()$plottype_RNA == "metacor") {
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2[1:2]
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
          }
          # par(mar = c(1, 1, 1, 1))
          # meta::forest(plotdata_pc()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
          
          if (is.na(plotdata_pc()$plot_RNA$pval.Q)) {
            meta::forest(plotdata_pc()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
          } else {
            if (plotdata_pc()$plot_RNA$pval.Q<0.05) {
              if (plotdata_pc()$plot_RNA$TE.random>0) {
                meta::forest(plotdata_pc()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              } else {
                meta::forest(plotdata_pc()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              }
            } else {
              if (plotdata_pc()$plot_RNA$TE.common>0) {
                meta::forest(plotdata_pc()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              } else {
                meta::forest(plotdata_pc()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              }
            }
          }          
        }
      })
      
      output$pdf_pc_RNA <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          if (plotdata_pc()$plottype_RNA %in% c("heatmap", "barplot")) {
            ggsave(filename = fname, plot = plotdata_pc()$plot_RNA, height = ifelse(test = plotdata_pc()$plottype_RNA == "heatmap", yes = max(nrow(plotdata_pc()$df_RNA)*5/17+45/17,5), no = 5), width = 10, units = "in", device = cairo_pdf)
          }  else if (plotdata_pc()$plottype_RNA == "metacor") {
            cairo_pdf(file = fname, width=10, height=length(plotdata_pc()$plot_RNA$studlab)*0.2+2.4)
            if (input$discrete_col == "default") {
              discrete_color <- portalcol2[1:2]
            } else {
              discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
            }
            if (is.na(plotdata_pc()$plot_RNA$pval.Q)) {
              meta::forest(plotdata_pc()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              if (plotdata_pc()$plot_RNA$pval.Q<0.05) {
                if (plotdata_pc()$plot_RNA$TE.random>0) {
                  meta::forest(plotdata_pc()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                } else {
                  meta::forest(plotdata_pc()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                }
              } else {
                if (plotdata_pc()$plot_RNA$TE.common>0) {
                  meta::forest(plotdata_pc()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                } else {
                  meta::forest(plotdata_pc()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                }
              }
            }
            dev.off()
          }
        }
      )
      
      Transcriptomics <- if (is.null(plotdata_pc()$plot_protein)) {
        NULL
      } else {
        "Transcriptomics:"
      }
      
      
      return(list(Transcriptomics,
                  plotOutput("plot_pc_RNA",height = ifelse(test = plotdata_pc()$plottype_RNA == "heatmap", 
                                                           yes = paste0(max(nrow(plotdata_pc()$df_RNA)*20+210,400),'px'), 
                                                           no = ifelse(test = plotdata_pc()$plottype_RNA == "metacor",
                                                                       yes = paste0(max(length(plotdata_pc()$plot_RNA$studlab)*15+150,400),'px'),
                                                                       no = "400px"))),
                  br(),
                  downloadButton(outputId = 'pdf_pc_RNA', label = "Download pdf")
      )
      )      
    } else {
      return(NULL)
    }
  }})
  output$ui_plot_pc_protein <- renderUI({if (input$SearchButton_pc) {
    if (plotdata_pc()$has_plot_protein) {
      output$plot_pc_protein <- renderPlot({
        plotdata_pc()$plot_protein
      })
      
      output$pdf_pc_protein <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          ggsave(filename = fname, plot = plotdata_pc()$plot_protein, height = 10, width = 10, units = "in", device = cairo_pdf)
        }
      )
      
      return(list("Proteomics:",
                  plotOutput("plot_pc_protein",height = ifelse(test = plotdata_pc()$plottype_protein == "heatmap", yes = paste0(max(nrow(plotdata_pc()$df_protein)*20+210,400),'px'), no = "400px")),
                  br(),
                  downloadButton(outputId = 'pdf_pc_protein', label = "Download pdf")
      )
      )      
    } else {
      NULL
    }
  }})
  output$ui_df_pc_RNA <- renderUI({if (input$SearchButton_pc) {
    if (plotdata_pc()$has_plot_RNA) {
      output$df_pc_RNA <- DT::renderDataTable(server = FALSE,{return(plotdata_pc()$df_RNA %>% mutate_if(is.numeric, round, digits = 4))},
                                              extensions = c('Buttons'),
                                              options = list(scrollX = TRUE,
                                                             pageLength = 10,
                                                             lengthMenu = c(10, 25, 50, 100),
                                                             dom = 'Blfrtip',
                                                             buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                              )
      )
      
      Transcriptomics <- if (is.null(plotdata_pc()$df_protein)) {
        NULL
      } else {
        "Transcriptomics:"
      }
      
      return(list(Transcriptomics,
                  DT::dataTableOutput("df_pc_RNA")))      
    } else {
      NULL
    }    
  }})
  output$ui_df_pc_protein <- renderUI({if (input$SearchButton_pc) {
    if (plotdata_pc()$has_plot_protein) {
      output$df_pc_protein <- DT::renderDataTable(server = FALSE,{return(plotdata_pc()$df_protein %>% mutate_if(is.numeric, round, digits = 4))},
                                                  extensions = c('Buttons'),
                                                  options = list(scrollX = TRUE,
                                                                 pageLength = 10,
                                                                 lengthMenu = c(10, 25, 50, 100),
                                                                 dom = 'Blfrtip',
                                                                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                                  )
      )
      return(list("Proteomics:",
                  DT::dataTableOutput("df_pc_protein")))      
    } else {
      NULL
    }
  }})
  output$ui_citation_pc <- renderUI({if (input$SearchButton_pc) {
    if (plotdata_pc()$has_plot_RNA) {
      citation_temp <- citation[citation$Cohort_ID %in% plotdata_pc()$cohort_RNA,3:9,drop=F] %>% unique()
      
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
  
  
  
  
  ## output Sex ----
  
  output$ui_plot_sex_RNA <- renderUI({if (input$SearchButton_sex) {
    if (plotdata_sex()$has_plot_RNA) {
      output$plot_sex_RNA <- renderPlot({
        if (plotdata_sex()$plottype_RNA %in% c("heatmap", "boxplot")) {
          plotdata_sex()$plot_RNA
        } else if (plotdata_sex()$plottype_RNA %in% c("metacont","metagen")) {
          # par(mar = c(1, 1, 1, 1))
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2[1:2]
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
          }
          if (is.na(plotdata_sex()$plot_RNA$pval.Q)) {
            meta::forest(plotdata_sex()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
          } else {
            if (plotdata_sex()$plot_RNA$pval.Q<0.05) {
              if (plotdata_sex()$plot_RNA$TE.random>0) {
                meta::forest(plotdata_sex()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              } else {
                meta::forest(plotdata_sex()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              }
            } else {
              if (plotdata_sex()$plot_RNA$TE.common>0) {
                meta::forest(plotdata_sex()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              } else {
                meta::forest(plotdata_sex()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              }
            }
          }
        }
      })
      output$pdf_sex_RNA <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          if (plotdata_sex()$plottype_RNA %in% c("heatmap", "boxplot")) {
            ggsave(filename = fname, plot = plotdata_sex()$plot_RNA, height = ifelse(test = plotdata_sex()$plottype_RNA == "heatmap", yes = max(nrow(plotdata_sex()$df_RNA)*5/17+4,5), no = 5), width = 12, units = "in", device = cairo_pdf)
          } else if (plotdata_sex()$plottype_RNA == "metacont") {
            cairo_pdf(file = fname, width=13, height=length(plotdata_sex()$plot_RNA$studlab)*0.2+2.4)
            if (input$discrete_col == "default") {
              discrete_color <- portalcol2[1:2]
            } else {
              discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
            }
            if (is.na(plotdata_sex()$plot_RNA$pval.Q)) {
              meta::forest(plotdata_sex()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              if (plotdata_sex()$plot_RNA$pval.Q<0.05) {
                if (plotdata_sex()$plot_RNA$TE.random>0) {
                  meta::forest(plotdata_sex()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                } else {
                  meta::forest(plotdata_sex()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                }
              } else {
                if (plotdata_sex()$plot_RNA$TE.common>0) {
                  meta::forest(plotdata_sex()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                } else {
                  meta::forest(plotdata_sex()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                }
              }
            }
            dev.off()
          }
        }
      )
      
      Transcriptomics <- if (is.null(plotdata_sex()$plot_protein)) {
        NULL
      } else {
        "Transcriptomics:"
      }
      
      return(list(Transcriptomics,
                  plotOutput("plot_sex_RNA", height = ifelse(test = plotdata_sex()$plottype_RNA == "heatmap", 
                                                             yes = paste0(max(nrow(plotdata_sex()$df_RNA)*20+260,400),'px'), 
                                                             no = ifelse(test = plotdata_sex()$plottype_RNA == "boxplot", 
                                                                         yes = paste0(plotdata_sex()$row_RNA*300,'px'), 
                                                                         no = ifelse(test = plotdata_sex()$plottype_RNA == "metacont",
                                                                                     yes = paste0(max(length(plotdata_sex()$plot_RNA$studlab)*15+150,400),'px'),
                                                                                     no = "400px")))),
                  br(),
                  downloadButton(outputId = 'pdf_sex_RNA', label = "Download pdf")
      )
      )       
    } else {
      NULL
    }
    
    
  }})
  output$ui_plot_sex_protein <- renderUI({if (input$SearchButton_sex) {
    if (plotdata_sex()$has_plot_protein) {
      output$plot_sex_protein <- renderPlot({
        plotdata_sex()$plot_protein
      })
      
      output$pdf_sex_protein <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          ggsave(filename = fname, plot = plotdata_sex()$plot_protein, height = 10, width = 10, units = "in", device = cairo_pdf)
        }
      )
      
      return(list("Proteomics:",
                  plotOutput("plot_sex_protein", height = paste0(plotdata_sex()$row_protein*300,'px')),
                  br(),
                  downloadButton(outputId = 'pdf_sex_protein', label = "Download pdf")
      )
      )      
    } else {
      NULL
    }
    
  }})
  output$ui_df_sex_RNA <- renderUI({if (input$SearchButton_sex) {
    if (plotdata_sex()$has_plot_RNA) {
      output$df_sex_RNA <- DT::renderDataTable(server = FALSE,{return(plotdata_sex()$df_RNA %>% mutate_if(is.numeric, round, digits = 4))},
                                               extensions = c('Buttons'),
                                               options = list(scrollX = TRUE,
                                                              pageLength = 10,
                                                              lengthMenu = c(10, 25, 50, 100),
                                                              dom = 'Blfrtip',
                                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                               )
      )
      
      Transcriptomics <- if (is.null(plotdata_sex()$df_protein)) {
        NULL
      } else {
        "Transcriptomics:"
      }
      
      return(list(Transcriptomics,
                  DT::dataTableOutput("df_sex_RNA")))
    } else {
      NULL
    }
    
  }})
  output$ui_df_sex_protein <- renderUI({if (input$SearchButton_sex) {
    if (plotdata_sex()$has_plot_protein) {
      output$df_sex_protein <- DT::renderDataTable(server = FALSE,{return(plotdata_sex()$df_protein %>% mutate_if(is.numeric, round, digits = 4))},
                                                   extensions = c('Buttons'),
                                                   options = list(scrollX = TRUE,
                                                                  pageLength = 10,
                                                                  lengthMenu = c(10, 25, 50, 100),
                                                                  dom = 'Blfrtip',
                                                                  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                                   )
      )
      
      
      return(list("Proteomics:",
                  DT::dataTableOutput("df_sex_protein")))      
    } else {
      NULL
    }
  }})
  output$ui_citation_sex <- renderUI({if (input$SearchButton_sex) {
    if (plotdata_sex()$has_plot_RNA) {
      citation_temp <- citation[citation$Cohort_ID %in% plotdata_sex()$cohort_RNA,3:9] %>% unique()
      
      for (i in 1:nrow(citation_temp)) {
        if (i==1) {
          citation_list <- list()
        }
        
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
      NULL
    }
  }})
  
  
  
  
  ## output BMI cross-sectional----
  
  output$ui_plot_bmicross_RNA <- renderUI({if (input$SearchButton_bmicross) {
    
    if (plotdata_bmicross()$has_plot_RNA) {
      output$plot_bmicross_RNA <- renderPlot({
        if (plotdata_bmicross()$plottype_RNA %in% c("heatmap", "boxplot")) {
          plotdata_bmicross()$plot_RNA
        } else if (plotdata_bmicross()$plottype_RNA %in% c("metacont","metagen")) {
          # par(mar = c(1, 1, 1, 1))
          if (input$discrete_col == "default") {
            discrete_color <- portalcol2[1:2]
          } else {
            discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
          }
          if (is.na(plotdata_bmicross()$plot_RNA$pval.Q)) {
            meta::forest(plotdata_bmicross()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
          } else {
            if (plotdata_bmicross()$plot_RNA$pval.Q<0.05) {
              if (plotdata_bmicross()$plot_RNA$TE.random>0) {
                meta::forest(plotdata_bmicross()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              } else {
                meta::forest(plotdata_bmicross()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              }
            } else {
              if (plotdata_bmicross()$plot_RNA$TE.common>0) {
                meta::forest(plotdata_bmicross()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              } else {
                meta::forest(plotdata_bmicross()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
              }
            }
          }
        }
      })
      output$pdf_bmicross_RNA <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          if (plotdata_bmicross()$plottype_RNA %in% c("heatmap", "boxplot")) {
            ggsave(filename = fname, plot = plotdata_bmicross()$plot_RNA, height = ifelse(test = plotdata_bmicross()$plottype_RNA == "heatmap", yes = max(nrow(plotdata_bmicross()$df_RNA)*5/17+4,5), no = 5), width = 12, units = "in", device = cairo_pdf)
          } else if (plotdata_bmicross()$plottype_RNA == "metacont") {
            cairo_pdf(file = fname, width=13, height=length(plotdata_bmicross()$plot_RNA$studlab)*0.2+2.4)
            if (input$discrete_col == "default") {
              discrete_color <- portalcol2[1:2]
            } else {
              discrete_color <- colorRampPalette(brewer.pal(n = 7, name = input$discrete_col_CB2))(2)
            }
            if (is.na(plotdata_bmicross()$plot_RNA$pval.Q)) {
              meta::forest(plotdata_bmicross()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            } else {
              if (plotdata_bmicross()$plot_RNA$pval.Q<0.05) {
                if (plotdata_bmicross()$plot_RNA$TE.random>0) {
                  meta::forest(plotdata_bmicross()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                } else {
                  meta::forest(plotdata_bmicross()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                }
              } else {
                if (plotdata_bmicross()$plot_RNA$TE.common>0) {
                  meta::forest(plotdata_bmicross()$plot_RNA, sortvar = -TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                } else {
                  meta::forest(plotdata_bmicross()$plot_RNA, sortvar = TE, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
                }
              }
            }
            dev.off()
          }
        }
      )
      
      Transcriptomics <- if (is.null(plotdata_bmicross()$plot_protein)) {
        NULL
      } else {
        "Transcriptomics:"
      }
      
      return(list(Transcriptomics,
                  plotOutput("plot_bmicross_RNA",height = ifelse(test = plotdata_bmicross()$plottype_RNA == "heatmap", 
                                                                 yes = paste0(max(nrow(plotdata_bmicross()$df_RNA)*20+260,400),'px'), 
                                                                 no = ifelse(test = plotdata_bmicross()$plottype_RNA == "boxplot", 
                                                                             yes = paste0(plotdata_bmicross()$row_RNA*300,'px'), 
                                                                             no = ifelse(test = plotdata_bmicross()$plottype_RNA == "metacont",
                                                                                         yes = paste0(max(length(plotdata_bmicross()$plot_RNA$studlab)*15+150,400),'px'),
                                                                                         no = "400px")))),
                  br(),
                  downloadButton(outputId = 'pdf_bmicross_RNA', label = "Download pdf")
      )
      )       
    } else {
      NULL
    }
  }})
  output$ui_plot_bmicross_protein <- renderUI({if (input$SearchButton_bmicross) {
    
    if (plotdata_bmicross()$has_plot_protein) {
      output$plot_bmicross_protein <- renderPlot({
        plotdata_bmicross()$plot_protein
      })
      
      output$pdf_bmicross_protein <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          ggsave(filename = fname, plot = plotdata_bmicross()$plot_protein, height = 10, width = 10, units = "in", device = cairo_pdf)
        }
      )
      
      return(list("Proteomics:",
                  plotOutput("plot_bmicross_protein", height = paste0(plotdata_bmicross()$row_protein*300,'px')),
                  br(),
                  downloadButton(outputId = 'pdf_bmicross_protein', label = "Download pdf")
      )
      )      
    } else {
      NULL
    }
  }})
  output$ui_df_bmicross_RNA <- renderUI({if (input$SearchButton_bmicross) {
    if (plotdata_bmicross()$has_plot_RNA) {
      output$df_bmicross_RNA <- DT::renderDataTable(server = FALSE,{return(plotdata_bmicross()$df_RNA %>% mutate_if(is.numeric, round, digits = 4))},
                                                    extensions = c('Buttons'),
                                                    options = list(scrollX = TRUE,
                                                                   pageLength = 10,
                                                                   lengthMenu = c(10, 25, 50, 100),
                                                                   dom = 'Blfrtip',
                                                                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                                    )
      )
      
      Transcriptomics <- if (is.null(plotdata_bmicross()$plot_protein)) {
        NULL
      } else {
        "Transcriptomics:"
      }
      
      return(list(Transcriptomics,
                  DT::dataTableOutput("df_bmicross_RNA")))       
    } else {
      NULL
    }
  }})
  output$ui_df_bmicross_protein <- renderUI({if (input$SearchButton_bmicross) {
    
    if (plotdata_bmicross()$has_plot_protein) {
      NULL
    } else {
      output$df_bmicross_protein <- DT::renderDataTable(server = FALSE,{return(plotdata_bmicross()$df_protein %>% mutate_if(is.numeric, round, digits = 4))},
                                                        extensions = c('Buttons'),
                                                        options = list(scrollX = TRUE,
                                                                       pageLength = 10,
                                                                       lengthMenu = c(10, 25, 50, 100),
                                                                       dom = 'Blfrtip',
                                                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                                        )
      )
      
      
      return(list("Proteomics:",
                  DT::dataTableOutput("df_bmicross_protein")))
    }
    
  }})
  output$ui_citation_bmicross <- renderUI({if (input$SearchButton_bmicross) {
    
    if (plotdata_bmicross()$has_plot_RNA) {
      citation_temp <- citation[citation$Cohort_ID %in% plotdata_bmicross()$cohort_RNA,3:9] %>% unique()
      
      for (i in 1:nrow(citation_temp)) {
        if (i==1) {
          citation_list <- list()
        }
        
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
  
  
  
  
  ## output BMI longitudinal----
  plot_ggsave <- reactiveVal()
  
  output$ui_plot_bmilong <- renderUI({if (input$SearchButton_bmilong) {
    if (plotdata_bmilong()$has_plot_RNA) {
      output$plot_bmilong <- renderPlot({
        if (class(plotdata_bmilong()$plot_RNA)[1] == "list") {
          # plotdata_bmilong()$plot_RNA
          egg::ggarrange(plots = plotdata_bmilong()$plot_RNA, ncol = plotdata_bmilong()$col_RNA)
          plot_ggsave(egg::ggarrange(plots = plotdata_bmilong()$plot_RNA, ncol = plotdata_bmilong()$col_RNA))
        } else if (class(plotdata_bmilong()$plot_RNA)[1] == "metacor") {
          # par(mar = c(1, 1, 1, 1))
          meta::forest(plotdata_bmilong()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
        }
      })
      
      output$pdf_bmilong <- downloadHandler(
        filename = function() {
          paste0(Sys.time(),".pdf")
        },
        content = function(fname) {
          if (class(plotdata_bmilong()$plot_RNA)[1] == "list") {
            ggsave(filename = fname, plot = plot_ggsave(), height = 10, width = 10, units = "in", device = cairo_pdf)
          } else if (class(plotdata_bmilong()$plot_RNA)[1] == "metacor") {
            cairo_pdf(file = fname, width=13, height=length(plotdata_bmilong()$plot_RNA$studlab)*0.2+2.4)
            meta::forest(plotdata_bmilong()$plot_RNA, test.overall.fixed = T, test.overall.random = T, squaresize = 1, digits.pval = 4, fontfamily = "Red Hat Display" ,col.square = discrete_color[1], col.diamond.common = discrete_color[2], col.diamond.random = discrete_color[2], col.square.lines = "#000000")
            dev.off()
          }
        }
      )
      
      return(list(plotOutput("plot_bmilong",height = paste0(plotdata_bmilong()$row_RNA*300,"px")),
                  br(),
                  downloadButton(outputId = 'pdf_bmilong', label = "Download pdf")
      )
      )       
    } else {
      NULL
    }
    
  }})
  output$ui_df_bmilong <- renderUI({if (input$SearchButton_bmilong) {
    if (plotdata_bmilong()$has_plot_RNA) {
      output$df_bmilong <- DT::renderDataTable(server = FALSE,{return(plotdata_bmilong()$df_RNA %>% mutate_if(is.numeric, round, digits = 4))},
                                               extensions = c('Buttons'),
                                               options = list(scrollX = TRUE,
                                                              pageLength = 10,
                                                              lengthMenu = c(10, 25, 50, 100),
                                                              dom = 'Blfrtip',
                                                              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                                               )
      )
      
      return(DT::dataTableOutput("df_bmilong"))
    } else {
      NULL
    }
  }})
  output$ui_citation_bmilong <- renderUI({if (input$SearchButton_bmilong) {
    if (plotdata_bmilong()$has_plot_RNA) {
      citation_temp <- citation[citation$Cohort_ID %in% plotdata_bmilong()$cohort_RNA,3:9,drop=F] %>% unique()
      
      for (i in 1:nrow(citation_temp)) {
        if (i==1) {
          citation_list <- list()
        }
        
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
  
}) # shinyServer

