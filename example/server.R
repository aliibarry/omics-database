##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

# load data
data_dir = "data/database.RData"
scRNA_dir = "data/drg.combined.rds"
theme_dir = "labels.R"

load(data_dir)
source(theme_dir) # include the theme and labels

df_list = list(TPM_subtype,TPM_mouse,TPM_rat, TPM_ipsc, TPM_HS_diabetes, TPM_HS_CTS,TPM_zheng, TPM_humandrg)
col_list = list(bulkseq_colData,TPM_mouse_colData,TPM_rat_colData,ipsc_colData,db_colData,skin_colData,zheng_colData, humandrg_colData) 
pain_list = list("ipsi", "SNI", "SNT", "patient", "PDPN", "pre_Surgery", "", "P")
species_list = list("mouse", "mouse", "rat", "human", "human", "human", "mouse", "human")

credentials <- data.frame(
  user = c("shiny", "ndcn"),
  password = c("azerty", "ndcn-rnaseq"),
  stringsAsFactors = FALSE
)


generate_degdf <- function(deg_df_list, pop_list, dataset) {
  i=1
  final_df = data.frame()
  for (df in deg_df_list) {
    # for the young and old dataset
    if ("hsymbol" %in% colnames(df)){
      df$symbol = df$hsymbol
    }
    df = df[c("log2FoldChange", "padj", "symbol")]
    colnames(df) = c("LFC", "padj", "symbol")
    df$Population = rep(pop_list[i], nrow(df))
    final_df = rbind(final_df, df)
    i=i+1
  }
  if (!is.null(dataset)) {
    final_df$Dataset = rep(dataset, nrow(final_df))
  }
  final_df = mutate(final_df, sig=ifelse(final_df$padj<0.05, "SIG", "NS"))
  return(final_df)
}



generate_combine_dataset <- function(df_list, genes, sex, col_list, pain_list, species_list) {
  final_df = data.frame()
  i = 1
  for (df in df_list) {
    
    # filter
    if (('mgi_symbol' %in% colnames(df))) {
      filt = df[df[,'mgi_symbol'] %in% genes,]
    }
    else{
      filt = df[df[,'symbol'] %in% genes,]
    }
    
    # remove unwanted columns
    if ('mgi_symbol' %in% colnames(filt)) {
      matfilt = subset(filt, select=-c(symbol, mgi_symbol))
    }
    else{
      matfilt = subset(filt, select=-c(symbol))
    }
    
    # merge
    tcounts <- t(matfilt) %>%
      base::merge(as.data.frame(col_list[i]), ., by="row.names") %>%
      gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
    
    species = species_list[i]
    if (species == 'mouse') {
      tcounts$symbol <- mouse_gene_data[tcounts$gene,]$mgi_symbol
    }
    if (species == "rat") {
      tcounts$symbol <- rat_gene_data[tcounts$gene,]$rgd_symbol
    }
    if (species == "human") {
      tcounts$symbol <- human_gene_data[tcounts$gene,]$hgnc_symbol
    }
    ifelse(sex=="Both", tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Dataset, Species) %>%
             dplyr::summarise(expression=median(as.double(expression))), tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Sex, Dataset, Species) %>%
             dplyr::summarise(expression=median(as.double(expression))))
    tcounts_med <- tcounts_med[!tcounts_med$Condition %in% pain_list[i], ]
    final_df = rbind(tcounts_med, final_df)
    i = i+1
  }
  return(final_df)
}


shinydegplot <- function(input, output, session, deg_df_list, pop_list, dataset, genes) {
  
  deg_df = generate_degdf(deg_df_list, pop_list, dataset)
  
  deg <- reactive({
    df <- deg_df %>% filter(symbol %in% genes)
    return(df)
  }) %>% bindCache(genes)
  
  output$deg_plot <- renderPlotly({
    ggplot(data = deg(), aes(x=interaction(Population), y = symbol,
                             text = paste('padj: ',padj))) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      geom_point(aes(col=LFC, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
      th + labs(shape = "", size = "") + facet_wrap(~Dataset, scales = "free_x") + scale_x_discrete(labels=subpopulation_labels)
  })
  
}

# plot the combined deg plots
deg_combine_server <- function(id, datar, datam, datah) {
  moduleServer(id, function(input, output, session){
    df_rat = datar()
    df_mouse = datam()
    df_human = datah()
    
    g1 = ggplot(data = df_rat, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) +
      geom_point(aes(col=LFC, shape=sig, size=0.3)) + ggtitle("Rat") +
      labs(shape = "", size = "") + facet_grid(.~Dataset, scale = "free", space='free')  +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) +
      scale_colour_continuous(limits=c(0,1)) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    
    g2 = ggplot(data = df_human, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) + ggtitle("Human") +
      geom_point(aes(col=LFC, shape=sig, size=0.3)) + facet_grid(.~Dataset, scale = "free", space='free') +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    g3 = ggplot(data = df_mouse, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) +
      ggtitle("Mouse") + geom_point(aes(col=LFC, shape=sig, size=0.3)) +
      facet_grid(.~Dataset, scale = "free", space='free') + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "") + guides(size=FALSE)
    
    combine_degplot = plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
    
    output$deg <- renderCachedPlot({combine_degplot},
                                   cacheKeyExpr = {list(datar, datah, datam)})
    return(combine_degplot)
    
    
  })
}

plotcombine_server <- function(id, df, sex, genes) {
  moduleServer(id, function(input, output, session) {
    df_rat = df[df$Species == "Rat (DRG)",]
    df_mouse = df[df$Species == "Mouse (DRG)",]
    df_human = df[df$Species == "human",]
    
    g1 = ggplot(df_rat, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) + scale_size_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(col = FALSE, size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Rat") + scale_x_discrete(labels=population_labels)
    g2 = ggplot(df_human, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(col = FALSE, size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Human") + scale_x_discrete(labels=population_labels)
    g3 = ggplot(df_mouse, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Mouse") +
      scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
    
    plot = plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
    
    if (sex == "Separate") {
      g1 = g1 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g2 = g2 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g3 = g3 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      plot = plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
    }
    
    output$dot <- renderCachedPlot({plot}, cacheKeyExpr = {list(df,genes)})
    
    return(plot)
    
  })
}

shinyscrna <- function(input, output, session,genes,pbmc, parent=NULL) {
  
  if (is.null(parent) == FALSE) {
    observeEvent(input$link_to_home, {
      newvalue <- "tabhome"
      observe({
        updateTabItems(session=parent, "tabs", newvalue)
      })
    })
  }
  
  observe({
    updateSelectizeInput(session,
                         inputId = "scgeneid",
                         label = NULL,
                         choices = rownames(pbmc),
                         server = TRUE,
                         selected = "ATF3"
    )
  })
  
  observeEvent(input$load, {
    hg = subset(human_gene_data, mgi_symbol %in% genes())$hgnc_symbol
    dotplot = DotPlot(pbmc, features = hg) +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(size=8, angle = 45, hjust= 1),
            axis.text.y = element_text(size=8),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(), legend.justification = c(0,0.3),
            legend.title = element_text(size=10), legend.key.size = unit(0.4, "cm"))
    output$scrna_dots <- renderPlot({
      dotplot
    })
    
    dimplot = DimPlot(pbmc,label=TRUE)
    output$scrna_umap <- renderCachedPlot({dimplot}, cacheKeyExpr = {list(pbmc)})
    
    # download plots and tables in a zip file   
    output$scrna_plots <- downloadHandler(
      filename = function() {
        paste("scrna_plots", "zip", sep=".")
      },
      
      content = function(file){
        # Set temporary working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        fs = c('dotplot.png', 'dimplot.png')
        
        # Save the plots 
        ggsave('dotplot.png', plot = dotplot, width = 6, height = 6, dpi = 300, units = "in", device='png')
        ggsave('dimplot.png', plot = dimplot, width = 6, height = 6, dpi = 300, units = "in", device='png')

        # Zip them up
        zip(file, fs)
      }
    )
    
  })
  
  observeEvent(input$scgeneid,{
    output$scrna_feature <- renderPlotly({
      FeaturePlot(pbmc, features = input$scgeneid)
    })
  })
  
}

plothomescdot_server <- function(id, pbmc, genes) {
  moduleServer(id, function(input, output, session) {
    hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
    g = DotPlot(pbmc, features = hg) + thc + ggtitle("Human DRG subtypes", subtitle="Spatial-seq (Tavares-Ferreira 2021)")
    output$home_scrna_dots <- renderCachedPlot({g}, cacheKeyExpr = {list(genes, pbmc)})
    return(g)
    
  })
}

shinytab <- function(input, output, session, sex, count_df, colData, gene_data, species, paintype, genes,
                     include_deg=FALSE, deg_df_list=NULL, pop_list=NULL, dataset=NULL, include_subtype=FALSE,
                     injury_data = NULL,parent=NULL) {
  if (is.null(parent) == FALSE) {
    observeEvent(input$link_to_home, {
      newvalue <- "tabhome"
      observe({
        updateTabItems(session=parent, "tabs", newvalue)
      })
    })
  }
  output$report <- downloadHandler(
    
    # For PDF output, change this to "report.pdf"
    filename = "RNAseq_plot.html",
    content = function(file) {
      
      if (species == 'mouse') {
        gene_data = mouse_gene_data
      }
      if (species == "rat") {
        gene_data = rat_gene_data
      }
      if (species == "human") {
        gene_data = human_gene_data
      }
      
      params <- list(matrix = count_df, colD = colData, species = species, paintype = paintype, genes = genes(), sex= sex, injury_data = injury_data, 
                     include_deg = include_deg, deg_df_list = deg_df_list, pop_list = pop_list, gene_data = gene_data)
      
      shiny::withProgress(value = 0,
                          message = 'Rendering plotting report',
                          detail =  'This might take several minutes.',{
                            rmarkdown::render(input = "report/RNAseq_plot.Rmd", 
                                              params = params, 
                                              output_file = file,
                                              envir = new.env(parent = globalenv()))
                          })
      
    }
  )
  
  observeEvent(input$load, {
    
    if ('mgi_symbol' %in% colnames(count_df)) {
      df = count_df[count_df[,'mgi_symbol'] %in% genes(),]
      matfilt = subset(df, select=-c(symbol, mgi_symbol))
    }
    else{
      df = count_df[count_df[,'symbol'] %in% genes(),]
      matfilt = subset(df, select=-c(symbol))
    }
    tcounts <- t(matfilt) %>%
      base::merge(colData, ., by="row.names") %>%
      gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
    if (species == 'mouse') {
      tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol # add gene symbols for goi
    }
    if (species == "rat") {
      tcounts$symbol <- gene_data[tcounts$gene,]$rgd_symbol
    }
    if (species == "human") {
      tcounts$symbol <- gene_data[tcounts$gene,]$hgnc_symbol
    }
    
    # for dotplot
    ifelse(sex=="Both", tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Dataset, Species) %>%
             dplyr::summarise(expression=median(as.double(expression))), tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Sex, Dataset, Species) %>%
             dplyr::summarise(expression=median(as.double(expression))))
    tcounts_med <- tcounts_med[!tcounts_med$Condition %in% paintype, ]
    
    dotplot = ggplot(data = tcounts_med, aes(x=Population, y = symbol)) +
      scale_colour_viridis_c(option = "magma", end = .90) +
      geom_point(aes(col=log(expression), size=log(expression))) + th +
      facet_grid(.~Dataset, scales = "free", space="free") +
      scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
    
    if (sex == "Separate") {
      dotplot = dotplot + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales ="free_x")
    }
    
    output$bulkseq_dots <- renderPlotly({dotplot})
    
    
    # for line plot
    if (is.null(injury_data) == FALSE) {
      count_df = injury_data
      if ('mgi_symbol' %in% colnames(count_df)) {
        df = count_df[count_df[,'mgi_symbol'] %in% genes(),]
        matfilt = subset(df, select=-c(symbol, mgi_symbol))
      }
      else{
        df = count_df[count_df[,'symbol'] %in% genes(),]
        matfilt = subset(df, select=-c(symbol))
      }
      tcounts <- t(matfilt) %>%
        base::merge(colData, ., by="row.names") %>%
        gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
      if (species == 'mouse') {
        tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol # add gene symbols for goi
      }
      if (species == "rat") {
        tcounts$symbol <- gene_data[tcounts$gene,]$rgd_symbol
      }
      if (species == "human") {
        tcounts$symbol <- gene_data[tcounts$gene,]$hgnc_symbol
      }
    }
    
    ifelse(sex=="Both", tcounts_m <- tcounts %>% dplyr::group_by(Condition, symbol, Timepoint, Dataset, Species) %>%
             dplyr::summarise(expression=median(as.double(expression))),
           tcounts_m <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Timepoint, Dataset, Species) %>%
             dplyr::summarise(expression=median(as.double(expression))))
    tcounts_m = tcounts_m[!tcounts_m$Condition %in% "Undetermined", ]
    
    lineplot = ggplot(data=tcounts_m, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) +
      scale_colour_viridis(discrete=TRUE, end = .80) +
      geom_line(aes(color=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + theme_line + ylab("Expression") +
      guides(linetype=guide_legend(label.hjust=0.5))
    if (sex =='Separate'){
      lineplot = lineplot + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales = "free_x") + labs(col="")
    }
    output$bulkseq_lines <- renderPlot({lineplot})
    
    # for subtype, optional
    if (include_subtype == TRUE) {
      ifelse(sex == "Both", tcounts_sub <- tcounts %>% dplyr::group_by(Condition, symbol, Timepoint, Population, Dataset, Species) %>%
               dplyr::summarise(expression=median(as.double(expression))),
             tcounts_sub <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Timepoint, Population, Dataset, Species) %>%
               dplyr::summarise(expression=median(as.double(expression))))
      
      pop_num = n_distinct(tcounts_sub$Population)
      subtypeplot =  ggplot(data=tcounts_sub, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) +
        scale_colour_viridis(discrete=TRUE, end = .80) +
        geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) +
        facet_wrap(~Population, ncol = pop_num, labeller=labeller(Population=population_labels))+
        theme_line + ylab("Expression")
      if (sex == "Separate") {
        subtypeplot = subtypeplot + facet_grid(Sex~Population, labeller=labeller(Population=population_labels, Sex=sexlabels)) + labs(col="")
      }
      output$bulkseq_lines_subtype <- renderPlot({subtypeplot})
    }
    
    # first, process data
    if ('mgi_symbol' %in% colnames(count_df)) {
      goi_df = count_df[count_df[,'mgi_symbol'] %in% genes(),]
    }
    else{
      goi_df = count_df[count_df[,'symbol'] %in% genes(),]
    }
    
    output$goi_table <- DT::renderDataTable({
      datatable = goi_df
      datatable$geneID <- rownames(datatable)
      DT::datatable(
        datatable,
        rownames=datatable$symbol,
        style="default",
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 5)
      )
    })
    
    
    # plot deg plots
    if (include_deg == TRUE) {
      
      geneids = genes()
      if (species == "rat") {
        geneids <- gene_data[gene_data$mgi_symbol %in% genes(),]$rgd_symbol
      }
      if (species == "human") {
        geneids <- gene_data[gene_data$mgi_symbol %in% genes(),]$hgnc_symbol
      }
      
      # generate deg df
      i=1
      final_df = data.frame()
      for (df in deg_df_list) {
        if ("hsymbol" %in% colnames(df)){
          df$symbol = df$hsymbol
        }
        df = df[c("log2FoldChange", "padj", "symbol")]
        colnames(df) = c("LFC","padj", "symbol")
        df$Population = rep(pop_list[i], nrow(df))
        final_df = rbind(final_df, df)
        i=i+1
      }
      final_df = mutate(final_df, sig=ifelse(final_df$padj<0.05, "SIG", "NS"))
      
      deg <- reactive({
        df <- final_df %>% filter(symbol %in% geneids)
        return(df)
      })
      
      # render contrast table
      output$deg_table <- DT::renderDataTable({
        DT::datatable(
          deg(),
          width = 12,
          class = 'nowrap',
          options = list(scrollX = TRUE, pageLength = 10)
        )
      })
      
      degplot = ggplot(data = deg(), aes(x=interaction(Population), y = symbol,
                                         text = paste('padj: ',padj))) +
        scale_colour_viridis_c(option = "viridis", end = .90) +
        geom_point(aes(col=LFC, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
        th + labs(shape = "", size = "") + scale_x_discrete(labels=subpopulation_labels)
      output$deg_plot <- renderPlotly({
        degplot
      })
    }
    
    # download plots and tables in a zip file   
    output$plots <- downloadHandler(
      filename = function() {
        paste("output", "zip", sep=".")
      },
      
      content = function(file){
        # Set temporary working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        fs = c('dotplot.png','lineplot.png')
        
        # Save the plots 
        ggsave('dotplot.png', plot = dotplot, width = 6, height = 6, dpi = 300, units = "in", device='png')
        ggsave('lineplot.png', plot = lineplot, width = 6, height = 6, dpi = 300, units = "in", device='png')
        if (include_subtype == TRUE){
          ggsave("subtype.png", plot=subtypeplot, width = 12, height = 6, dpi = 300, units = "in", device='png')
          fs = c(fs, "subtype.png")
        }
        if(include_deg == TRUE){
          ggsave("degplot.png", plot=degplot, width = 6, height = 6, dpi = 300, units = "in", device='png')
          fs = c(fs, "degplot.png")
        }
        
        # csv files for data 
        write.csv(goi_df,"result_table.csv", row.names = FALSE, col.names = TRUE)
        
        # for deg result table 
        if (include_deg == TRUE){
          write.csv(deg(),"deg_table.csv", row.names = FALSE, col.names = TRUE)
        }
        
        fs = c(fs, "result_table.csv")
        fs = c(fs, "deg_table.csv")
        # Zip them up
        zip(file, fs)
      }
    )
    
  })
  
  if (include_deg == TRUE){
    # for volcano
    generate_dict = function(names, dfs) {
      d = dict()
      for (i in c(1:length(names))) {
        d[[names[i]]] = dfs[i]
      }
      return(d)
    }
    
    pop_df_dic = generate_dict(pop_list, deg_df_list)
    
    # ploting volcano plots
    observeEvent(input$plotvolc, {

      geneids = genes()
      if (species == "rat") {
        geneids <- gene_data[gene_data$mgi_symbol %in% genes(),]$rgd_symbol
      }
      if (species == "human") {
        geneids <- gene_data[gene_data$mgi_symbol %in% genes(),]$hgnc_symbol
      }

      res = as.data.frame(pop_df_dic[[input$volc_pop]])
      res = mutate(res, log10fdr = -log10(padj))
      df = res

      output$volcanoplot <- renderCachedPlot({
        il_genes <- df %>% filter(symbol %in% geneids)
        ggplot(data = df, aes(x = log2FoldChange, y = log10fdr)) +
          geom_point(colour = "grey", alpha = 0.5) +
          geom_point(data = il_genes, # New layer containing data subset il_genes
                     size = 2,
                     shape = 21,
                     fill = "firebrick",
                     colour = "black") +
          geom_hline(yintercept = -log10(0.05),
                     linetype = "dashed") +
          geom_vline(xintercept = c(log2(0.5), log2(2)),
                     linetype = "dashed") +
          geom_label_repel(data = il_genes,
                           aes(label = symbol),
                           force = 2,
                           nudge_y = 1) +
          scale_colour_manual(values = cols) +
          scale_x_continuous(breaks = c(seq(-10, 10, 2)),
                             limits = c(-10, 10)) + theme_line +
          labs(y="-log10(FDR)")}, 
        cacheKeyExpr = {list(df, genes())})
    })
    
    # for contrast table!
    usercontrast <- reactive({
      res = as.data.frame(pop_df_dic[[input$contrast]])
      return(res)
    }) %>% bindCache(input$contrast)
    
    # render contrast table
    output$contrast_table <- DT::renderDataTable({
      DT::datatable(
        usercontrast(),
        width = 12,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })
  }
}



################################################ SERVER ##############################################################################
shinyServer(function(input, output, session) {
  
  pbmc = readRDS(scRNA_dir)
  options(shiny.maxRequestSize=60*1024^2)
  # res_auth <- secure_server(check_credentials = check_credentials(credentials))
  
  # select genes
  updateSelectizeInput(session,
                       inputId = "geneid",
                       label = "Search Genes:",
                       choices = TPM_subtype[,80],
                       server = TRUE,
                       selected = c("Trpv1", "Scn10a","Atf3")
  )
  
  # get a list of genes, separated by comma
  updateTextInput(session,
                  inputId = "genes",
                  label = "Search Genes:"
  )
  
  # linking to datasets
  observeEvent(input$linksubtype, {
    newvalue = "tabsubtype"
    updateTabItems(session, "tabs", newvalue)
  })
  
  ### linking to other tabs
  observeEvent(input$linkmouse, {
    newvalue <- "tabmouse"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkrat, {
    newvalue <- "tabrat"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkscrna, {
    newvalue <- "tabspat"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkipsc, {
    newvalue <- "tabhuman"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkhdrg, {
    newvalue <- "tabhdrg"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkzheng, {
    newvalue <- "tabzheng"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkdb, {
    newvalue <- "tabdb"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$linkcts, {
    newvalue <- "tabcts"
    updateTabItems(session, "tabs", newvalue)
  })
  
  # data preprocessing for deg plots
  mouse_deg_df = generate_degdf(list(b10d2, balb), c("b10d2", "balb"), "Mouse DRG")
  rat_deg_df = generate_degdf(list(ratdrg), c("DRG"), "Rat DRG")
  human_deg_df = generate_degdf(list(young, old), c("young", "old"), "iPSC HSN1")
  db_deg_df = generate_degdf(list(DB, DB_female, DB_male), c("Diabetes", "Diabetes_male","Diabetes_female"), "DPN (Skin)")
  cts_deg_df = generate_degdf(list(HS), c("skin"), "CTS (Skin)")
  
  subpopulations_df = list(TDNV_3D,TDNV_4W,CGRT_3D,CGRT_4W,MRTD_3D, MRTD_4W,CRTH_3D,CRTH_4W,TBAC_3D,TBAC_4W)
  labels = c("Nociceptors 3D","Nociceptors 4W","PEP 3D","PEP 4W","NP 3D","NP 4W","C-LTMR 3D",
             "C-LTMR 4W","Ad- AB-RA LTMRs 3D","Ad- AB-RA LTMRs 4W")
  
  subtype_deg_df = generate_degdf(subpopulations_df, labels, "Subtype DRG (Barry)")
  subtype_deg_df = subtype_deg_df[c("LFC", "padj", "symbol", "Population","Dataset", "sig")]
  
  mouse_all_deg_df = rbind(mouse_deg_df, subtype_deg_df)
  human_all_deg_df = rbind(human_deg_df, db_deg_df, cts_deg_df)
  combined_deg_df = rbind(mouse_all_deg_df, rat_deg_df, human_all_deg_df)
  
  # a reactive variable that records whether a file is uploaded.
  rv <- reactiveValues(
    clear = FALSE,
    data = FALSE
  )
  
  observeEvent(input$file, {
    rv$clear <- FALSE
    rv$data = TRUE
  }, priority = 1000)
  
  observeEvent(input$reset, {
    rv$clear = TRUE
    rv$data = FALSE
  })
  
  # first, process data
  file_input <- reactive({
    if (rv$clear == TRUE) {
      return(NULL)
    }
    if(rv$clear==FALSE && rv$data == TRUE) {
      goi = read.table(input$file$datapath)
      rownames(goi) <- goi[,1]
      goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
      return(goi)}
  })
  
  genes = reactive({
    if (is.null(file_input())) {
      genes = input$geneid
    }
    else {
      genes = file_input()
    }
  })
  
  # plot homepage plots
  observeEvent(input$load, {
    
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
    }
    else {
      genes = file_input()
    }
    
    final_df = generate_combine_dataset(df_list, genes, input$sex, col_list, pain_list, species_list)
    combine_dot = plotcombine_server("dot", final_df, input$sex, genes) #dotplot
    scrna_dot = plothomescdot_server("homespat", pbmc, genes)
    
    degr <- reactive({
      rg = subset(rat_gene_data, mgi_symbol %in% genes)$rgd_symbol
      df <- rat_deg_df %>% filter(symbol %in% rg)
      return(df)
    }) %>% bindCache(genes)
    
    degm <- reactive({
      df <- mouse_all_deg_df %>% filter(symbol %in% genes)
      return(df)
    }) %>% bindCache(genes)
    
    degh <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
      df <- human_all_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(genes)
    
    combine_degplot = deg_combine_server("deg_plot", reactive({degr()}), reactive({degm()}), reactive({degh()}))
    
    # download all plots into a zip 
    output$combineplots <- downloadHandler(
      filename = function() {
        paste("combined", "zip", sep=".")
      },
      
      content = function(file){
        # Set temporary working directory
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        fs = c()
        
        # Save the plots 
        ggsave('combine_deg.png', plot = combine_degplot, width = 14, height = 6, dpi = 300, units = "in", device='png')
        ggsave('combine_dot.png', plot = combine_dot, width = 14, height = 6, dpi = 300, units = "in", device='png')
        ggsave('scrna_dot.png', plot = scrna_dot, width = 10, height = 8, dpi = 300, units = "in", device='png')
        fs = c(fs, 'combine_deg.png')
        fs = c(fs, 'combine_dot.png')
        fs = c(fs, 'scrna_dot.png')
        zip(file, fs)
      }
    )
    
  })
  
  output$combineplot <- downloadHandler(
    
    # For PDF output, change this to "report.pdf"
    filename = "RNAseq_homeplot.html",
    content = function(file) {
      # data pre-processing 
      file_input <- reactive({
        if (rv$clear == TRUE) {
          return(NULL)
        }
        if(rv$clear==FALSE && rv$data == TRUE) {
          goi = read.table(input$file$datapath)
          rownames(goi) <- goi[,1]
          goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
          return(goi)}
      })
      
      if (is.null(file_input())) {
        genes = input$geneid
      }
      else {
        genes = file_input()
      }
      
      final_df = generate_combine_dataset(df_list, genes, input$sex, col_list, pain_list, species_list)
      
      degr <- reactive({
        rg = subset(rat_gene_data, mgi_symbol %in% genes)$rgd_symbol
        df <- rat_deg_df %>% filter(symbol %in% rg)
        return(df)
      }) %>% bindCache(genes)
      
      degm <- reactive({
        df <- mouse_all_deg_df %>% filter(symbol %in% genes)
        return(df)
      }) %>% bindCache(genes)
      
      degh <- reactive({
        hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
        df <- human_all_deg_df %>% filter(symbol %in% hg)
        return(df)
      }) %>% bindCache(genes)
      
      params <- list(matrix = final_df, genes = genes, sex = input$sex, scRNA = pbmc, ratdeg = reactive({degr()}),
                     mousedeg = reactive({degm()}), humandeg = reactive({degh()}), 
                     human_gene_data = human_gene_data)
      
      shiny::withProgress(value = 0,
                          message = 'Rendering plotting report',
                          detail =  'This might take several minutes.',{
                            rmarkdown::render(input = "report/combineplot.Rmd", 
                                              params = params, 
                                              output_file = file,
                                              envir = new.env(parent = globalenv()))
                          })
      
    }
  )
  
  
  
  callModule(shinytab, "subtypetab", input$sex, TPM_subtype, bulkseq_colData,mouse_gene_data,
             "mouse", "ipsi", reactive({genes()}), include_deg = TRUE, list(TDNV_3D, TDNV_4W, CGRT_3D, CGRT_4W,
                                                                                 MRTD_3D, MRTD_4W, CRTH_3D, CRTH_4W, TBAC_3D,
                                                                                 TBAC_4W),
             pop_list = c("Nociceptors 3D","Nociceptors 4W","PEP 3D","PEP 4W","NP 3D",
                          "NP 4W", "C-LTMR 3D","C-LTMR 4W","Ad- AB-RA LTMRs 3D",
                          "Ad- AB-RA LTMRs 4W"), dataset = "Subtype DRG(Barry)", include_subtype = TRUE,injury_data = bulkseq_mat, parent = session)
  
  callModule(shinytab, "mousetab", input$sex, TPM_mouse, TPM_mouse_colData, mouse_gene_data,
             "mouse", "SNI", reactive({genes()}), include_deg = TRUE,
             list(b10d2, balb), pop_list = c("B10D2", "BALB"), dataset = "Bulk DRG (Mouse)",
             include_subtype = FALSE,parent = session)
  
  callModule(shinytab, "rattab", input$sex, TPM_rat, TPM_rat_colData, rat_gene_data,
             "rat", "SNT", reactive({genes()}), include_deg = TRUE,
             list(ratdrg), pop_list = c("rat"), dataset = "Bulk DRG (Rat)",
             include_subtype = FALSE, parent = session)
  
  callModule(shinytab, "ipsctab", input$sex, TPM_ipsc, ipsc_colData, human_gene_data,
             "human", "patient", reactive({genes()}), include_deg = TRUE,
             list(young, old), pop_list = c("iPSCDN_young","iPSCDN_old"), dataset = "iPSC DSN",
             include_subtype = FALSE, parent = session)
  
  callModule(shinytab, "dbtab", input$sex, TPM_HS_diabetes, db_colData, human_gene_data,
             "human", "PDPN", reactive({genes()}), include_deg = TRUE,
             list(DB, DB_female, DB_male), pop_list = c("Diabetes_skin","Diabetes_skin_females",
                                                        "Diabetes_skin_male"), dataset = "Diabetes",
             include_subtype = FALSE, parent = session)
  
  callModule(shinytab, "ctstab", input$sex, TPM_HS_CTS, skin_colData, human_gene_data,
             "human", "pre_Surgery", reactive({genes()}), include_deg = TRUE,
             deg_df_list = list(HS), pop_list = c("Skin HS"), dataset = "Skin CTS",
             include_subtype = FALSE, parent = session)
  callModule(shinyscrna, "spattab", reactive({genes()}), pbmc, parent = session)
  callModule(shinytab, "zhengtab", input$sex, TPM_zheng, zheng_colData,gene_data,"mouse","", reactive({genes()}), parent = session)
  callModule(shinytab, "hdrgtab", input$sex, TPM_humandrg, humandrg_colData,human_gene_data,"human",c("P","Undetermined"), reactive({genes()}),parent = session)
  
  PlotHeight = reactive(
    return(length(data()))
  )
  ### leaflet map for contact details
  output$myMap <- renderLeaflet({
    m <- leaflet() %>% addTiles()
    m <- m %>% setView( -1.238233, 51.756192, zoom = 13)
    m %>% addPopups(-1.2217, 51.76529, "Neural Injury Group")
  })
  
})

options(warn=-1) # remove warnings

