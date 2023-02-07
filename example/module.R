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
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
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
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + ggtitle("Rat (DRG)") +
      labs(shape = "", size = "") + facet_grid(.~Dataset, scale = "free", space='free')  +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) +
      scale_colour_continuous(limits=c(0,1)) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    
    g2 = ggplot(data = df_human, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) + ggtitle("Human") +
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + facet_grid(.~Dataset, scale = "free", space='free') +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    g3 = ggplot(data = df_mouse, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) +
      ggtitle("Mouse (DRG)") + geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) +
      facet_grid(.~Dataset, scale = "free", space='free') + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "") + guides(size=FALSE)
    
    output$deg <- renderCachedPlot({plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))},
                                   cacheKeyExpr = {list(datar, datah, datam)})
    
    output$downloadDeg <- downloadHandler(
      filename = function() {
        paste("deg_plot", ".png", sep = "")
      },
      
      content = function(file) {
        p <- plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
        ggsave(p, filename = file, width = 14, height = 4, dpi = 300, units = "in", device='png')
      })
    
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
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Rat (DRG)") + scale_x_discrete(labels=population_labels)
    g2 = ggplot(df_human, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(col = FALSE, size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Human") + scale_x_discrete(labels=population_labels)
    g3 = ggplot(df_mouse, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Mouse (DRG)") +
      scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
    
    if (sex == "Both"){
      output$dot <- renderCachedPlot({
        plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
      }, cacheKeyExpr = {list(df, genes)})
      
      
    }
    if (sex == "Separate") {
      g1 = g1 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g2 = g2 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g3 = g3 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      output$dot <- renderCachedPlot({
        plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
      }, cacheKeyExpr = {list(df,genes)})
      
    }
    output$downloadDot <- downloadHandler(
      filename = function() {
        paste("plot1", ".png", sep = "")
      },
      
      content = function(file) {
        p <- plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
        ggsave(p, filename = file, width = 14, height = 4, dpi = 300, units = "in", device='png')
      }
    )
    
  })
}

shinyscrna <- function(input, output, session,genes,pbmc) {
  
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
    output$scrna_dots <- renderPlot({
      
      hg = subset(human_gene_data, mgi_symbol %in% genes())$hgnc_symbol
      DotPlot(pbmc, features = hg) +
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size=8, angle = 45, hjust= 1),
              axis.text.y = element_text(size=8),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(), legend.justification = c(0,0.3),
              legend.title = element_text(size=10), legend.key.size = unit(0.4, "cm"))
    })
    
    g = DimPlot(pbmc,label=TRUE)
    output$scrna_umap <- renderCachedPlot({g}, cacheKeyExpr = {list(pbmc)})
    
    output$downloadumap <- downloadHandler(
      filename = function() {
        paste("umap_plot", ".png", sep = "")
      },
      
      content = function(file) {
        ggsave(g, filename = file, width = 10, height = 8, dpi = 300, units = "in", device='png')
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
    g = DotPlot(pbmc, features = hg) + thc + ggtitle("Human DRG subtypes", subtitle="Spatial-seq (Tavares-Ferrelra 2021)")
    output$home_scrna_dots <- renderCachedPlot({g}, cacheKeyExpr = {list(genes, pbmc)})
    
    output$downloadscDot <- downloadHandler(
      filename = function() {
        paste("scplot", ".png", sep = "")
      },
      
      content = function(file) {
        ggsave(g, filename = file, width = 10, height = 8, dpi = 300, units = "in", device='png')
      }
    )
  })
}

shinytab <- function(input, output, session, sex, count_df, colData, gene_data, species, paintype, genes,
                     include_deg=FALSE, deg_df_list=NULL, pop_list=NULL, dataset=NULL, include_subtype=FALSE,
                     injury_data = NULL) {
  observeEvent(input$link_to_home, {
    newvalue <- "tabhome"
    observe({
      updateTabItems(session, "tabs", newvalue)
    })
  })

  
  observeEvent(input$load, {
    # first, process data
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
    
    output$table <- DT::renderDataTable({
      DT::datatable(
        tcounts_med,
        width = 12,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })
    
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
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("drgdirectory_search", ".csv", sep = "")
      },
      
      content = function(file) {
        datatable <- data()
        write.csv(datatable, file, row.names = FALSE)
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
      
      output$deg_plot <- renderPlotly({
        ggplot(data = deg(), aes(x=interaction(Population), y = symbol,
                                 text = paste('padj: ',padj))) +
          scale_colour_viridis_c(option = "viridis", end = .90) +
          geom_point(aes(col=log2FoldChange, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
          th + labs(shape = "", size = "") + scale_x_discrete(labels=subpopulation_labels)
      })
    }
    
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
      res = mutate(res, log10fdr=-log10(padj))
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
          labs(y="-log10(FDR)")}, cacheKeyExpr = {list(df, genes())})
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


