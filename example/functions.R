

## for Database
{
  conn <- dbConnect(RSQLite::SQLite(), "test.db")
  human_gene_data = dbReadTable(conn, "human_gene_data")
  mouse_gene_data = dbReadTable(conn, "mouse_gene_data")
  rat_gene_data = dbReadTable(conn, "rat_gene_data")
}


## for R Data
{
  # load("data/gene.RData")
}

## for network data
load("data/network.RData")

#' This function accepts a list of count matrices, computes median expression value for goi, and combines the results
#' into a meta data frame encompassing all experiments with their pain, species and experimental conditions information
#' included.
#'
#' @param df_list a list of count matrices from different experiments
#' @param genes genes of interest
#' @param sex sex option; 'Both' for plotting both sexes; 'Separate' for dividing datas and graphs based on sex.
#' @param col_list a list of experiment meta data that describes various experimental conditions and treatments
#' @param pain_list a list of strings describing the pain model used in each experiment
#' @param species_list a list of species for each experiment
#' @param var the variable that denotes expression
#' @param ... variables used to group by
#'
#' @return an integrated dataset containing median expression of genes of interest from all experiments
generate_combine_dataset <- function(df_list, genes, sex, col_list, pain_list, species_list, var, ...) {
  final_df = data.frame()
  i = 1
  for (d in names(df_list)) {
    
    # remove unwanted columns
    filt = as.data.frame(df_list[[d]])
    colD = col_list[[d]]
    
   
    matfilt = filt[,!(names(filt) %in% c("symbol", "mgi_symbol"))] # expression data
    colnames(matfilt) = c("gene", colnames(matfilt)[2:length(colnames(matfilt))])
    rownames(matfilt) = matfilt$gene
  
    
    # merge
    tcounts <- t(matfilt) %>%
      base::merge(as.data.frame(colD), ., by="row.names") %>%
      gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
    
    species = species_list[i]
    tcounts$symbol = filt[filt$X %in% tcounts$gene,]$symbol
    
    ifelse(sex=="Both", tcounts_med <- tcounts %>% dplyr::group_by(...) %>%
             dplyr::summarise(expression=median(as.double({{var}}))), tcounts_med <- tcounts %>% dplyr::group_by(..., Sex) %>%
             dplyr::summarise(expression=median(as.double({{var}}))))
    tcounts_med <- tcounts_med[!tcounts_med$Condition %in% pain_list[i], ]
    final_df = rbind(tcounts_med, final_df)
    
    i = i+1
  }
  return(final_df)
}


#' This server generates a combined dot plot for differential gene analysis using
#' data from all experiments. Three independent graphs using data from different
#' species are arranged to allow cross-species visualisation.
#'
#' @param df_rat a DESeq2 result file-like dataset containing 'log2FoldChange', 'padj', 'baseMean', etc.
#' @param df_human a DESeq2 result file-like dataset containing 'log2FoldChange', 'padj', 'baseMean', etc.
#' @param df_mouse a DESeq2 result file-like dataset containing 'log2FoldChange', 'padj', 'baseMean', etc.
#' @param genes query genes
#' @param x_var x axis
#' @param y_var y axis
#' @param color_var variable that determine the color
#' @param shape_var variable that determine the shape
#' @param var the variable to do facet
#'
#' @return None
deg_combine_server <- function(id, df_rat, df_mouse, df_human, x_var,y_var, color_var,shape_var, var) {
  moduleServer(id, function(input, output, session){
    
    x_var = enquo(x_var)
    y_var = enquo(y_var)
    color_var = enquo(color_var)
    shape_var = enquo(shape_var)
    
    # plot for rat
    g1 = ggplot2::ggplot(data = df_rat, aes(x=interaction(!!x_var), y = !!y_var, text = paste('padj: ',padj))) +
      geom_point(aes(col=!!color_var, shape=!!shape_var, size=0.3)) + ggtitle("Rat") +
      labs(shape = "", size = "") + facet_grid(.~.data[[var]], scale = "free", space='free')  +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) +
      scale_colour_continuous(limits=c(0,1)) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    # dot plot for human
    g2 = ggplot2::ggplot(data = df_human, aes(x=interaction(!!x_var), y = !!y_var, text = paste('padj: ',padj))) + ggtitle("Human") +
      geom_point(aes(col=!!color_var, shape=!!shape_var, size=0.3)) + facet_grid(.~.data[[var]], scale = "free", space='free') +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    # dot plot for mouse
    g3 = ggplot2::ggplot(data = df_mouse, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) +
      ggtitle("Mouse") + geom_point(aes(col=!!color_var, shape=!!shape_var, size=0.3)) +
      facet_grid(.~.data[[var]], scale = "free", space='free') + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "") + guides(size=FALSE)
    
    combine_degplot = cowplot::plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
    
    output$deg <- shiny::renderCachedPlot({combine_degplot},
                                          cacheKeyExpr = {list(df_rat, df_mouse, df_human)})
    return(combine_degplot)
    
    
  })
}


#' This server function plots an integrated dot plot showing gene expression levels across
#' different experiments (the DEG plots);
#'
#' @param df an integrated dataset containing median expression of genes of interest from all experiments
#' @param sex sex option; 'Both' for plotting both sexes; 'Separate' for dividing datas and graphs based on sex.
#' @param genes query genes
#' @param x_var Population
#' @param y_var symbol
#' @param col_var expression
#' @param facet_var Dataset # a string
#'
#' @return
plotcombine_server <- function(id, df, sex, genes, x_var, y_var, col_var, facet_var) {
  moduleServer(id, function(input, output, session) {
    
    # a final_df 
    x_var = enquo(x_var)
    y_var = enquo(y_var)
    col_var = enquo(col_var)
    # plot = ggplot(data = df, aes(x = !!x_var, y = !!y_var, fill=!!col_var))  + geom_tile() +
    #   facet_wrap(~Species, scales = "free") + th
    
    df_rat = df[df$Species == "Rat (DRG)",]
    df_mouse = df[df$Species == "Mouse (DRG)",]
    df_human = df[df$Species == "human",]
    
    
    plots = list()
    if (nrow(df_rat) != 0) {
      g1 = ggplot2::ggplot(df_rat, aes(x= !!x_var, y = !!y_var)) + scale_colour_viridis_c(option = "magma", end = .90) +
        scale_fill_continuous(limits=c(-10,10)) + scale_size_continuous(limits=c(-10,10)) +
        geom_point(aes(col=log(!!col_var), size=log(!!col_var))) + thc + guides(col = FALSE, size = FALSE) +
        facet_grid(.~.data[[facet_var]], scale = "free", space='free') + ggtitle("Rat") + scale_x_discrete(labels=population_labels)
      plots = append(plots, list(g1))
      
    }
    
    if (nrow(df_human) != 0) {
      g2 = ggplot2::ggplot(df_human, aes(x= !!x_var, y = !!y_var)) + scale_colour_viridis_c(option = "magma", end = .90) +
        scale_fill_continuous(limits=c(-10,10)) +
        geom_point(aes(col=log(!!col_var), size=log(!!col_var))) + thc + guides(col = FALSE, size = FALSE) +
        facet_grid(.~.data[[facet_var]], scale = "free", space='free') + ggtitle("Human") + scale_x_discrete(labels=population_labels)
      plots = append(plots, list(g2))
      
    }
    
    if (nrow(df_mouse) != 0){
      g3 = ggplot2::ggplot(df_mouse, aes(x= !!x_var, y = !!y_var)) + scale_colour_viridis_c(option = "magma", end = .90) +
        scale_fill_continuous(limits=c(-10,10)) +
        geom_point(aes(col=log(!!col_var), size=log(!!col_var))) + thc + guides(size = FALSE) +
        facet_grid(.~.data[[facet_var]], scale = "free", space='free') + ggtitle("Mouse") +
        scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
      plots = append(plots, list(g3))
    }
    
    
    # function to adjust rel_width! 
    n = length(plots)
    if (n==1){
      plot = cowplot::plot_grid(plots[[1]],ncol=1)
    }
    if(n==2){
      plot = cowplot::plot_grid(plots[[1]], plots[[2]],ncol=2, rel_widths = c(4/11,7/11))
    }
    if (n == 3){
      plot = cowplot::plot_grid(plots[[1]], plots[[2]], plots[[3]],ncol=3, rel_widths = c(1/11,4/11,6/11))
    }
    
    if (sex == "Separate") {
      stopifnot("Sex" %in% colnames(df))
      g1 = g1 + facet_grid(.~.data[[facet_var]]+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g2 = g2 + facet_grid(.~.data[[facet_var]]+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g3 = g3 + facet_grid(.~.data[[facet_var]]+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      plot = plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
    }
    
    output$dot <- shiny::renderCachedPlot({plot}, cacheKeyExpr = {list(df,genes)})
    
    return(plot)
    
  })
}


#' This server function generates snRNA dotplots and dimplot using a Seurat object
#' 
#' @param pbmc a Seurat object 
#' @param parent Default = NULL; a variable used to inherit environment from parent namespace
#' @param genes query genes 
#' 
#' @return None 
shinyscrna <- function(input, output, session,genes,pbmc, parent=NULL) {
  
  if (is.null(parent) == FALSE) {
    shiny::observeEvent(input$link_to_home, {
      newvalue <- "tabhome"
      shiny::observe({
        shinydashboard::updateTabItems(session=parent, "tabs", newvalue)
      })
    })
  }
  
  # update the search bar
  observe({
    shiny::updateSelectizeInput(session,
                                inputId = "scgeneid",
                                label = NULL,
                                choices = rownames(pbmc),
                                server = TRUE,
                                selected = "ATF3"
    )
  })
  
  observeEvent(input$load, {
    
    ## for database 
    {
      sql = "SELECT * FROM human_gene_data WHERE mgi_symbol = ?"
      hg = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol = genes())))$hgnc_symbol
    }
    
    
    ## for R Data
    {
      # hg = subset(human_gene_data, mgi_symbol %in% genes())$hgnc_symbol
    }
    
    dotplot = Seurat::DotPlot(pbmc, features = hg) +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(size=8, angle = 45, hjust= 1),
            axis.text.y = element_text(size=8),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(), legend.justification = c(0,0.3),
            legend.title = element_text(size=10), legend.key.size = unit(0.4, "cm"))
    output$scrna_dots <- shiny::renderPlot({
      dotplot
    })
    
    dimplot = Seurat::DimPlot(pbmc,label=TRUE)
    output$scrna_umap <- shiny::renderCachedPlot({dimplot}, cacheKeyExpr = {list(pbmc)})
    
    # download plots and tables in a zip file   
    output$scrna_plots <- shiny::downloadHandler(
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
  
  shiny::observeEvent(input$scgeneid,{
    output$scrna_feature <- plotly::renderPlotly({
      FeaturePlot(pbmc, features = input$scgeneid)
    })
  })
  
}

#' This server function generates an integrated dotplot containing all scRNA data in the homepage 
#' 
#' @param pbmc a Seurat object 
#' @param genes query genes 
#' 
#' @return None 
plothomescdot_server <- function(id, pbmc, genes) {
  moduleServer(id, function(input, output, session) {
    ## For Database
    {
      sql = "SELECT * FROM human_gene_data WHERE mgi_symbol = ?"
      hg = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol = genes)))$hgnc_symbol
    }

    
    ## for RData
    {
      # hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
    }
    
    g = Seurat::DotPlot(pbmc, features = hg) + thc + ggtitle("Human DRG subtypes", subtitle="Spatial-seq (Tavares-Ferreira 2021)")
    output$home_scrna_dots <- renderCachedPlot({g}, cacheKeyExpr = {list(genes, pbmc)})
    return(g)
    
  })
}

#' @param sex sex option; 'Both' for plotting both sexes; 'Separate' for dividing datas and graphs based on sex.
#' @param count_df count matrix for each gene
#' @param colData experimental conditions and treatement of each dataset
#' @param gene_data mapping of gene symbols to the specified species
#' @param species the species of the input dataset
#' @param paintype the paintype of input experiment
#' @param genes genes searched in the search bar
#' @param include_deg default = FALSE; set to TRUE when the dataset has result from deseq2
#' @param deg_df_list not NULL when include_deg == TRUE; a list of DESeq2 result file-like dataframes
#' @param pop_list default = NULL; set to not NULL when include_deg == TRUE; a list of population names in the experiment
#' @param dataset default = NULL; set to not NULL when include_deg == TRUE;
#' @param include_subtype default = FALSE; set to TRUE when the dataset has different sub-populations
#' @param injury_data default = NULL; a count matrix containing both injury and naive RNAseqs.
#' @param parent default = NULL; a parameter that passes the global environment
shinypage <- function(input, output, session, sex, count_df, colData, species, paintype, genes,
                      include_deg=FALSE, deg_df_name = NULL, dataset=NULL, include_subtype=FALSE,
                      injury_data = NULL,include_count = TRUE, parent=NULL) {
  
  # generate R markdown reports
  output$report <- shiny::downloadHandler(
    
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
      
      count_matrix = NULL
      if (include_count == TRUE) {
        count_matrix = dbReadTable(conn, count_df)
      }
      
      
      injury_matrix = NULL
      
      if (is.null(injury_data) == FALSE) {
        injury_matrix = dbReadTable(conn, injury_data)
      }
      
      deg_df = NULL
      
      # available for deg analysis 
      if (include_deg == TRUE){
        deg_df = dbReadTable(conn, deg_df_name)
        if (!"Population" %in% names(deg_df)){
          deg_df$Population = rep("", nrow(deg_df))
        }
      }
  
      params <- list(matrix = count_matrix, colD = colData, species = species, paintype = paintype, genes = genes(), sex= sex, injury_data = injury_matrix,
                     include_deg = include_deg, deg_df = deg_df, gene_data = gene_data, include_count = include_count)
      
      shiny::withProgress(value = 0,
                          message = 'Rendering plotting report',
                          detail =  'This might take several minutes.',{
                            rmarkdown::render(input = "reports/RNAseq_plot.Rmd",
                                              params = params,
                                              output_file = file,
                                              envir = new.env(parent = globalenv()))
                          })
      
    }
  )
  
  shiny::observeEvent(input$load, {
    
    if (include_count == TRUE){
      # query count datasets
      
      # for database 
      {
        if (species == "human") {
          sql = paste("SELECT * FROM", count_df, "WHERE mgi_symbol = ?")
        }
        else{
          sql = paste("SELECT * FROM", count_df, "WHERE symbol = ?")
        }
        filt = RSQLite::dbGetPreparedQuery(conn, sql, bind.data=data.frame(symbol=genes()))
        filt = as.data.frame(filt)
      }
     
      
      
      # using R Data
      {
        # if (species == "human"){
        #   filt = count_df[count_df$mgi_symbol %in% genes(),]
        #   X = rownames(filt)
        #   filt = cbind(X, filt)
        #   
        # }
        # else {
        #   filt = count_df[count_df$symbol %in% genes(),]
        #   X = rownames(filt)
        #   filt = cbind(X, filt)
        # }
      }
      
      
      matfilt = filt[,!(names(filt) %in% c("symbol", "mgi_symbol"))] # expression data
      colnames(matfilt) = c("gene", colnames(matfilt)[2:length(colnames(matfilt))])
      rownames(matfilt) = matfilt$gene

      tcounts <- t(matfilt) %>%
        base::merge(colData, ., by="row.names") %>%
        tidyr::gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
      
      tcounts$symbol = filt[filt$X %in% tcounts$gene,]$symbol
      
      # get median expression
      tcounts_med = get_median(tcounts, sex, expression, Condition, Population, symbol, Dataset, Species)
      tcounts_med <- tcounts_med[!tcounts_med$Condition %in% paintype, ] # remove pain samples
      
      dotplot = plot_dotplot(tcounts_med, Population, symbol, expression, "Dataset", sex)
      output$bulkseq_dots <- plotly::renderPlotly({dotplot})
      
      
      # plot line plot
      if (is.null(injury_data) == FALSE) {
        # # for database
        {
          count_df = injury_data
          sql = paste("SELECT * FROM", count_df, "WHERE symbol = ?")
          filt = RSQLite::dbGetPreparedQuery(conn, sql, bind.data=data.frame(symbol=genes()))
          filt = as.data.frame(filt)
        }
        
        
        ## for R Data
        {
          # filt = injury_data[injury_data$symbol %in% genes(),]
        }
        
        rownames(filt) = filt$gene
        matfilt = filt[,!(names(filt) %in% c("symbol", "mgi_symbol", "gene", "description"))]
        tcounts <- t(matfilt) %>%
          base::merge(colData, ., by="row.names") %>%
          tidyr::gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
        
        tcounts$symbol = filt[filt$gene %in% tcounts$gene,]$symbol
      } # if provided injury data
      
      
      # data for lineplot
      tcounts_m = get_median(tcounts, sex, var = expression, Condition, symbol, Timepoint, Dataset, Species)
      tcounts_m = tcounts_m[!tcounts_m$Condition %in% "Undetermined", ]
      lineplot = plot_lineplot(tcounts_m, sex, Condition, expression, symbol, Timepoint, symbol, Timepoint)
      output$bulkseq_lines <- shiny::renderPlot({lineplot})
      
      # plotting line graphs for each population
      if (include_subtype == TRUE) {
        tcounts_sub = get_median(tcounts, sex, expression, Condition, symbol, Timepoint, Population, Dataset, Species)
        pop_num = n_distinct(tcounts_sub$Population) # number of columns
        subtypeplot = plot_subtype(tcounts_sub, sex, pop_num, Condition, expression, symbol, "Population", Timepoint, symbol, Timepoint)
        output$bulkseq_lines_subtype <- shiny::renderPlot({subtypeplot})
      }
      
      # subset the count data for genes of interest
      output$goi_table <- DT::renderDataTable({
        datatable = filt
        datatable$geneID <- rownames(datatable)
        DT::datatable(
          datatable,
          rownames=datatable$symbol,
          style="default",
          class = 'nowrap',
          options = list(scrollX = TRUE, pageLength = 5)
        )
      })
    }
    
    
    
    # if deg_df data is included, plot deg plots
    if (include_deg == TRUE) {
      
      geneids = genes()
      # convert geneids to species-specific
    
      ## for Database 
      {
        if (species == "rat") {
          sql = "SELECT * FROM rat_gene_data WHERE mgi_symbol = ?"
          geneids = as.data.frame(RSQLite::dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol = genes())))$rgd_symbol
        }
        if (species == "human") {
          sql = "SELECT * FROM human_gene_data WHERE mgi_symbol = ?"
          geneids = as.data.frame(RSQLite::dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol = genes())))$hgnc_symbol
        }

        sql = paste("SELECT * from", deg_df_name, "WHERE symbol = ?")
        final_df = as.data.frame(RSQLite::dbGetPreparedQuery(conn, sql, bind.data = data.frame(symbol = geneids)))
        final_df = mutate(final_df, sig=ifelse(final_df$padj<0.05, "SIG", "NS"))
      }
      
      
      ## for Rdata
      {
        # if (species == "rat") {
        #   geneids = rat_gene_data[rat_gene_data$mgi_symbol %in% genes(),]$rgd_symbol
        #   
        # }
        # if (species == "human"){
        #   geneids = human_gene_data[human_gene_data$mgi_symbol %in% genes(),]$hgnc_symbol
        # }
        # 
        # final_df = deg_df_name[deg_df_name$symbol %in% geneids,]
        # final_df = mutate(final_df, sig=ifelse(final_df$padj<0.05, "SIG", "NS"))
      }

      
      if (!"Population" %in% names(final_df)){
        final_df$Population = rep(dataset, nrow(final_df))
      }
    
      
      if (nrow(final_df) != 0) {
        degplot = plot_degplot(final_df, Population, symbol, log2FoldChange, sig)
        output$deg_plot <- plotly::renderPlotly({degplot})
      }
      
      
      if (nrow(final_df) == 0) {
        output$textone <- renderUI({
          p(strong("Plots unavailable. Selected genes are absent from this dataset. Please select other gene symbols."))
        })
      }
      
      # render contrast table
      output$deg_table <- DT::renderDataTable({
        DT::datatable(
          final_df,
          width = 12,
          class = 'nowrap',
          options = list(scrollX = TRUE, pageLength = 10)
        )
      })
      
      ## for database 
      {
        combine_df = dbReadTable(conn, deg_df_name)
        combine_df = mutate(combine_df, sig=ifelse(combine_df$padj<0.05, "SIG", "NS"))
      }
     
      
      ## for RData
      {
        # combine_df = deg_df_name
        # combine_df = mutate(combine_df, sig=ifelse(combine_df$padj<0.05, "SIG", "NS"))
      }
      
      
      if (!"Population" %in% names(combine_df)){
        combine_df$Population = rep(dataset, nrow(combine_df))
      }
      
      shiny::updateSelectizeInput(session,
                                  inputId = "volc_pop",
                                  label = NULL,
                                  choices = unique(combine_df$Population),
                                  server = TRUE,
                                  selected = NULL
      )
      
      shiny::observeEvent(input$plotvolc, {
        
        geneids = genes()
        # convert geneids to species-specific
        
        ## for Database
        {
          if (species == "rat") {
            sql = "SELECT * FROM rat_gene_data WHERE mgi_symbol = ?"
            geneids = as.data.frame(RSQLite::dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol = genes())))$rgd_symbol
          }
          if (species == "human") {
            sql = "SELECT * FROM human_gene_data WHERE mgi_symbol = ?"
            geneids = as.data.frame(RSQLite::dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol = genes())))$hgnc_symbol
          }
        }
       
        
        ## for RData
        {
          # if (species == "rat") {
          #   geneids = rat_gene_data[rat_gene_data$mgi_symbol %in% genes(),]$rgd_symbol
          #   
          # }
          # if (species == "human"){
          #   geneids = human_gene_data[human_gene_data$mgi_symbol %in% genes(),]$hgnc_symbol
          # }
        }
       
        
        # get the df for the selected population
        res = combine_df[combine_df$Population == input$volc_pop,]
        res = mutate(res, log10fdr=-log10(padj))
        il_genes <- res %>% filter(symbol %in% geneids)
        volcanoplot = plot_volcanoplot(res, il_genes, log2FoldChange, log10fdr, symbol)
        output$volcanoplot <- renderCachedPlot({
          volcanoplot}, cacheKeyExpr = {list(res, genes())})
      })
      
      shiny::updateSelectizeInput(session,
                                  inputId = "contrast",
                                  label = NULL,
                                  choices = unique(combine_df$Population),
                                  server = TRUE,
                                  selected = NULL
      )
      
      # for contrast table!
      usercontrast <- reactive({
        res = combine_df[combine_df$Population == input$contrast,]
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
    
    # download plots and tables in a zip file
    output$plots <- shiny::downloadHandler(
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
        
        # csv files for gene count data
        write.csv(count_df,"result_table.csv", row.names = FALSE, col.names = TRUE)
        
        fs = c(fs, "result_table.csv")
        
        # for deg result table
        if (include_deg == TRUE){
          write.csv(deg(),"deg_table.csv", row.names = FALSE, col.names = TRUE)
        }
        
        
        fs = c(fs, "deg_table.csv")
        # Zip them up
        zip(file, fs)
      }
    )
    
  })
  
}

#' # add network functions 
#' #' This function retrieves the information about a list of proteins
#' #'
#' #' @param proteins
#' #'
#' #' @return
#' query_information = function(proteins){
#'   resultTable <- biomaRt::getBM(attributes = c("hgnc_symbol", "description"),
#'                                 filters    = "hgnc_symbol",
#'                                 values     = proteins,
#'                                 mart       = human) # pass mart object to network.Rdata
#'   
#'   Description = c()
#'   for (i in c(1:nrow(resultTable))){
#'     df = resultTable[i,]
#'     protein = df$hgnc_symbol
#'     des = strsplit(df$description, split = "\\[")[[1]][1]
#'     Description = c(Description, as.character(des))
#'   }
#'   result = resultTable
#'   result$Description = Description
#'   return(result[,c("hgnc_symbol", "Description")])
#' }

# add network functions 
#' This function retrieves the information about a list of proteins
#'
#' @param proteins
#'
#' @return
query_information = function(proteins){
  resultTable <- biomaRt::getBM(attributes = c("hgnc_symbol", "description"),
                                filters    = "hgnc_symbol",
                                values     = proteins,
                                mart       = human) # pass mart object to network.Rdata
  
  Description = c()
  for (i in c(1:nrow(resultTable))){
    df = resultTable[i,]
    protein = df$hgnc_symbol
    des = strsplit(df$description, split = "\\[")[[1]][1]
    Description = c(Description, as.character(des))
  }
  result = resultTable
  result$Description = Description
  return(result[,c("hgnc_symbol", "Description")])
}


#' This function plots volcano and deg plots based on pop selected
#'
#' @param interactions interaction data obtained from STRING database
#' @param gene_list a list of query genes
#' @param deg_df a data frame depending on the mode; for mode == "meta", a data frame containing
#' Pain Enrichment Score for all genes are given; for mode == "indiv", a data frame containing
#' differential gene expression information is given.
#' @param mode the mode of network. if mode == "meta", then the color will be generated based on
#' their Pain Enrichment Score
#'
#' @return
get_nodes = function(interactions, gene_list, deg_df, mode, snps, pg1, pg2){
  deg_df = deg_df[!is.na(deg_df$symbol),]
  proteins = unique(c(unique(interactions$preferredName_A), unique(interactions$preferredName_B)))
  nodes = data.frame()
  id = unique(c(unique(interactions$stringId_A), unique(interactions$stringId_B)))
  nodes= as.data.frame(id)
  nodes$label = proteins
  
  if (mode == "indiv"){
    deg_df = deg_df[!is.na(deg_df$log2FoldChange),]
    # deg_df = mutate(deg_df, log2FoldChange = ifelse(deg_df$padj < 0.05, deg_df$log2FoldChange, 0))
    deg_df = deg_df %>% arrange(log2FoldChange)
    # deg_df$color.background = c(viridis::viridis(nrow(deg_df), begin = 0.6, end = 1.0))
    
    high = deg_df[deg_df$log2FoldChange > 0,]
    low = deg_df[deg_df$log2FoldChange < 0,]
    zero = deg_df[deg_df$log2FoldChange == 0,]
    
    high = high %>% arrange(log2FoldChange)
    colfunc <- colorRampPalette(c("#FBFFFF","#961252"))
    high$color = c(colfunc(nrow(high)))
    
    colfunc <- colorRampPalette(c("#0277bd","#FBFFFF"))
    low = low %>% arrange(log2FoldChange)
    low$color = c(colfunc(nrow(low)))
    zero$color = rep("#FBFFFF", nrow(zero))
    
    deg_df = rbind(high, low, zero)
    
    colors = c()
    font_colors = c()
    q1 = high[(3 * nrow(high))/4,]$log2FoldChange
    q3 = low[(nrow(low))/4,]$log2FoldChange
    
    print(q1)
    print(q3)
    for (protein in nodes$label){
      fc = "black"
      if (protein %in% deg_df$symbol) {
        df = deg_df[deg_df$symbol == protein,]
        if (df$padj < 0.05){
          colors = c(colors, df$color)
          if ((df$log2FoldChange > q1) | (df$log2FoldChange < q3)) {
            fc = "white"
          }
        }
        else {
          colors = c(colors, "white")
        }
      }
      if (!protein %in% deg_df$symbol)  {
        colors = c(colors, "lightgrey")
      }
      font_colors = c(font_colors, fc)
    }
    
    nodes$color.background = colors
    nodes$font.color = font_colors
    
  }
  
  # pain enrichment score
  if (mode == "pes") {
    deg_df = deg_df %>% arrange(score)
    high = deg_df[deg_df$score > 0.60,]
    low = deg_df[deg_df$score <= 0.60,]
    
    colfunc <- colorRampPalette(c("#FBFFFF","#961252"))
    high = high %>% arrange(score)
    high$color = c(colfunc(nrow(high)))
    
    colfunc <- colorRampPalette(c("#0277bd","#FBFFFF"))
    low = low %>% arrange(score)
    low$color = c(colfunc(nrow(low)))
    
    deg_df = rbind(high, low)
    # colfunc <- colorRampPalette(c("#4E0707", "#104a8e"))
    # deg_df = rbind(deg_df, data.frame(symbol = "NA", score = -10.7, ID = "NA"))
    # deg_df$color = c(colfunc(nrow(deg_df)))
    
    colors = c() 
    font_colors = c()
    q1 = high[(3 * nrow(high))/4,]$score
    q3 = low[(nrow(low))/4,]$score
    
    for (protein in nodes$label){
      fc = "black"
      if (protein %in% deg_df$symbol) {
        c = deg_df[deg_df$symbol == protein,]$color
        colors = c(colors,c)
        
        if ((deg_df[deg_df$symbol == protein,]$score > q1) | (deg_df[deg_df$symbol == protein,]$score < q3)) {
          fc = "white"
        }
      }
      else {
        colors = c(colors, "lightgrey")
      }
      font_colors = c(font_colors, fc)
    }
    
    nodes$color.background = colors
    nodes$font.color = font_colors
    
  }
  
  # get the global pain score, scale based on its value
  if (mode == "meta") {
    deg_df = deg_df %>% arrange(score)
    
    high = deg_df[deg_df$score > 0,]
    low = deg_df[deg_df$score < 0,]
    zero = deg_df[deg_df$score == 0,]
    
    colfunc <- colorRampPalette(c("#FBFFFF","#961252"))
    high = high %>% arrange(score)
    high$color = c(colfunc(nrow(high)))
    
    colfunc <- colorRampPalette(c("#0277bd","#FBFFFF"))
    low = low %>% arrange(score)
    low$color = c(colfunc(nrow(low)))
    zero$color = "#FBFFFF"
    
    deg_df = rbind(high, low, zero)
    
    # colfunc <- colorRampPalette(c("#4E0707", "#104a8e"))
    # deg_df = rbind(deg_df, data.frame(symbol = "NA", score = -10.7, ID = "NA"))
    # deg_df$color = c(colfunc(nrow(deg_df)))
    
    colors = c() 
    font_colors = c()
    q1 = high[(3 * nrow(high))/4,]$score
    q3 = low[(nrow(low))/4,]$score
    for (protein in nodes$label){
      fc = "black"
      if (protein %in% deg_df$symbol) {
        c = deg_df[deg_df$symbol == protein,]$color
        colors = c(colors,c)
        if ((deg_df[deg_df$symbol == protein,]$score > q1) | (deg_df[deg_df$symbol == protein,]$score < q3)) {
          fc = "white"
        }
      }
      else {
        colors = c(colors, "lightgrey")
      }
      font_colors = c(font_colors, fc)
    }
    
    nodes$color.background = colors
    nodes$font.color = font_colors
    
  }
  
  nodes = mutate(nodes, font.strokeColor = "#343434") # border colour for text
  nodes = mutate(nodes, font.strokeWidth = ifelse(nodes$label %in% gene_list, 1, 0)) # bold for selected genes
  nodes = mutate(nodes, borderWidth = ifelse(nodes$label %in% snps$symbol, 1, 0.25)) #node border for HPGD (snps)
  nodes = mutate(nodes, color.border = ifelse(nodes$label %in% snps$symbol, "black", "#444444")) # red border for HPGD (snps)
  nodes = mutate(nodes, shape = ifelse(nodes$label %in% pg1$symbol, "box", "circle")) # shape for genes in pain gene database (pg1)
  nodes = mutate(nodes, shadow = ifelse(nodes$label %in% pg2$symbol, TRUE, FALSE)) # shadow for genes in dolorisk priority list
  return(nodes)
  
  # # nodes = mutate(nodes, font.strokeColor = "#FFFFFF") # border colour for text
  # # nodes = mutate(nodes, font.strokeWidth = 2) # border colour for text
  # nodes = mutate(nodes, font.background = ifelse(nodes$label %in% pg2$symbol, "#597D35", "undefined")) 
  # # nodes = mutate(nodes, font.strokeWidth = ifelse(nodes$label %in% gene_list, 2, 1)) # bold for selected genes
  # nodes = mutate(nodes, borderWidth = ifelse(nodes$label %in% snps$symbol, 1.5, 0.25)) #node border for HPGD (snps)
  # nodes = mutate(nodes, color.border = ifelse(nodes$label %in% snps$symbol, "black", "#444444")) # red border for HPGD (snps)
  # nodes = mutate(nodes, shape = ifelse(nodes$label %in% pg1$symbol, "box", "circle")) # shape for genes in pain gene database (pg1)
  # nodes = mutate(nodes, shadow = ifelse(nodes$label %in% gene_list, TRUE, FALSE)) # shadow for genes in dolorisk priority list
  # return(nodes)
}



#' #' This function plots volcano and deg plots based on pop selected
#' #'
#' #' @param interactions interaction data obtained from STRING database
#' #' @param gene_list a list of query genes
#' #' @param deg_df a data frame depending on the mode; for mode == "meta", a data frame containing
#' #' Pain Enrichment Score for all genes are given; for mode == "indiv", a data frame containing
#' #' differential gene expression information is given.
#' #' @param mode the mode of network. if mode == "meta", then the color will be generated based on
#' #' their Pain Enrichment Score
#' #'
#' #' @return
#' get_nodes = function(interactions, gene_list, deg_df, mode, snps, pg1, pg2){
#'   deg_df = deg_df[!is.na(deg_df$symbol),]
#'   pain_score = pain_score[!is.na(pain_score$symbol),]
#'   proteins = unique(c(unique(interactions$preferredName_A), unique(interactions$preferredName_B)))
#'   nodes = data.frame()
#'   id = unique(c(unique(interactions$stringId_A), unique(interactions$stringId_B)))
#'   nodes= as.data.frame(id)
#'   nodes$label = proteins
#'   
#'   if (mode == "comp") {
#'     score_df = score_df %>% arrange(score)
#'     
#'     high = score_df[score_df$score > median(score_df$score),]
#'     low  = score_df[score_df$score < median(score_df$score),]
#'     zero = score_df[score_df$score == median(score_df$score),]
#'     
#'     colfunc <- colorRampPalette(c("#FBFFFF","#961252"))
#'     high = high %>% arrange(score)
#'     high$color = c(colfunc(nrow(high)))
#'     
#'     colfunc <- colorRampPalette(c("#0277bd","#FBFFFF"))
#'     low = low %>% arrange(score)
#'     low$color = c(colfunc(nrow(low)))
#'     zero$color = "#FBFFFF"
#'     
#'     score_df = rbind(high, low, zero)
#'     
#'     # colfunc <- colorRampPalette(c("#4E0707", "#104a8e"))
#'     # deg_df = rbind(deg_df, data.frame(symbol = "NA", score = -10.7, ID = "NA"))
#'     # deg_df$color = c(colfunc(nrow(deg_df)))
#'     
#'     colors = c() 
#'     for (protein in nodes$label){
#'       if (protein %in% score_df$symbol) {
#'         c = score_df[score_df$symbol == protein,]$color
#'         colors = c(colors,c)
#'       }
#'       else {
#'         colors = c(colors, "lightgrey")
#'       }
#'     }
#'     
#'     nodes$color.background = colors
#'     
#'   }
#'   
#'   if (mode == "indiv"){
#'     # deg_df = deg_df[!is.na(deg_df$symbol),]
#'     deg_df = deg_df[!is.na(deg_df$log2FoldChange),]
#'     # deg_df = mutate(deg_df, log2FoldChange = ifelse(deg_df$padj < 0.05, deg_df$log2FoldChange, 0))
#'     deg_df = deg_df %>% arrange(log2FoldChange)
#'     # deg_df$color.background = c(viridis::viridis(nrow(deg_df), begin = 0.6, end = 1.0))
#'     
#'     high = deg_df[deg_df$log2FoldChange > 0,]
#'     low = deg_df[deg_df$log2FoldChange < 0,]
#'     zero = deg_df[deg_df$log2FoldChange == 0,]
#'     
#'     high = high %>% arrange(log2FoldChange)
#'     colfunc <- colorRampPalette(c("#FBFFFF","#961252"))
#'     high$color = c(colfunc(nrow(high)))
#'     
#'     colfunc <- colorRampPalette(c("#0277bd","#FBFFFF"))
#'     low = low %>% arrange(log2FoldChange)
#'     low$color = c(colfunc(nrow(low)))
#'     zero$color = rep("#FBFFFF", nrow(zero))
#'       
#'     deg_df = rbind(high, low, zero)
#'     
#'     colors = c()
#'     for (protein in nodes$label){
#'       if (protein %in% deg_df$symbol) {
#'         df = deg_df[deg_df$symbol == protein,]
#'         if (df$padj < 0.05){
#'           colors = c(colors, df$color)
#'         }
#'         else {
#'           colors = c(colors, "white")
#'         }
#'       }
#'       if (!protein %in% deg_df$symbol)  {
#'         colors = c(colors, "lightgrey")
#'       }
#'     }
#'     
#'     nodes$color.background = colors
#'     
#'    
#'   }
#'   
#'   # get the global pain score, scale based on its value
#'   if (mode == "meta") {
#'     pain_score = pain_score %>% arrange(score)
#'    
#'     high = pain_score[pain_score$score > median(pain_score$score),]
#'     low = pain_score[pain_score$score < median(pain_score$score),]
#'     zero = pain_score[pain_score$score == median(pain_score$score),]
#'     
#'     colfunc <- colorRampPalette(c("#FBFFFF","#961252"))
#'     high = high %>% arrange(score)
#'     high$color = c(colfunc(nrow(high)))
#'     
#'     colfunc <- colorRampPalette(c("#0277bd","#FBFFFF"))
#'     low = low %>% arrange(score)
#'     low$color = c(colfunc(nrow(low)))
#'     zero$color = "#FBFFFF"
#'     
#'     pain_score = rbind(high, low, zero)
#'     
#'     colors = c() 
#'     for (protein in nodes$label){
#'       if (protein %in% pain_score$symbol) {
#'         c = pain_score[pain_score$symbol == protein,]$color
#'         colors = c(colors,c)
#'       }
#'       else {
#'         colors = c(colors, "lightgrey")
#'       }
#'     }
#'     
#'     nodes$color.background = colors
#'     
#'   }
  
  
  # nodes = mutate(nodes, font.strokeColor = "#343434") # border colour for text 
  # nodes = mutate(nodes, font.strokeWidth = ifelse(nodes$label %in% gene_list, 1, 0)) # bold for selected genes
  # nodes = mutate(nodes, borderWidth = ifelse(nodes$label %in% snps$symbol, 1, 0.25)) #node border for HPGD (snps)
  # nodes = mutate(nodes, color.border = ifelse(nodes$label %in% snps$symbol, "black", "#444444")) # red border for HPGD (snps)
  # nodes = mutate(nodes, shape = ifelse(nodes$label %in% pg1$symbol, "box", "circle")) # shape for genes in pain gene database (pg1)
  # nodes = mutate(nodes, shadow = ifelse(nodes$label %in% pg2$symbol, TRUE, FALSE)) # shadow for genes in dolorisk priority list
  # return(nodes)
# }

# Function to query gene interactions from String API
query_interactions <- function(symbols, add_color_nodes) {
  url <- paste0("https://string-db.org/api/tsv/interactions?identifiers=", symbols, "&add_color_nodes=",add_color_nodes)
  resp <- httr::GET(url)
  tsv <- httr::content(resp, "text")
  data <- read.delim(text = tsv, sep = "\t")
  return(data)
}

# plot functions 

#' a function that plot dotplots
#'
#' @param data
#' @param x_var
#' @param y_var
#' @param col_var
#' @param facet_var
#' @param sex
#'
#' @return a dotplot
plot_dotplot <- function(data, x_var, y_var, col_var, facet_var, sex){
  x_var = enquo(x_var)
  y_var = enquo(y_var)
  col_var = enquo(col_var)
  
  g = ggplot2::ggplot(data = data, aes(x=!!x_var, y = !!y_var)) +
    scale_colour_viridis_c(option = "magma", end = .90) +
    geom_point(aes(col=log(!!col_var), size=log(!!col_var))) + th +
    facet_grid(.~.data[[facet_var]], scales = "free", space="free") +
    scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
  
  if (sex == "Separate") {
    g = g + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales ="free_x")
  }
  return(g)
}

#' @param group_vars symbol, Timepoint
#' @param linetype_var Timepoint
#'
plot_lineplot <- function(data,sex, x_var, y_var, col_var, linetype_var, ...){
  x_var = enquo(x_var)
  y_var = enquo(y_var)
  col_var = enquo(col_var)
  group_vars = enquos(...)
  linetype_var = enquo(linetype_var)
  
  lineplot = ggplot2::ggplot(data, aes(x=!!x_var, y=!!y_var, group=interaction(!!!group_vars))) +
    scale_colour_viridis(discrete=TRUE, end = .80) +
    geom_line(aes(color=!!col_var, linetype=!!linetype_var)) + geom_point(aes(col=!!col_var)) + theme_line + ylab("Expression") +
    guides(linetype=guide_legend(label.hjust=0.5))
  if (sex =='Separate'){
    lineplot = lineplot + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales = "free_x") + labs(col="")
  }
  
  return(lineplot)
}

plot_subtype <- function(data, sex, n_col, x_var, y_var, col_var, facet_var, linetype_var, ...){
  x_var = enquo(x_var)
  y_var = enquo(y_var)
  col_var = enquo(col_var)
  group_vars = enquos(...)
  linetype_var = enquo(linetype_var)
  
  subtypeplot =  ggplot2::ggplot(data=data, aes(x=!!x_var, y=!!y_var, group=interaction(!!!group_vars))) +
    scale_colour_viridis(discrete=TRUE, end = .80) +
    geom_line(aes(col=!!col_var, linetype=!!linetype_var)) + geom_point(aes(col=!!col_var)) +
    facet_wrap(~.data[[facet_var]], ncol = n_col, labeller=labeller(Population=population_labels))+
    theme_line + ylab("Expression")
  if (sex == "Separate") {
    subtypeplot = subtypeplot + facet_grid(Sex~.data[[facet_var]], labeller=labeller(Population=population_labels, Sex=sexlabels)) + labs(col="")
  }
  
  return(subtypeplot)
}

plot_degplot <- function(data, x_var, y_var, col_var, shape_var){
  
  x_var = enquo(x_var)
  y_var = enquo(y_var)
  col_var = enquo(col_var)
  shape_var = enquo(shape_var)
  
  degplot = ggplot2::ggplot(data, aes(x=interaction(!!x_var), y = !!y_var,
                                      text = paste('padj: ',padj))) +
    scale_colour_viridis_c(option = "viridis", end = .90) +
    geom_point(aes(col=!!col_var, shape=!!shape_var, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
    th + labs(shape = "", size = "") + scale_x_discrete(labels=subpopulation_labels)
  
  return(degplot)
}

plot_volcanoplot <- function(data, il_genes, x_var, y_var, label_var){
  
  x_var = enquo(x_var)
  y_var = enquo(y_var)
  label_var = enquo(label_var)
  
  volcanoplot = ggplot2::ggplot(data = data, aes(x = !!x_var, y = !!y_var)) +
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
                     aes(label = !!label_var),
                     force = 2,
                     nudge_y = 1) +
    scale_colour_manual(values = cols) +
    scale_x_continuous(breaks = c(seq(-10, 10, 2)),
                       limits = c(-10, 10)) + theme_line +
    labs(y="-log10(FDR)")
  
  return(volcanoplot)
  
}


#' a function that gets the median expression
get_median <- function(df, sex, var, ...){
  ifelse(sex=="Both", data <- df %>% dplyr::group_by(...) %>%
           dplyr::summarise(expression=median(as.double({{var}}))), df <- tcounts %>% dplyr::group_by(..., Sex) %>%
           dplyr::summarise(expression=median(as.double({{var}}))))
  return(data)
}


####################### UI functions ########################################
#' The UI function for each individual dataset page
#'
#' @param datasetTitle
#' @param includedeg
#' @param volc_title
#' @param include_subtype
#' @param des_dir
#' @param image_dir
#' @param includegoi
shinypageUI <- function(id, datasetTitle, includedeg = FALSE, volc_title = NULL,
                        include_subtype = FALSE, des_dir = NULL, image_dir = NULL, includegoi = TRUE, 
                        include_count = TRUE) {
  
  shiny::fluidPage(
    tags$head(
      tags$style(
        HTML("
        .container{
      width: 100%;
      height: 100%;
      margin: 0;
      padding: 0;
    }
      ")
      )
    ),
    
    shiny::fluidRow(
      class = "container",
      shinydashboard::box(title = datasetTitle,
                          width = 12,
                          status = "primary",
                          solidHeader = TRUE,
                          if (is.null(des_dir) == FALSE) {includeHTML(des_dir)},
                          if (is.null(image_dir) == FALSE) {img(src = image_dir, height = 150, width = 400)}
      )
    ),
    br(), # only appear if data is loaded
    fluidRow(class = "container",
             column(width = 10, 
                    shiny::actionButton(NS(id, "load"), "Plot Graphs", icon = icon("play-circle")),
                    shiny::downloadButton(NS(id,"report"), "Generate Code"),
                    shiny::downloadButton(NS(id,"plots"), "Download Data"),
                    shiny::helpText(em("Load data before plotting."))
                    )
                   
      
    ),
    
    if (include_count == TRUE){
      shinydotUI(id)
    },
    br(),
    if (include_subtype == TRUE) {
      shinysubtypeUI(id)
    },br(),
    if (includedeg == TRUE) {
      shinydegplotUI(id, "Volcano Plot")
    },br(),
    if(includegoi == TRUE) {
      shinygoitableUI(id)
    },br(),
    if(includedeg == TRUE) {
      shinycontrasttabUI(id)
    }
  ) # main panel closure
}


shinydotUI <- function(id) {
  shiny::fluidRow(
    class = "container",
    shinydashboard::box(title = "Naive",
                        status = "primary",
                        solidHeader = TRUE,
                        plotly::plotlyOutput(NS(id, "bulkseq_dots")),
                        height = "36em"
    ),
    shinydashboard::box(title = "Injury", status = "primary", solidHeader = TRUE,
                        shiny::plotOutput(NS(id, "bulkseq_lines")),
                        height = "36em"
    ), br(), br()
  )
}

shinyscrnapageUI <- function(id, datasetTitle, des_dir = NULL, image_dir = NULL) {
  shiny::fluidPage(
    tags$head(
      tags$style(
        HTML("
        .container{
      width: 100%;
      height: 100%;
      margin: 0;
      padding: 0;
    }
      ")
      )
    ),
    shiny::fluidRow(class = "container",
      box(title = datasetTitle,
                        width = 12,
                        status = "primary",
                        solidHeader = TRUE,
                        if (is.null(des_dir) == FALSE) {includeHTML(des_dir)},
                        if (is.null(image_dir) == FALSE) {img(src = image_dir, height = 150, width = 400)}
    )),br(),
    
    fluidRow(class = "container",
      column(width = 10,
             shiny::actionButton(NS(id, "load"), "Plot Graphs", icon = icon("play-circle")),
             shiny::downloadButton(NS(id,"scrna_plots"), "Download Plots"),
             shiny::helpText(em("Load data before plotting.")))
    ),
    shiny::fluidRow(class = "container",
      shinydashboard::box(width=6,
                          solidHeader = TRUE,
                          title = "Naive", status = "primary",
                          shiny::plotOutput(NS(id, "scrna_dots"))
      ),
      shinydashboard::box(width=6,
                          solidHeader = TRUE,
                          title = "UMAP", status = "primary",
                          fluidRow(
                            column(12,
                                   shiny::plotOutput(NS(id, "scrna_umap")))
                          )
      )
    ), br(),br(),
    shiny::fluidRow(class = "container",
      shinydashboard::box(width=6,
                          height = 550,
                          title = "Feature Plot", status = "primary",
                          solidHeader = TRUE,
                          selectizeInput(
                            inputId = NS(id,"scgeneid"),
                            label = NULL,
                            multiple = FALSE,
                            choices = NULL,
                            selected = ""
                          ),
                          plotly::plotlyOutput(NS(id, "scrna_feature"))
      )
    ),
    shiny::fluidRow(class = "container",
      shiny::plotOutput(NS(id, "scrna_violin"))
    )
  )
}

shinysubtypeUI <- function(id) {
  shiny::fluidRow(
    class = "container",
    shinydashboard::box(width = 12,
                        title = "Subtype Results", status = "primary",
                        solidHeader = TRUE,
                        shiny::plotOutput(NS(id, "bulkseq_lines_subtype"))
    )
  )
}

shinydegplotUI <- function(id, volc_title) {
  shiny::fluidRow(
    class = "container",
    shinydashboard::box(width = 6,
                        title = "Differential Gene Analysis",
                        status = "primary",
                        solidHeader = TRUE,
                        uiOutput(NS(id,"textone")),
                        plotly::plotlyOutput(NS(id, "deg_plot")),
                        height = "38em"),
    shinydashboard::box(width = 6,
                        title = volc_title, status = "primary",
                        solidHeader = TRUE,
                        actionButton(NS(id, "plotvolc"), "Plot Volcano Graphs"),
                        selectInput(NS(id, "volc_pop"), "",
                                    choices = NULL,
                                    selected = ""),
                        shiny::plotOutput(NS(id, "volcanoplot"), height = "26em"), height = "38em"
    )
    
  )
}

shinycontrasttabUI<- function(id) {
  shiny::fluidRow(
    class = "container",
    shinydashboard::box(width = 12,
                        title = "Differential Analysis Table",
                        status = "primary",
                        solidHeader = TRUE,
                        selectInput(NS(id, "contrast"), "",
                                    choices = NULL,
                                    selected = ""),
                        DT::dataTableOutput(NS(id, "contrast_table"))
    )
  )
}

shinygoitableUI <- function(id) {
  shiny::fluidRow(
    class = "container",
    shinydashboard::box(
      width = 12,
      title = "Result Table",
      status = "primary",
      solidHeader = TRUE,
      DT::dataTableOutput(NS(id,"goi_table"))
    )
  )
}

plothomescdot_ui <- function(id, dataset, combined=FALSE) {
  shiny::fluidRow(
    column(8, offset = 2,
           shiny::plotOutput(NS(id, "home_scrna_dots")))
  )
}

plotcombine_ui <- function(id) {
  shiny::fluidRow(
    column(12,
           shiny::plotOutput(NS(id,"dot"))
    )
  )
}

deg_combine_ui <- function(id) {
  shiny::fluidRow(
    column(12,
           shiny::plotOutput(NS(id,"deg")))
  )
}



