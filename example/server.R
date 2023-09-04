##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

scRNA_dir   = "data/drg.combined.rds"
theme_dir   = "global.R"
network_dir = "data/network.RData"
# data_dir = "data/drg.RData"
# gene_dir = "data/gene.RData"

des_df = data.frame(
  Dataset = c("Subtype DRG (Barry)", "Mouse DRG Bulk", "Rat DRG Bulk",
              "Human iPSC HSN1", "Human Skin Diabetes", "Human skin Carpal Tunnel (CTS)",
              "Mouse DRG Subtype (Zheng)","Human DRG RNAseq", "Human DRG spatial-seq",
              "Rat DRG HIV", "Rat DRG SNT (L5)", "Rat DRG SNI", "Rat DRG bone-cancer"),
  pain = c("ipsi", "SNI", "SNT", "patient", "PDPN", "pre_Surgery", "", "P", "", "HIV", "SNT_L5vsSHAM", "SNI_adult", "bone_cancer"),
  Pain_Model = c("Spared Nerve Injury", "Spared Nerve Injury", "Spared Nerve Injury", "Neuropathic Pain", "Diabetic Peripheral Neuropathy", "Neuropathic Pain",
                 "N/A", "N/A", "Neuropathic Pain", "HIV", "Spinal Nerve Transection", "Spared Nerve Injury", "Bone Cancer"),
  Species = c("mouse", "mouse", "rat", "human", "human", "human", "mouse", "human", "human", "rat", "rat", "rat", "mouse"),
  count = c("TPM_subtype", "TPM_mouse", "TPM_rat", "TPM_ipsc", "TPM_HS_diabetes", "TPM_HS_CTS", "TPM_zheng", "TPM_humandrg", "", "", "", "", ""),
  include_subtypes = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,  FALSE, FALSE, FALSE, FALSE),
  include_degs = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
  Name = c("subtype", "mouse", "rat", "ipsc", "HS_diabetes", "HS_CTS", "zheng", "humandrg", "scrna", "hiv", "sntl5", "sni", "bone_cancer"),
  injury_data = c("bulkseq_mat", "", "", "", "", "", "", "", "", "", "", "", ""),
  deg_df_name = c("subtype_deg_df", "mouse_deg_df", "rat_deg_df", "ipsc_deg_df", "db_deg_df", "cts_deg_df", "", "", "", 
                  "HIV_deg_df", "SNT_L5vsSHAM_deg_df", "SNI_adult_deg_df", "bone_cancer_deg_df"),
  include_gois = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  colDatas = c("bulkseq_colData", "TPM_mouse_colData", "TPM_rat_colData", "ipsc_colData", "db_colData", "skin_colData",
               "zheng_colData", "humandrg_colData", "", "", "", "", ""),
  type = c("bulk","bulk","bulk","bulk","bulk","bulk","bulk","bulk", "single", "bulk", "bulk", "bulk", "bulk"),
  Tissue = c("DRG", "DRG", "DRG", "iPSC", "Skin", "Skin", "DRG", "DRG", "DRG", "DRG", "DRG","DRG", "Bone"), 
  include_count = c(TRUE, TRUE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE,TRUE, FALSE, FALSE, FALSE, FALSE))


source("functions.R")


################################################ SERVER ##############################################################################
#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
shinyServer(function(input, output, session) {
  
  
  
  # for RData
  {
    # load(data_dir)
    # load(gene_dir)
    # allgenes = GeneData[["genebank"]]$all_genes
  }
  
  
  # For Database
  {
    conn <- dbConnect(RSQLite::SQLite(), "test.db")

    if (attempt::is_try_error(conn)){
      # Notify the user
      send_notification("Could not connect")
    } else {
      # Continue computing if the connection was successful
      human_gene_data = dbReadTable(conn, "human_gene_data")
      mouse_gene_data = dbReadTable(conn, "mouse_gene_data")
      rat_gene_data = dbReadTable(conn, "rat_gene_data")
    }
    allgenes = as.data.frame(dbReadTable(conn, "genebank"))$all_genes
  }
 
  
  
  bulk_df = des_df[des_df$type == "bulk",]
  
  output$shinypages  = shiny::renderUI({
    j = input$dataset_table_rows_selected
    if (is.null(j) == FALSE) {
      select = des_df[j,]
      
      if (select$type == "single"){
        shinyscrnapageUI("spatpage", "Human DRG spatial-seq", des_dir="des/spatial.Rhtml")
      }
      
      else {
        name = select$Name
        des_dir = paste("des/",select$Name, ".Rhtml", sep = "")
        shinypageUI(name, select$Dataset,
                    select$include_degs,
                    include_subtype = select$include_subtypes, des_dir = des_dir, includegoi = select$include_gois, 
                    include_count = select$include_count)
      }
    }
    

  })
  
  variables = reactiveValues(
    pbmc = NULL,
    indrows = NULL,
    option = NULL
  )
  
  load(network_dir)
  
  if (scRNA_dir != ""){
    variables$pbmc = readRDS(scRNA_dir) # the scRNA data
  }
  
  
  # read col_lists 

  ## for RData
  {
    # col_list = Seq_Data[[2]]
  }
  
  
  
  ## for Database
  {
    col_list <- vector(mode = "list", nrow(des_df))
    names(col_list) = des_df$colDatas

    for (i in c(1:nrow(des_df))) {
      row = des_df[i,]
      cd = row$colDatas
      ds_name = row$Dataset

      if (cd != ""){
        colD = dbReadTable(conn, cd)
        rownames(colD) = colD$row_names
        colD = colD[,-1] # remove the first col
        colD$Dataset = rep(ds_name, nrow(colD))

        if (!"Population" %in% colnames(colD)){
          colD$Population = rep(ds_name, nrow(colD))
        }

        if (!"Sex" %in% colnames(colD)) {
          colD$Sex = rep("mixed", nrow(colD))
        }

        if (!"Species" %in% colnames(colD)) {
          colD$Species = rep(row$species, nrow(colD))
        }

        col_list[[i]] = colD

      }

    }

    col_list = Filter(Negate(is.null), col_list) # remove null values
  }
  
  
  # gene search bar
  shiny::updateSelectizeInput(session,
                              inputId = "geneid",
                              label = "Search Genes:",
                              choices = allgenes,
                              server = TRUE,
                              selected = c("Trpv1", "Scn10a","Atf3")
  )
  
  # change to the plotting 
  shiny::observeEvent(input$metaa, {
    
    Newtab <- switch(input$tabs,
                     "tabhome" = "tabmeta"
                     
    )
    
    # if length(col_list) != 0, then the data loading has finished
    if (length(col_list) != 0) {
      shinydashboard::updateTabItems(session, "tabs", Newtab) 
    }
    
  })
  
  
  observeEvent(input$citation, {
    newvalue <- "tabcode"
    observe({
      updateTabItems(session, "tabs", newvalue)
    })
  })
  
  
  
  # change to the plotting 
  shiny::observeEvent(input$inda, {
    
    Newtab <- switch(input$tabs,
                     "tabhome" = "tabdata"
                     
    )
    if (length(col_list) != 0) {
      shinydashboard::updateTabItems(session, "tabs", Newtab) 
    }
    
  })
  
  # change to the plotting 
  shiny::observeEvent(input$neta, {
    
    Newtab <- switch(input$tabs,
                     "tabhome" = "tabnet"
                     
    )
    if (length(col_list) != 0) {
      shinydashboard::updateTabItems(session, "tabs", Newtab) 
    }
    
  })
  
  # a reactive variable that records whether a file is uploaded.
  rv <- reactiveValues(
    clear = FALSE,
    data = FALSE
  )
  
  shiny::observeEvent(input$file, {
    rv$clear <- FALSE
    rv$data = TRUE
  }, priority = 1000)
  
  # clear the file space when files are removed
  shiny::observeEvent(input$reset, {
    rv$clear = TRUE
    rv$data = FALSE
  })
  
  # read into gene ids from files
  file_input <- reactive({
    if (rv$clear == TRUE) {
      return(NULL)
    }
    if(rv$clear==FALSE && rv$data == TRUE) {
      goi = read.table(input$file$datapath)
      rownames(goi) <- goi[,1]
      goi <- goi[which(rownames(goi) %in% allgenes==TRUE),]
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
        goi <- goi[which(rownames(goi) %in% allgenes==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
    }
    else {
      genes = file_input()
    }
    
    # retrieve count data of selected datasets 
    
    # for database
    {
      df_list = vector(mode = "list", nrow(bulk_df))
      names(df_list) = bulk_df$colDatas
    
      selected_data = c(1:nrow(bulk_df))
      
      # selected_data = input$meta_table_selected_rows
      for (i in selected_data){
        curr = bulk_df[i,]
        d = curr$count
        if (d!=""){
          if (curr$Species == "human"){
            sql = paste("SELECT * FROM", d, "WHERE mgi_symbol = ?")
          }
          else {
            sql = paste("SELECT * FROM", d, "WHERE symbol = ?")
          }
          count_data = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(symbol=genes)))
          df_list[[i]] = count_data
        }
      }
      df_list = Filter(Negate(is.null), df_list)
    }
    
    # for RData
    {
      # df_list = list()
      # for (i in c(1:nrow(bulk_df))){
      #   curr = bulk_df[i,]
      #   count_data = Seq_Data[[1]][[curr$count]]
      #   if (curr$Species == "human"){
      #     count_data = count_data[count_data$mgi_symbol %in% genes,]
      #   }
      #   else {
      #     count_data = count_data[count_data$symbol %in% genes,]
      #   }
      #   X = rownames(count_data)
      #   count_data = cbind(X, count_data)
      #   df_list[[curr$colDatas]] = count_data
      # }
    }

    # generate a combined data frame containing count data from all datasets
    final_df = generate_combine_dataset(df_list, genes, input$sex, col_list, bulk_df$pain[selected_data], bulk_df$Species[selected_data], 
                                        expression, Condition, Population, symbol, Dataset, Species)
    
    
    
    combine_dot = plotcombine_server("dot", final_df, input$sex, genes, Population, symbol, expression, "Dataset") #dotplot
    
    if (is.null(variables$pbmc) != TRUE){
      scrna_dot = plothomescdot_server("homespat", variables$pbmc, genes)
    }
    
    
    # query deg_dfs
    
    ## for Database 
    {
      degm = RSQLite::dbGetPreparedQuery(conn, "SELECT * FROM mouse_all_deg_df WHERE symbol = ?", bind.data=data.frame(symbol=genes))
      sql = "SELECT * FROM rat_gene_data WHERE mgi_symbol = ?"
      rg = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol= genes)))$rgd_symbol
      degr = RSQLite::dbGetPreparedQuery(conn, "SELECT * FROM rat_deg_df WHERE symbol = ?", bind.data=data.frame(symbol=rg))
      sql = "SELECT * FROM human_gene_data WHERE mgi_symbol = ?"
      hg = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol= genes)))$hgnc_symbol
      degh = RSQLite::dbGetPreparedQuery(conn, "SELECT * FROM human_deg_df WHERE symbol = ?", bind.data=data.frame(symbol=hg))
    }
   
    
    ## for RData
    {
      # mouse_degdf = Seq_Data[[3]][["mouse_all_deg_df"]]
      # human_degdf = Seq_Data[[3]][["human_all_deg_df"]]
      # rat_degdf = Seq_Data[[3]][["rat_all_deg_df"]]
      # 
      # rg = rat_gene_data[rat_gene_data$mgi_symbol %in% genes,]$rgd_symbol
      # hg = human_gene_data[human_gene_data$mgi_symbol %in% genes,]$hgnc_symbol
      # degm = mouse_degdf[mouse_degdf$symbol %in% genes,]
      # degr = rat_degdf[rat_degdf$symbol %in% rg,]
      # degh = human_degdf[human_degdf$symbol %in% hg,]
    }
   
    combine_degplot = deg_combine_server("deg_plot",degr, degm, degh, Population,symbol, log2FoldChange,sig, "Dataset")
    
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
        
        # csv files for gene count data
        write.csv(final_df,"result_table.csv", row.names = FALSE, col.names = TRUE)
        
        fs = c(fs, "result_table.csv")
        
        zip(file, fs)
      }
    )
    
    
    output$combineplot <- downloadHandler(
      filename = "RNAseq_homeplot.html",
      content = function(file) {
        
        params <- list(matrix = final_df, genes = genes, sex = input$sex, scRNA = variables$pbmc, ratdeg = degr,
                       mousedeg = degm, humandeg = degh,
                       human_gene_data = human_gene_data)
        
        shiny::withProgress(value = 0,
                            message = 'Rendering plotting report',
                            detail =  'This might take several minutes.',{
                              rmarkdown::render(input = "reports/combineplot.Rmd",
                                                params = params,
                                                output_file = file,
                                                envir = new.env(parent = globalenv()))
                            })
        
      }
    ) # closure for download handler
    
    
  }) # closure for input$load
  
  output$meta_table = DT::renderDataTable({
    DT::datatable(
      bulk_df[c("Dataset", "Pain_Model", "Species", "Tissue")],
      width = 6,
      class = 'nowrap',
      options = list(scrollX = TRUE, scrollY = TRUE, pageLength = 5), 
      selection = "none"
      # selection = list(mode = "multiple", selected = c(1:nrow(bulk_df)), target = "row")
    )
  })
  
  ############### Pain Networks ####################################################
  # select genes
  shiny::updateSelectizeInput(session,
                              inputId = "gene_symbols",
                              label = "Search Genes:",
                              choices = human_gene_data$hgnc_symbol,
                              server = TRUE,
                              selected = c("ADCYAP1")
  )
  
  observeEvent(input$submit, {
    # Query gene interactions
    genes_list = as.vector(input$gene_symbols)
    genes = paste(genes_list, collapse = "%0d")
    interactions <- query_interactions(genes, 15)
    proteins = unique(c(unique(interactions$preferredName_A), unique(interactions$preferredName_B)))
    nodes = data.frame()
    id = unique(c(unique(interactions$stringId_A), unique(interactions$stringId_B)))
    nodes= as.data.frame(id)
    nodes$label = proteins
    
    
    if (input$pop == "Composite LFC") {
      mode = "meta"
      metrics = score_df
      dt = score_df[score_df$symbol %in% proteins,]
      tit = "CLFC" # title for legend
      
      output$netlegend = renderUI({
        HTML(paste0(
          "<br><br>",
          "<p style = 'display: block; margin-left: 40px; margin-right: auto;'>&emsp;&emsp;&emsp;PES</p>
                      <img style = 'display: block; margin-left: auto; margin-right: auto;' src='CES.png' width = '50' height = '180'>
                      ",
          "<br>")
        )
      })
    }
    
    else if (input$pop == "Predicted Pain Score") {
      mode = "pes"
      metrics = pain_score
      dt = pain_score[pain_score$symbol %in% proteins,]
      tit = "PPS" # title for legend
      
      output$netlegend = renderUI({
        HTML(paste0(
          "<br><br>",
          "<p style = 'display: block; margin-left: 40px; margin-right: auto;'>&emsp;&emsp;&emsp;PES</p>
                      <img style = 'display: block; margin-left: auto; margin-right: auto;' src='CES.png' width = '50' height = '180'>
                      ",
          "<br>")
        )
      })
    }
    
    else {
      mode = "indiv"
      # the dict contains deg_df based on the name of datasets
      
      res = network_df[network_df$experiment_id == input$pop,]
      res = mutate(res, sig=ifelse((res$padj<0.05), "SIG", "NS"))
      metrics = res
      dt = filter(res, symbol %in% proteins)  # this means the proteins must be in human symbol
      rownames(dt) = dt$symbol
      dt = dt[c("log2FoldChange","padj", "sig")]
      tit = "log2FoldChange"
      
      output$netlegend = renderUI({
        HTML(paste0(
          "<br><br>",
          "<p style = 'display: block; margin-left: 40px; margin-right: auto;'>&emsp;&emsp;&emsp;LFC</p>
                      <img style = 'display: block; margin-left: auto; margin-right: auto;' src='CES.png' width = '50' height = '180'>
                      ",
          "<br>")
        )
      })
    }
    
    
    
    # network construction
    nodes = get_nodes(interactions, genes_list, metrics, mode=mode, snps, pg1, pg2)
    edges <- data.frame(from = interactions$stringId_A, to = interactions$stringId_B, width = interactions$score,
                        color = "#444444")
    
    lnodes <- data.frame(shape = c( "icon"), 
                         icon.color = c("#961252", "#d8aec4ff", "lightgrey", "#aed5eaff","#0277bd"),
                         icon.code  = c("f35b", "f111", "f111", "f111", "f358"),
                         icon.size = 30) 
    
    # network visualisation
    nw = visNetwork::visNetwork(nodes, edges, height = "500px", width = "100%") %>%
      visPhysics(enabled = FALSE,solver = "forceAtlas2Based",
                 forceAtlas2Based = list(avoidOverlap = 1)) %>%
      visEdges(smooth = FALSE) %>%
      visInteraction(navigationButtons = TRUE, zoomSpeed = 0.6) %>%
      visLegend(addNodes = lnodes, 
                #main = list(text = "Enrichment", style = "font-family:Gadugi"),
                useGroups = FALSE, 
                width = 0.2, 
                position = "right",
                ncol = 1,
                stepY = 50,
                #stepX = 50,
                zoom = FALSE)
    output$network <- renderVisNetwork({nw})
    
    output$downloadnet <- downloadHandler(
      filename = "graph.html",
      content = function(file) {
        visSave(nw, file)
      }
    )
    
    
    output$contrast_table <- DT::renderDataTable({a
      DT::datatable(
        dt,
        width = 6,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 8)
      )
    })
    
    # a table containing description and information of genes included
    output$protein_table <- DT::renderDataTable({
      DT::datatable(
        query_information(proteins),
        width = 6,
        class = "wrap",
        options = list(scrollX = TRUE, pageLength = 8)
      )
    })
    
  })
  ##################################################################################################
  
  
  
  
  
  
  
  # 
  # ############### Pain Networks ####################################################
  # # select genes
  # shiny::updateSelectizeInput(session,
  #                             inputId = "gene_symbols",
  #                             label = "Search Genes:",
  #                             choices = human_gene_data$hgnc_symbol,
  #                             server = TRUE,
  #                             selected = c("GRIN1", "LYN", "KCND2")
  # )
  # 
  # observeEvent(input$submit, {
  #   # Query gene interactions
  #   genes_list = as.vector(input$gene_symbols)
  #   genes = paste(genes_list, collapse = "%0d")
  #   interactions <- query_interactions(genes, 15)
  #   proteins = unique(c(unique(interactions$preferredName_A), unique(interactions$preferredName_B)))
  #   nodes = data.frame()
  #   id = unique(c(unique(interactions$stringId_A), unique(interactions$stringId_B)))
  #   nodes= as.data.frame(id)
  #   nodes$label = proteins
  #   
  # 
  #   if (input$pop == "Predicted Pain Score") {
  #     mode = "meta"
  #     metrics = pain_score
  #     dt = pain_score[pain_score$symbol %in% proteins,]
  #     tit = "PPS" # title for legend
  #     
  #     output$netlegend = renderUI({
  #       HTML(paste0(
  #         "<br><br>",
  #         "<p style = 'display: block; margin-left: 40px; margin-right: auto;'>&emsp;&emsp;&emsp;PES</p>
  #                     <img style = 'display: block; margin-left: auto; margin-right: auto;' src='CES.png' width = '50' height = '180'>
  #                     ",
  #         "<br>")
  #       )
  #     })
  #   }
  #   
  #   else if (input$pop == "Composite LFC") {
  #     mode = "comp"
  #     metrics = score_df
  #     dt = score_df[score_df$symbol %in% proteins,]
  #     tit = "Composite LFC" # title for legend
  #     
  #     output$netlegend = renderUI({
  #       HTML(paste0(
  #         "<br><br>",
  #         "<p style = 'display: block; margin-left: 40px; margin-right: auto;'>&emsp;&emsp;&emsp;PES</p>
  #                     <img style = 'display: block; margin-left: auto; margin-right: auto;' src='CES.png' width = '50' height = '180'>
  #                     ",
  #         "<br>")
  #       )
  #     })
  #   }
  #   
  #   else {
  #     mode = "indiv"
  #     # the dict contains deg_df based on the name of datasets
  #     
  #     res = network_df[network_df$experiment_id == input$pop,]
  #     res = mutate(res, sig=ifelse((res$padj<0.05), "SIG", "NS"))
  #     metrics = res
  #     dt = filter(res, symbol %in% proteins)  # this means the proteins must be in human symbol
  #     rownames(dt) = dt$symbol
  #     dt = dt[c("log2FoldChange","padj", "sig")]
  #     tit = "log2FoldChange"
  #     
  #     output$netlegend = renderUI({
  #       HTML(paste0(
  #         "<br><br>",
  #         "<p style = 'display: block; margin-left: 40px; margin-right: auto;'>&emsp;&emsp;&emsp;LFC</p>
  #                     <img style = 'display: block; margin-left: auto; margin-right: auto;' src='CES.png' width = '50' height = '180'>
  #                     ",
  #         "<br>")
  #       )
  #     })
  #   }
  #   
  #   
  #   
  #   # network construction
  #   nodes = get_nodes(interactions, genes_list, metrics, mode=mode, snps, pg1, pg2)
  #   edges <- data.frame(from = interactions$stringId_A, to = interactions$stringId_B, width = interactions$score,
  #                       color = "#444444")
  #   
  #   lnodes <- data.frame(shape = c( "icon"), 
  #                        icon.color = c("#961252", "#d8aec4ff", "lightgrey", "#aed5eaff","#0277bd"),
  #                        icon.code  = c("f35b", "f111", "f111", "f111", "f358"),
  #                        icon.size = 30) 
  # 
  #   # network visualisation
  #   nw = visNetwork::visNetwork(nodes, edges, height = "500px", width = "100%") %>%
  #     visPhysics(enabled = FALSE,solver = "forceAtlas2Based",
  #                forceAtlas2Based = list(avoidOverlap = 1)) %>%
  #     visEdges(smooth = FALSE) %>%
  #     visInteraction(navigationButtons = TRUE, zoomSpeed = 0.6) %>%
  #     visLegend(addNodes = lnodes, 
  #               #main = list(text = "Enrichment", style = "font-family:Gadugi"),
  #               useGroups = FALSE, 
  #               width = 0.2, 
  #               position = "right",
  #               ncol = 1,
  #               stepY = 50,
  #               #stepX = 50,
  #               zoom = FALSE)
  #   output$network <- renderVisNetwork({nw})
  #   
  #   output$downloadnet <- downloadHandler(
  #     filename = "graph.html",
  #     content = function(file) {
  #       visSave(nw, file)
  #     }
  #   )
  #  
  #   output$contrast_table <- DT::renderDataTable({a
  #     DT::datatable(
  #       dt,
  #       width = 6,
  #       class = 'nowrap',
  #       options = list(scrollX = TRUE, pageLength = 8)
  #     )
  #   })
  #   
  #   # a table containing description and information of genes included
  #   output$protein_table <- DT::renderDataTable({
  #     DT::datatable(
  #       query_information(proteins),
  #       width = 6,
  #       class = "wrap",
  #       options = list(scrollX = TRUE, pageLength = 8)
  #     )
  #   })
  #   
  # })
  # ##################################################################################################
  
  
  # update the table rows selected
  shiny::observe({
    variables$indrows =input$dataset_table_rows_selected
    
  })
  
  output$dataset_table = DT::renderDataTable({
    DT::datatable(
      des_df[c("Dataset", "Pain_Model", "Species", "Tissue")],
      width = 6,
      class = 'nowrap',
      options = list(scrollX = TRUE, scrollY = TRUE, pageLength = 5), selection = list(mode = "single", selected = c(1), target = 'row')
    )
  })
  

  
  # make sure, no such injury_table: NULL
  shiny::observe({
    i = input$dataset_table_rows_selected
    
    if (is.null(i) == FALSE) {
      selected = des_df[i,]
      name = selected$Name
      
      if (i==1){
        
        # for database
        {
          shiny::callModule(shinypage, name, input$sex, selected$count, col_list[[selected$colDatas]],
                            selected$Species, selected$pain, reactive({genes()}), include_deg = selected$include_degs, deg_df_name = selected$deg_df_name,
                            dataset = selected$Dataset, include_subtype = selected$include_subtypes, include_count = selected$include_count, injury_data = "bulkseq_mat", parent = session)

        }
        
        # for RData
        {
          # shiny::callModule(shinypage, name, input$sex, Seq_Data[[1]][[selected$count]], col_list[[selected$colDatas]],
          #                   selected$Species, selected$pain, reactive({genes()}), include_deg = selected$include_degs, deg_df_name = DEG_Data[[selected$deg_df_name]],
          #                   dataset = selected$Dataset, include_subtype = selected$include_subtypes, include_count = selected$include_count, injury_data = Seq_Data[[1]][["bulkseq_mat"]], parent = session)
          # 
        }
        
        }
      if (selected$Dataset == "Human DRG spatial-seq"){
        shiny::callModule(shinyscrna, "spatpage", reactive({genes()}), variables$pbmc)
      }
      if ((i!=1) & (i!=9)){
        
        # for Database
        {
          shiny::callModule(shinypage, name, input$sex, selected$count, col_list[[selected$colDatas]],
                            selected$Species, selected$pain, reactive({genes()}), include_deg = selected$include_degs, deg_df_name = selected$deg_df_name,
                            dataset = selected$Dataset, include_subtype = selected$include_subtypes, include_count = selected$include_count, parent = session)

        }
       
        # for Rdata
        {
          # shiny::callModule(shinypage, name, input$sex, Seq_Data[[1]][[selected$count]], col_list[[selected$colDatas]],
          #                   selected$Species, selected$pain, reactive({genes()}), include_deg = selected$include_degs, deg_df_name = DEG_Data[[selected$deg_df_name]],
          #                   dataset = selected$Dataset, include_subtype = selected$include_subtypes, include_count = selected$include_count, parent = session)
          # 
        }
       
        }
    }
    
  })
  
  ### leaflet map for contact details
  output$myMap <- renderLeaflet({
    m <- leaflet() %>% addTiles()
    m <- m %>% setView( -1.238233, 51.756192, zoom = 13)
    m %>% addPopups(-1.2217, 51.76529, "Neural Injury Group")
  })
  
  # disconnect db when stops
  onStop(function() {
    dbDisconnect(conn)
  })
  
  
})


options(warn=-1) # remove warnings



