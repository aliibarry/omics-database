##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

### Deployment, setRepositories may help streamline
# rsconnect::deployApp()
# setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.16/bioc"))

# add libraries 
{
  library(visNetwork)
  library(cowplot)
  library(shiny)
  library(shinythemes)
  library(shinydashboard)
  library(shinyFeedback)
  library(plotly)
  library(tidyverse)
  library("RColorBrewer")
  library(stringr)
  library(ggplot2)
  library(viridis)
  library("cowplot")
  library("fontawesome") 
  library(shinyalert)
  library(RSQLite)
  library(shinybusy)    # for the loader
}

#################################################################

### use as needed (## set data input in 'functions_base.R')
table_dir = "des/table.html"
sql_dir     = "./test.db"
scRNA_dir   = NULL # write dir for scRNA data
network_dir = "./network.RData" #if network data added

#################################################################

if (is.null(scRNA_dir) == FALSE) {
  library(Seurat) # for scRNA; only load if scRNA data is present 
}

source("functions_base.R")


{
  des_df = data.frame(
    Dataset      = c("trial A", "trial B"),
    Name         = c("trial A", "trial B"),
    colDatas     = c("cd_mouseA", "cd_mouseB"),
    Species      = c("mouse", "mouse"),
    type         = c("bulk", "bulk"), #bulk or single
    Tissue       = c("DRG", "brain"),
    Model        = c("injury", "naive"),
    count        = c("db_A", "db_B"), #count files
    include_degs = c(TRUE, FALSE), #boolean
    deg_df_name  = c("deg_A", ""),
    include_gois = c(TRUE, TRUE),  #boolean
    include_subtypes = c(TRUE, FALSE), #boolean
    extra_data   = c("",""), #leave "" where not needed
    pain         = c(" "))   #extra notes
}

ui = function(req){fluidPage(
  useShinyFeedback(),
  # includeCSS("www/style.css"), # can add style sheet to modify aes
  shinyFeedback::useShinyFeedback(),
  
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title="Title", titleWidth = 225
    ),
    shinydashboard::dashboardSidebar(width = 225,
                                     sidebarMenu(
                                       id = "tabs",
                                       
                                       shiny::selectizeInput(
                                         inputId = "geneid",
                                         label = "Search Genes:",
                                         multiple = TRUE,
                                         choices = NULL
                                       ),
                                       
                                       shiny::fileInput("file",
                                                        label = "Upload File:",
                                                        accept = c(
                                                          'text/csv',
                                                          'text/comma-separated-values',
                                                          'text/tab-separated-values',
                                                          'text/plain',
                                                          '.csv', '.txt',
                                                          '.tsv')),
                                       actionButton("reset", "Clear"),
                                       shiny::selectizeInput(
                                         inputId = "sex",
                                         label = "Select Sex:",
                                         choices = c('Both', 'Separate'),
                                         selected = 'Both'),
                                       
                                       # sidebar menu for tabs (pages)
                                       shinydashboard::menuItem("Home", tabName = "tabhome", icon = icon("home"), badgeLabel = "Start Here!", badgeColor = "green"),
                                       shinydashboard::menuItem("Grouped Analysis", tabName = "tabmeta", icon = icon("dashboard")),
                                       shinydashboard::menuItem("Individual Analysis", tabName = "tabdata", icon = icon("database")),
                                       shinydashboard::menuItem("Network Analysis", tabName = "tabnet", icon = icon("circle-nodes")),
                                       br(),
                                       br()
                                     ),
    ),
    
    # Main panels for output, each TabItem is a different page, same sidebar menu
    shinydashboard::dashboardBody(
      # add a loader
      shinybusy::add_busy_spinner(spin = "fading-circle", position = "bottom-right"),
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "tabhome", 
                                
                                shinydashboard::box(
                                  width = 12,
                                  status = "primary",
                                  column(width = 6, offset = 4,
                                         h1(strong("Title"))), ## add title
                                  column(width = 8, offset = 2,
                                         h4("text"), ## add description
                                         br(),br() #breaks for spacing
                                  ),
                                  
                                  column(width = 4, offset = 4, 
                                         actionButton("meta_tab","Grouped Analysis", class = "btn-success btn-lg btn-block", icon = icon("dashboard")),
                                         actionButton("ind_tab","Individual Analysis", class = "btn-success btn-lg btn-block", icon = icon("database")),
                                         actionButton("net_tab","Network Analysis", class = "btn-success btn-lg btn-block", icon = icon("circle-nodes")),
                                         br(),br()
                                  )
                                ),
                                br()
                                
        ),
        
        shinydashboard::tabItem(tabName="tabmeta",
                                br(),
                                shiny::fluidRow(
                                  shinydashboard::box(status = "primary",
                                                      width = 12,
                                                      title = "Datasets",
                                                      solidHeader = TRUE,
                                                      
                                                      
                                                      DT::dataTableOutput("meta_table"),
                                                      
                                                      collapsible = TRUE),
                                  
                                ), br(),
                                fluidRow(
                                  column(12,offset = 0,
                                         
                                         shiny::actionButton("load", "Plot Graphs", icon = icon("play-circle")),
                                         shiny::downloadButton("combineplot", "Generate code"),
                                         shiny::downloadButton("combineplots", "Download Plots")
                                         
                                  )
                                ), br(),
                                
                                shiny::fluidRow(
                                  shinydashboard::box(width = 12,
                                                      title = "Naive",
                                                      collapsible = TRUE,
                                                      solidHeader = TRUE,
                                                      status = "primary",
                                                      plotcombine_ui("dot"),
                                                      uiOutput("scRNA")
                                  )
                                ),
                                fluidRow(
                                  br(),
                                  column(width = 12,
                                         p("NS: Non-significant, SIG: significant. Defaults: FDR > 0.05, LFC > 1. ")
                                  ),
                                  shinydashboard::box(width = 12,
                                                      title = "Differential Gene Analysis",
                                                      collapsible = TRUE,
                                                      solidHeader = TRUE,
                                                      status = "primary",
                                                      deg_combine_ui("deg_plot"))
                                )
                                
                                
        ), # tabItem meta
        
        # UI for the network page
        shinydashboard::tabItem(tabName = "tabnet",
                                fluidRow(
                                  shinydashboard::box(title = "Network Analysis",
                                                      width = 12,
                                                      status = "primary",
                                                      solidHeader = TRUE,
                                                      p("Description")
                                                      # includeHTML("des/network.Rhtml")
                                  )
                                ),br(),
                                fluidRow(
                                  column(width = 3,
                                         selectizeInput(
                                           inputId = "gene_symbols",
                                           label = "Enter Query Gene Symbols:",
                                           multiple = TRUE,
                                           choices = NULL
                                         )
                                  ),
                                  column(width = 3,
                                         selectInput("pop", "Select Experiment:",
                                                     choices = c("Composite Enrichment Score", unique(network_df$experiment_id)),
                                                     selected = "Diabetes")
                                         
                                  ),
                                  column(width = 3,br(),
                                         actionButton("submit", "Construct Network")), style = "margin-top:20px;"
                                ), 
                                fluidRow(
                                  shinydashboard::box(width = 12, status = "primary",
                                                      solidHeader = TRUE, title = "Network",
                                                      column(width = 12,
                                                             shiny::downloadButton("downloadnet", "Download"),
                                                             visNetwork::visNetworkOutput("network", height = "435px")
                                                      ), height = "38em"
                                  )
                                ), br(),
                                fluidRow(
                                  shinydashboard::box(width = 12, status = "primary",
                                                      solidHeader = TRUE, title = "Result Table",
                                                      tabsetPanel(
                                                        tabPanel( "DEG Table",
                                                                  DT::dataTableOutput("contrast_table")
                                                        ),
                                                        tabPanel("Nodes",
                                                                 DT::dataTableOutput("protein_table")
                                                        )
                                                      ),
                                                      height = "38em"
                                  )
                                ), br(),
                                
                                # a table containing experiments involved and the treatment groups
                                fluidRow(
                                  shinydashboard::box(status = "primary",
                                                      solidHeader = TRUE,
                                                      width = 12,
                                                      title = "Datasets",
                                                      column(width = 12,
                                                             p("Table describing experiments")
                                                             # includeHTML("des/datatable.html")
                                                      ),
                                                      collapsible = TRUE),
                                )
        ),
        
        # UI for individual experiments
        shinydashboard::tabItem(tabName = "tabdata",
                                fluidRow(
                                  shinydashboard::box(status = "primary",
                                                      solidHeader = TRUE,
                                                      width = 12,
                                                      title = "Choose A Dataset for Individual Analysis",
                                                      DT::dataTableOutput("dataset_table"),
                                                      collapsible = TRUE)
                                ),  br(),# fluidRow closure
                                fluidRow(
                                  uiOutput("shinypages")
                                )
                                
        )
      ) # tabItems (all)
    )  #mainpanel dashboardbody close
  ) # dashboad page close
) #fluid page
}

server = function(input, output, session) {
  
  conn <- dbConnect(RSQLite::SQLite(), sql_dir)
  load(network_dir)
  
  variables = reactiveValues(
    pbmc = NULL,
    metarows = NULL,
    indrows = NULL,
    option = NULL,
    finished = FALSE
  )
  
  if (is.null(scRNA_dir) == FALSE) {
    variables$pbmc = readRDS(scRNA_dir) # the scRNA data
  }
  
  if (attempt::is_try_error(conn)){
    send_notification("Could not connect")
  } 
  
  else {
    # Continue computing if the connection was successful
    human_gene_data = as.data.frame(dbGetQuery(conn, "SELECT * FROM human_gene_data"))
    mouse_gene_data = as.data.frame(dbGetQuery(conn, "SELECT * FROM mouse_gene_data"))
    rat_gene_data   = as.data.frame(dbGetQuery(conn, "SELECT * FROM rat_gene_data"))
  }
  
  
  allgenes = as.data.frame(dbGetQuery(conn, "SELECT * FROM genebank"))$all_genes
  bulk_df  = des_df[des_df$type == "bulk",]
  
  col_list = list()
  
  coldata = as.vector(des_df$colDatas)# names for the meta data
  for (i in c(1:length(coldata))){
    cd = coldata[i]
    
    # will be "" for scRNA data
    if (cd != ""){
      sql = paste("SELECT * FROM", cd)
      colD = as.data.frame(dbGetQuery(conn, sql))
      
      # for each description, ensure that they have dataset, pain, etc.
      cnames = colnames(colD)
      if (!"Timepoint" %in% cnames) {
        colD$Timepoint = rep("NA", nrow(colD))
      }
      
      if (!"Sex" %in% cnames) {
        colD$Sex = rep("mixed", nrow(colD))
      }
      
      if (!"Population" %in% cnames) {
        colD$Population = rep(des_df[i,]$Dataset, nrow(colD))
      }
      
      if (!"Species" %in% cnames) {
        colD$Species = rep(des_df[i,]$species, nrow(colD))
      }
      
      if (!"Dataset" %in% cnames) {
        colD$Dataset = rep(des_df[i,]$Dataset, nrow(colD))
      }
      
      if (!"Condition" %in% cnames) {
        colD$Condition = rep(des_df[i,]$disease_condition, nrow(colD))
      }
      
      rownames(colD) = colD$row_names
      colD = colD[,-1] # remove the first col
      col_list = append(col_list, list(colD))
    }
  }
  
  variables$finished = TRUE
  
  
  # gene search bar
  shiny::updateSelectizeInput(session,
                              inputId = "geneid",
                              label = "Search Genes:",
                              choices = allgenes,
                              server = TRUE,
                              selected = c("Trpv1")
  )
  
  shiny::observeEvent(input$meta_tab, {
    Newtab <- switch(input$tabs,
                     "tabhome" = "tabmeta"
                     
    )
    if (variables$finished == TRUE) {
      shinydashboard::updateTabItems(session, "tabs", Newtab) 
    }
    
  })
  
  # change to the plotting 
  shiny::observeEvent(input$ind_tab, {
    
    Newtab <- switch(input$tabs,
                     "tabhome" = "tabdata"
                     
    )
    if (variables$finished == TRUE) {
      shinydashboard::updateTabItems(session, "tabs", Newtab) 
    }
    
  })
  
  # change to the plotting 
  shiny::observeEvent(input$net_tab, {
    Newtab <- switch(input$tabs,
                     "tabhome" = "tabnet"
                     
    )
    if (variables$finished == TRUE) {
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
  
  variables$option = "none"
  

  
  # remove modal 
  shiny::observeEvent(input$dismiss_modal, {
    shiny::removeModal()
  })
  
  # only show scRNA plots if there is scRNA data 
  output$scRNA <- renderUI({
    if (is.null(variables$pbmc) != TRUE){
      plothomescdot_ui("homespat")
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
    df_list = list()
    
    if (variables$option == "none") {
      selected_data = c(1:nrow(bulk_df))
    }
    
    else {
      selected_data = variables$metarows 
    }
    
    for (i in selected_data){
      curr = bulk_df[i,]
      d = curr$count
      if (curr$Species == "human"){
        sql = paste("SELECT * FROM", d, "WHERE mgi_symbol = ?")
      }
      else {
        sql = paste("SELECT * FROM", d, "WHERE symbol = ?")
      }
      count_data = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(symbol=genes)))
      df_list = append(df_list, list(count_data))
    }
    
    
    final_df = generate_combine_dataset(df_list, genes, input$sex, col_list[selected_data], des_df$pain[selected_data], des_df$Species[selected_data], 
                                        expression, Condition, Population, symbol, Dataset, Species)
    
    
    
    # homepage dotplot
    combine_dot = plotcombine_server("dot", final_df, input$sex, genes, Population, symbol, expression, "Dataset") #dotplot
    if (is.null(variables$pbmc) != TRUE){
      scrna_dot = plothomescdot_server("homespat", variables$pbmc, genes)
    }
    
    
    # query deg_dfs
    degm = RSQLite::dbGetPreparedQuery(conn, "SELECT * FROM mouse_all_deg_df WHERE symbol = ?", bind.data=data.frame(symbol=genes))
    sql = "SELECT * FROM rat_gene_data WHERE mgi_symbol = ?"
    rg = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol= genes)))$rgd_symbol
    degr = RSQLite::dbGetPreparedQuery(conn, "SELECT * FROM rat_deg_df WHERE symbol = ?", bind.data=data.frame(symbol=rg))
    
    sql = "SELECT * FROM human_gene_data WHERE mgi_symbol = ?"
    hg = as.data.frame(dbGetPreparedQuery(conn, sql, bind.data=data.frame(mgi_symbol= genes)))$hgnc_symbol
    degh = RSQLite::dbGetPreparedQuery(conn, "SELECT * FROM human_deg_df WHERE symbol = ?", bind.data=data.frame(symbol=hg))
    
    
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
      bulk_df[c("Dataset", "Model", "Species", "Tissue")],
      width = 6,
      class = 'nowrap',
      options = list(scrollX = TRUE, scrollY = TRUE, pageLength = 5), selection = variables$option
    )
  })
  
  # individual analysis
  output$shinypages  = shiny::renderUI({
    j = input$dataset_table_rows_selected
    select = des_df[j,]
    if (select$type == "single"){
      shinyscrnapageUI("spatpage", "Human DRG spatial-seq")
    }
    
    else {
      name = select$Name
      shinypageUI(name, select$Dataset,
                  select$include_degs,
                  select$include_subtypes, des_dir = NULL, includegoi = select$include_gois)
    }
    
  })
  
  ############### Pain Networks ####################################################
  # select genes
  shiny::updateSelectizeInput(session,
                              inputId = "gene_symbols",
                              label = "Search Genes:",
                              choices = human_gene_data$hgnc_symbol,
                              server = TRUE,
                              selected = c("ATF3")
  )
  
  observeEvent(input$submit, {
    # Query gene interactions
    genes_list = as.vector(input$gene_symbols)
    genes = paste(genes_list, collapse = "%0d")
    interactions <- query_interactions(genes, 10)
    proteins = unique(c(unique(interactions$preferredName_A), unique(interactions$preferredName_B)))
    nodes = data.frame()
    id = unique(c(unique(interactions$stringId_A), unique(interactions$stringId_B)))
    nodes= as.data.frame(id)
    nodes$label = proteins
    nodes = mutate(nodes, shape = ifelse(nodes$label %in% input$gene_symbols, "diamond", "circle"))
    
    
    if (input$pop == "Composite Enrichment Score") {
      mode = "meta"
      metrics = score_df
      dt = score_df[score_df$symbol %in% proteins,]
    }
    else {
      mode = "indiv"
      # the dict contains deg_df based on the name of datasets
      res = data[data$experiment_id == input$pop,]
      res = mutate(res, sig=ifelse((res$padj<0.05), "SIG", "NS"))
      metrics = res
      dt = filter(res, symbol %in% proteins)  # this means the proteins must be in human symbol
      rownames(dt) = dt$symbol
      dt = dt[c("log2FoldChange","padj", "sig")]
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
    
    output$contrast_table <- DT::renderDataTable({
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
  
  # update the table rows selected
  shiny::observe({
    variables$metarows = as.vector(input$meta_table_rows_selected)
    variables$indrows =input$dataset_table_rows_selected
    
  })
  
  output$dataset_table = DT::renderDataTable({
    DT::datatable(
      des_df[c("Dataset", "Model", "Species", "Tissue")],
      width = 6,
      class = 'nowrap',
      options = list(scrollX = TRUE, scrollY = TRUE, pageLength = 5), 
      selection = list(mode = "single", selected = c(1), target = 'row')
    )
  })
  
  
  # individual dataset analysis 
  shiny::observe({
    i = input$dataset_table_rows_selected
    
    if (is.null(i) == FALSE) {
      selected = des_df[i,]
      name = selected$Name
      # this dataset contains an additional injury file 
      if (selected$Dataset == "Subtype DRG (Barry)"){
        shiny::callModule(shinypage, name, input$sex, selected$count, as.data.frame(col_list[i]),
                          selected$Species, selected$pain, reactive({genes()}), include_deg = selected$include_degs, deg_df_name = selected$deg_df_name,
                          dataset = selected$Dataset, include_subtype = selected$include_subtypes, extra_data = " ", parent = session)
      }
      
      if (selected$Dataset == "Human DRG spatial-seq"){
        shiny::callModule(shinyscrna, "spatpage", reactive({genes()}), variables$pbmc)
      }
      
      # for all other bulk data 
      if ((selected$Dataset != "Subtype DRG (Barry)") & (selected$type != "single")){
        shiny::callModule(shinypage, name, input$sex, selected$count, as.data.frame(col_list[i]),
                          selected$Species, selected$pain, reactive({genes()}), include_deg = selected$include_degs, deg_df_name = selected$deg_df_name,
                          dataset = selected$Dataset, include_subtype = selected$include_subtypes, parent = session)
      }
    }
    
  })
  
  # disconnect db when stops
  onStop(function() {
    dbDisconnect(conn)
  })
  
  
}

shiny::shinyApp(ui = ui, server = server)