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

members <- data.frame(name=c("App1", "App2"), nr=c("App1", "App2"))

# add libraries 
{
  
  library(visNetwork) # for plotting network 
  library(httr) # for connecting to api of string database
  library(cowplot)
  library(shiny)
  library(shinythemes) 
  library(shinydashboard) # for ui structure 
  library(shinyFeedback)
  library(plotly)
  library(tidyverse) # for data cleaning + setup 
  library("RColorBrewer") # for colored legends 
  library(stringr)
  library(ggplot2) # for plotting 
  library(viridis) # for colored legends
  library("fontawesome") # for icons
  library(Seurat) # for scRNA seq data analysis 
  library(shinyalert) # for modal boxes 
  library(RSQLite) # for database 
  library(shinybusy)    # for the loader
  library("leaflet") # for map 
  library(markdown)
}




table_dir = "des/table.html"

source("functions.R")
load("data/network.RData")

ui = function(req){fluidPage(
  tags$head(includeHTML("www/google-analytics.html")),
  useShinyFeedback(),
  includeCSS("www/style.css"),
  shinyFeedback::useShinyFeedback(),
  
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title="PRH", titleWidth = 225
    ),
    shinydashboard::dashboardSidebar(width = 225,
                                     sidebarMenu(
                                       id = "tabs",
                                       
                                       HTML(paste0( # oxford logo + ndcn link
                                         "<br><br>",
                                         "<a href='https://www.ndcn.ox.ac.uk/research/neural-injury-group' target='_blank'> <img style = 'display: block; margin-left: auto; margin-right: auto;' src='oxfordlogo2.png' width = '90'></a>",
                                         "<br>")
                                       ),
                                       
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
                                       shinydashboard::menuItem("Home", tabName = "tabhome", icon = icon("home"), badgeLabel = "Start Here!", badgeColor = "teal"),
                                       shinydashboard::menuItem("Grouped Analysis ", tabName = "tabmeta", icon = icon("dashboard")),
                                       shinydashboard::menuItem("Individual Analysis", tabName = "tabdata", icon = icon("database")),
                                       shinydashboard::menuItem("Network Analysis", tabName = "tabnet", icon = icon("circle-nodes")),
                                       shinydashboard::menuItem("User Guide", tabName = "tabcode", icon = icon("folder-open")),
                                       shinydashboard::menuItem("Contact", tabName = "tabguide", icon = icon("address-book")),
                                       #shinydashboard::menuItem("test", tabName = "tabtest", icon = icon("address-book")),
                                       br(),
                                       br()
                                     ),
                                     #Footer (icons + social media links)
                                     HTML(paste0(
                                       "<table style='clear: both;padding: 0;text-align: center;vertical-align: middle;line-height: normal;
                       
                        margin: 30px;position: fixed;bottom: 0px;width: 165px;'>", # start table 
                                       # icons placed in 6 column table
                                       "<tr >",
                                       "<td style='padding: 10px;'></td>",
                                       "<td style='padding: 5px;'><a href='https://twitter.com/aliibarry' target='_blank'><i class='fab fa-twitter fa-lg'></i></a></td>",
                                       "<td style='padding: 5px;'><a href='https://github.com/aliibarry/shiny' target='_blank'><i class='fab fa-github fa-lg'></i></a></td>",
                                       "<td style='padding: 5px;'><a href='https://www.ndcn.ox.ac.uk/research/neural-injury-group' target='_blank'><i class='fas fa-brain fa-lg'></i></a></td>",
                                       
                                       "<td style='padding: 10px;'></td>",
                                       "</tr>",
                                       
                                       # second table row, merged columns    
                                       "<tr>",
                                       "<script>","var today = new Date();","var yyyy = today.getFullYear();","</script>",
                                       "<td  colspan='6', style = 'text-align: center;background-color:#222d32ff;'><small>&copy; 
                        <a href='https://github.com/aliibarry?tab=shiny' target='_blank'>github.com</a> - <script>document.write(yyyy);</script></small></td>",
                                       "</tr>",
                                       "</table>", #end table
                                       "<br>")
                                     )
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
                                  column(width = 8, offset = 0,
                                         h1(strong("Pain RNA-seq Hub"))),
                                  column(width = 8, offset = 0,
                                         h4("This resource from the Neural Injury Group (NDCN) generates visualisation of pain-related genes both in the context of gene expression and their network associations."),
                                         h4("This currently includes a curated selection of bulk RNA-seq and scRNA-seq studies. This resource is still being updated (as of April 2024). Queries, including citations details, should be directed to Ali Barry (allison.barry@ndcn.ox.ac.uk)) "),
                                         
                                         br(),br()
                                  ),
                                  
                                  column(width = 4, offset = 4, 
                                         actionButton("metaa","Grouped Analysis ", class = "btn-success btn-lg btn-block", icon = icon("dashboard")),
                                         actionButton("inda","Individual Analysis", class = "btn-success btn-lg btn-block", icon = icon("database")),
                                         actionButton("neta","Network Analysis", class = "btn-success btn-lg btn-block", icon = icon("circle-nodes")),
                                         br(),br()
                                  ),
                                  HTML(paste0(
                                    "<table style='margin-left:auto; margin-right:auto; bottom=0;'>",
                                    "<style>
                                  		.logo-row {
                                  			display: flex; /* Use flexbox layout to arrange the logos */
                                  			justify-content: space-between; /* Distribute the logos evenly along the row */
                                  			align-items: center; /* Center the logos vertically in the row */
                                  			padding: 20px; /* Add some padding around the row */
                                  		}
                                  	</style>",
                                    "<tr class = 'logo-row'>",
                                    "<td><a href='https://www.gtc.ox.ac.uk/' target='_blank'> <img style= 'display: center;' src='gtclogo.png', width = '100'></a></td>",
                                    "<td><a href='https://wellcome.org/' target='_blank'> <img style = 'display: center;' src='wt-logo.svg'. width = '100'></a></td>",
                                    "<td><a href='https://www.ukri.org/councils/mrc/' target='_blank'> <img style = 'display: center;' src='ukri.png'. width = '90' height='90'></a></td>",
                                    "<td><img style = 'display: center;' src='painstorm.png'. width = '260' height='90'></td>",
                                    "</tr>",
                                    "</table>"
                                  )
                                  )
                                  
                                ),
                                
                                
                                br() 
        ),
        
        shinydashboard::tabItem(tabName="tabmeta",
                                br(),
                                
                                shiny::fluidRow(
                                  shinydashboard::box(status = "primary", 
                                                      width = 12, 
                                                      title = "User Guide", 
                                                      solidHeader = TRUE, 
                                                      includeMarkdown('des/usernotes.Rmd'))
                                ),
                                  br(),
                                shiny::fluidRow(
                                  shinydashboard::box(status = "primary",
                                                      width = 12,
                                                      title = "Datasets in Grouped Analysis ",
                                                      solidHeader = TRUE,
                        #                               p("All results are plotted as median transcripts per million (TPM).
                        # Search result data for each dataset is available for download
                        # on their respective pages (see: Datasets). Rodent raw data is available
                        # on GEO. Human data from the Neural Injury Group (Oxford) is available on request.
                        # Differentially expressed genes (DEGs) are defined as an FDR < 0.05 and an absolute
                        # log2 fold change (LFC) > 1. "),
                                                      hr(),
                                                      # checkboxInput("selectall", label = "Include all datasets. ", value = TRUE),
                                                      DT::dataTableOutput("meta_table"),
                        p("Find more information about papers at: ", actionLink("citation", "Citations")),
                                                      # p("All datasets are included in default. Select multiple datasets for Grouped Analysis  by 
                                                      #      clicking on rows. Double click to deselect a dataset."),
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
                                                      HTML("<p style = 'font-size:12px'>*Ray et al. used quantile normalized TPM counts</p>"),
                                                      plothomescdot_ui("homespat")
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
                                  shinydashboard::box(title = "Pain Network Analysis",
                                                      width = 12,
                                                      status = "primary",
                                                      solidHeader = TRUE,
                                                      includeHTML("des/network.Rhtml")
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
                                                     choices = c("Composite LFC", "Predicted Pain Score", unique(network_df$experiment_id)),
                                                     selected = "Predicted Pain Score")
                                         
                                  ),
                                  column(width = 3,br(),
                                         actionButton("submit", "Construct Network")), 
                                  style = "margin-top:20px;"
                                ), 
                                
                                fluidRow(
                                  shinydashboard::box(width = 12, status = "primary",
                                                      solidHeader = TRUE, title = "Network",
                                                      column(width = 12,
                                                             shiny::downloadButton("downloadnet", "Download"),
                                                             visNetwork::visNetworkOutput("network", height = "435px")
                                                      ), 
                                                      height = "38em"

                                  ),
                                ), br(),
                                fluidRow(
                                  shinydashboard::box(width = 12, status = "primary",
                                                      solidHeader = TRUE, title = "Result Table",
                                                      tabsetPanel(
                                                        tabPanel( "Table",
                                                                  DT::dataTableOutput("contrast_table")
                                                        ),
                                                        tabPanel("Nodes",
                                                                 DT::dataTableOutput("protein_table")
                                                        )
                                                      ),
                                                      height = "38em"
                                  )
                                ),br(),
                                # a table containing experiments involved and the treatment groups
                                fluidRow(
                                  shinydashboard::box(status = "primary",
                                                      solidHeader = TRUE,
                                                      width = 12,
                                                      title = "Datasets",
                                                      column(width = 12,
                                                             includeHTML("des/datatable.html")
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
                                                      helpText("Select one dataset for analysis by 
                           clicking on rows. Click again to deselect a dataset."),
                                                      collapsible = TRUE)
                                ),  br(),# fluidRow closure
                                fluidRow(
                                  uiOutput("shinypages")
                                )
                                
        ), 
        
        tabItem(tabName="tabcode",
                h4("User Guide"),
                includeMarkdown('des/codedata.md'),
                p(strong("Citations")),
                includeHTML("des/table.html")
                
        ), #tabitem code
        
        tabItem(tabName="tabguide",
                h4("Contact"),
                includeMarkdown('des/userguide.md'),
                
                # Lab location map, in a column solely for aesthetic
                column(10, offset = 1,
                       leafletOutput('myMap', width = "100%", height = 350)
                )
                
        ) # tabItem help guide
        
        #################AMB
        
      #tabitem code
      
      # tabItem(tabName="tabtest",
      #         h4("Testing"),
      #         #SselectInput("Member", label=h5("Choose a option"),choices=c('App1','App2')),
      #         # Lab location map, in a column solely for aesthetic
      #         htmlOutput("frame")
      #        
      #        
      #) # tabtest help guide
        
        
      ) # tabItems (all)
    )  #mainpanel dashboardbody close
  ) # dashboad page close
) #fluid page
}


# ui = secure_app(ui) #adds password protection



