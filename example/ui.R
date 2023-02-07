##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

#rsconnect::deployApp()
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.16/bioc"))

library(cowplot)
library(shiny)
library(data.table)
library(shinythemes)
library(shinydashboard)
library(shinyFeedback)
library(plotly)
library(tidyverse)
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library(ggrepel)
library(stringr)
library(ggplot2)
library(viridis)
library("pacman")
library("leaflet")
library(gridExtra)
library("cowplot")
library("plyr")
library(shinymanager)
library("fontawesome")
library(Seurat)
library("hover")
# library(omics)
library(profvis)
library(dict)

# source("shinyui.R")

shinytabUI <- function(id, tabName, datasetTitle, includedeg = FALSE, volc_title = NULL,pop_list = NULL,
                       include_subtype = FALSE, des_dir = NULL, image_dir = NULL, includegoi = NULL) {
  tabItem(tabName=tabName,
          fluidRow(
            box(title = datasetTitle,
                width = 12,
                status = "primary",
                solidHeader = TRUE,
                if (is.null(des_dir) == FALSE) {includeHTML(des_dir)},
                if (is.null(image_dir) == FALSE) {img(src = image_dir, height = 150, width = 400)}
            )
          ),
          br(),
          actionLink(NS(id, "link_to_home"), "Home", icon = icon("home")),
          br(),br(),
          actionButton(NS(id, "load"), "Plot Graphs"),
          downloadButton(NS(id,"report"), "Generate Code"),
          downloadButton(NS(id,"plots"), "Download Data"),
          br(),br(),
          shinydotUI(id),
          br(),
          if (include_subtype == TRUE) {
            shinysubtypeUI(id)
          },br(), 
          if (includedeg == TRUE) {
            shinydegplotUI(id, pop_list, volc_title)
          },br(),
          if(is.null(includegoi) == TRUE) {
            shinygoitableUI(id)
          },br(),
          if(includedeg == TRUE) {
            shinycontrasttabUI(id, pop_list)
          }
  )
}

shinydotUI <- function(id) {
  fluidRow(
    box(title = "Naive",
        status = "primary",
        solidHeader = TRUE,
        plotlyOutput(NS(id, "bulkseq_dots")),
        height = "36em"
    ),
    box(title = "Injury", status = "primary", solidHeader = TRUE, 
        plotOutput(NS(id, "bulkseq_lines")), 
        height = "36em"
    ), br(), br()
  )
}

shinyscrnaUI <- function(id, tabName, datasetTitle, des_dir = NULL, image_dir = NULL) {
  tabItem(tabName = tabName,
          fluidRow(box(title = datasetTitle,
                       width = 12,
                       status = "primary",
                       solidHeader = TRUE,
                       if (is.null(des_dir) == FALSE) {includeHTML(des_dir)},
                       if (is.null(image_dir) == FALSE) {img(src = image_dir, height = 150, width = 400)}
          )),br(),
          actionLink(NS(id, "link_to_home"), "Home", icon = icon("home")),br(),br(),
          actionButton(NS(id, "load"), "Plot Graphs"),
          downloadButton(NS(id,"scrna_plots"), "Download Plots"),
          br(), br(),
          fluidRow(
            box(width=6,
                solidHeader = TRUE,
                title = "Naive", status = "primary",
                plotOutput(NS(id, "scrna_dots"))
            ),
            box(width=6,
                solidHeader = TRUE,
                title = "UMAP", status = "primary",
                fluidRow(
                  column(12,
                         plotOutput(NS(id, "scrna_umap")))
                )
            )
          ), br(),br(),
          fluidRow(
            box(width=6,
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
                plotlyOutput(NS(id, "scrna_feature"))
            )
          ),
          fluidRow(
            plotOutput(NS(id, "scrna_violin"))
          )
  )
}


shinysubtypeUI <- function(id) {
  fluidRow(
    box(width = 12,
        title = "Subtype Results", status = "primary",
        solidHeader = TRUE, 
        plotOutput(NS(id, "bulkseq_lines_subtype"))
    )
  )
}

shinydegplotUI <- function(id, pop_list, volc_title) {
  fluidRow(
    box(width = 6,
        title = "Differential Gene Analysis",
        status = "primary",
        solidHeader = TRUE, 
        plotlyOutput(NS(id, "deg_plot")), 
        height = "38em"),
    box(width = 6,
        title = volc_title, status = "primary",
        solidHeader = TRUE, 
        actionButton(NS(id, "plotvolc"), "Plot Volcano Graphs"),
        selectInput(NS(id, "volc_pop"), "",
                    choices = pop_list,
                    selected = ""),
        plotOutput(NS(id, "volcanoplot"), height = "26em"), height = "38em"
    )
    
  )
}

shinycontrasttabUI<- function(id, pop_list) {
  fluidRow(
    box(width = 12,
        title = "Differential Analysis Table",
        status = "primary",
        solidHeader = TRUE, 
        selectInput(NS(id, "contrast"), "",
                    choices = pop_list,
                    selected = ""),
        DT::dataTableOutput(NS(id, "contrast_table"))
    )
  )
}

shinygoitableUI <- function(id) {
  fluidRow(
    box(
      width = 12,
      title = "Result Table",
      status = "primary",
      solidHeader = TRUE, 
      DT::dataTableOutput(NS(id,"goi_table"))
    )
  )
}


plothomescdot_ui <- function(id, dataset, combined=FALSE) {
  fluidRow(
    column(8, offset = 2,
           plotOutput(NS(id, "home_scrna_dots")))
  )
}

plotcombine_ui <- function(id) {
  fluidRow(
    column(12,
           plotOutput(NS(id,"dot"))
    )
  )
}

deg_combine_ui <- function(id) {
  fluidRow(
    column(12,
           plotOutput(NS(id,"deg")))
  )
}

ui = function(req) {fluidPage(
  
  includeCSS("www/style.css"),
  
  shinyFeedback::useShinyFeedback(),
  
  dashboardPage(
    dashboardHeader(title="DRG Directory", titleWidth = 225
    ),
    dashboardSidebar(width = 225,
                     sidebarMenu(
                       id = "tabs",
                       
                       HTML(paste0( # oxford logo + ndcn link
                         "<br><br>",
                         "<a href='https://www.ndcn.ox.ac.uk/research/neural-injury-group' target='_blank'> <img style = 'display: block; margin-left: auto; margin-right: auto;' src='oxfordlogo2.png' width = '105'></a>",
                         "<br><br>")
                       ),
                       
                       selectizeInput(
                         inputId = "geneid",
                         label = "Search Genes:",
                         multiple = TRUE,
                         choices = NULL
                       ),
                       
                       fileInput("file",
                                 label = "Upload File:",
                                 accept = c(
                                   'text/csv',
                                   'text/comma-separated-values',
                                   'text/tab-separated-values',
                                   'text/plain',
                                   '.csv', '.txt',
                                   '.tsv')),
                       actionButton("reset", "Clear"),
                       
                       selectizeInput(
                         inputId = "sex",
                         label = "Select Sex:",
                         choices = c('Both', 'Separate'),
                         selected = 'Both'),
                       
                       # sidebar menu for tabs (pages)
                       menuItem("Home", tabName = "tabhome", icon = icon("home")),
                       menuItem("Datasets", tabName = "tabdatasets", icon = icon("database"), 
                                menuItem("Mouse DRG Subtype (Barry)", tabName = "tabsubtype"),
                                menuItem("Mouse DRG Subtype (Zheng)", tabName = "tabzheng"),
                                menuItem("Mouse DRG Bulk", tabName = "tabmouse"), 
                                menuItem("Rat DRG Bulk", tabName = "tabrat"), 
                                menuItem("Human iPSC HSN1", tabName = "tabhuman"),
                                menuItem("Human skin Diabetic NP", tabName = "tabdb"),
                                menuItem("Human skin CTS", tabName = "tabcts"), 
                                menuItem("Human DRG spatial-seq", tabName = "tabspat"),
                                menuItem("Human DRG Bulk", tabName = "tabhdrg")
                       ),
                       menuItem("User Guide", tabName = "tabcode", icon = icon("folder-open")),
                       menuItem("Contact", tabName = "tabguide", icon = icon("info-circle")),
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
                       "<td style='padding: 5px;'><a href='https://orcid.org/0000-0002-6787-6889' target='_blank'><i class='fab fa-orcid fa-lg'></i></a></td>",
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
    dashboardBody(
      tabItems(
        tabItem(tabName="tabhome", 
                #h4("Home"),
                fluidRow(
                  tableOutput("contents"),
                  # text summary for page
                  box(width=12,
                      status = "primary", 
                      solidHeader = TRUE,
                      title = "Overview",
                      br(),
                      p("This database provides an interface to explore various -omics datasets.
                        Full descriptions and project-specific analyses are highlighed in the 'Datasets' pages. We 
                        strongly encourage users to familiarise themselves with the methodology of each study 
                        individually. Citations are provided below."),
                      br(),
                      p("All results are plotted as median transcripts per million (TPM).
                        Search result data for each dataset is available for download
                        on their respective pages (see: Datasets). Rodent raw data is available
                        on GEO. Human data from the Neural Injury Group (Oxford) is available on request.
                        Differentially expressed genes (DEGs) are defined as an FDR < 0.05 and an absolute 
                        log2 fold change (LFC) > 1. ")
                  )
                ), #fluidrow
                br(),
                fluidRow(
                  box(status = "primary",
                      solidHeader = TRUE,
                      collapsed = TRUE,
                      width = 12, 
                      title = "Dataset Overview",
                      column(width = 12, 
                             includeHTML("des/table.html"), 
                             HTML("<p style = 'font-size:12px'>*Browse individual datasets through buttons under 'Plots'</p>")
                             # p("*Browse individual datasets through buttons under 'Plots'")
                      ),
                      collapsible = TRUE),
                ),
                fluidRow(
                  br(),
                  column(12,offset = 0, 
                         actionButton("load", "Plot Graphs", class = "btn-primary-1"), 
                         downloadButton("combineplot", "Generate code"),
                         downloadButton("combineplots", "Download Plots"),
                         br()
                  ), 
                  br()
                ),
                
                fluidRow(
                  br(),
                  column(width = 12,
                         p("Gene search can be modified in the side bar. ")
                  ),
                  box(width = 12,
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
                  box(width = 12, 
                      title = "Differential Gene Analysis",
                      collapsible = TRUE,   
                      solidHeader = TRUE,
                      status = "primary", 
                      deg_combine_ui("deg_plot"))
                )
                
                
        ), # tabItem HOME
        
        # tab subtype
        shinytabUI("subtypetab", "tabsubtype", "Subtype DRG (Barry)", TRUE, "Ipsilateral vs Contralateral",
                   pop_list = c("Nociceptors 3D","Nociceptors 4W","PEP 3D","PEP 4W","NP 3D",
                                "NP 4W", "C-LTMR 3D","C-LTMR 4W","Ad- AB-RA LTMRs 3D",
                                "Ad- AB-RA LTMRs 4W"), include_subtype = TRUE, des_dir = "des/datasetsummary.Rhtml",
                   image_dir = "schematic.png"),
        
        # for mouse data page
        shinytabUI("mousetab", "tabmouse", "Mouse DRG Bulk", TRUE, "SHAM vs SNI",
                   pop_list = c("B10D2", "BALB"), des_dir = "des/mousedrg.Rhtml"),
        
        shinytabUI("rattab", "tabrat", "Rat DRG Bulk", TRUE, "SHAM vs SNT",
                   pop_list = c("rat"), des_dir = "des/mousedrg.Rhtml"),
        
        # need two: by differentiation or by cell line
        shinytabUI("ipsctab", "tabhuman", "Human iPSC HSN1", TRUE, "Healthy vs HSN1",
                   pop_list = c("iPSCDN_young","iPSCDN_old"), des_dir = "des/ipsc_notes.Rhtml", image_dir = NULL),
        
        shinytabUI("dbtab", "tabdb", "Human skin Diabetes", TRUE, "Painful vs Painless",
                   pop_list= c("Diabetes_skin","Diabetes_skin_females",
                               "Diabetes_skin_male"), des_dir = "des/diabetes.Rhtml", includegoi = FALSE),
        
        shinytabUI("ctstab", "tabcts", "Human skin Carpal Tunnel (CTS)", TRUE, "Painful vs Painless",
                   pop_list = c("Skin HS"), des_dir = "des/humanskincts.Rhtml", includegoi = FALSE),
        
        shinyscrnaUI("spattab", "tabspat", "Human DRG spatial-seq", des_dir="des/spatial.Rhtml"),
        shinytabUI("zhengtab", "tabzheng", "Mouse DRG Subtype (Zheng)", des_dir = "des/zheng.Rhtml"),
        shinytabUI("hdrgtab", "tabhdrg", "Human DRG RNAseq", des_dir = "des/humandrg.Rhtml"),
        
        ## Supply simple links for each paper + supplementary repository
        tabItem(tabName="tabcode",
                h4("Data Access"),
                includeMarkdown('des/codedata.md'),
                br(),
                includeHTML("des/table.html"),
                br(),
                br(),
                
                # HTML links to funders, embedded in a table for alignment
                HTML(paste0(
                  "<table style='margin-left:auto; margin-right:auto; bottom=0;'>",
                  
                  "<tr>",
                  "<td><a href='https://www.gtc.ox.ac.uk/' target='_blank'> <img style= 'display: center;' src='gtclogo.png', width = '100'></a></td>",
                  "<td><a href='https://wellcome.org/' target='_blank'> <img style = 'display: center;' src='wt-logo.svg'. width = '100'></a></td>",
                  "<td><a href='https://www.ukri.org/councils/mrc/' target='_blank'> <img style = 'display: center;' src='ukri.png'. width = '90' height='90'></a></td>",
                  "</tr>",
                  "</table>"
                )
                )
                
        ), #tabitem code
        
        tabItem(tabName="tabguide",
                h4("Contact"),
                includeMarkdown('des/userguide.md'),
                
                # Lab location map, in a column solely for aesthetic
                column(10, offset = 1,
                       leafletOutput('myMap', width = "100%", height = 350)
                )
                
        ) # tabItem help guide
      ) # tabItems (all)
    )  #mainpanel dashboardbody close
  ) # dashboad page close
) #fluid page
}#ui

# ui = secure_app(ui)
