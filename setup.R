##########################
# Info: You can store datafiles in two ways: either through a R object  (a list of lists)  or a SQLite database. 
# To create these, you need to store your count data as .csv format in data/CountTables/ folder 
# you also need to store your ColData in data/ColData/ folder 
# you will need to store your DegData in data/DegData/ folder 
# Then, run code in chunks to complete the creation of a R object or the database. 
######################

#1. Install necessary packages
{

  #install.packages("stringr")
  library(stringr)
  library(RSQLite)
  library(dplyr)
  library(tidyr)

}

################################ Option 1: Create a sql database ############################################################
# create the database 

# 1. create a database by entering a directory in ".db" format. e.g. "drg.db" 

db_dir = "test.db" 
{
  conn <- dbConnect(RSQLite::SQLite(), db_dir)
  h = read.csv("./data/GeneData/human_gene_data.csv", header = TRUE)
  m = read.csv("./data/GeneData/mouse_gene_data.csv", header = TRUE)
  r = read.csv("./data/GeneData/rat_gene_data.csv", header = TRUE)
  genes = read.csv("./data/GeneData/genebank.csv", header = TRUE)
  dbWriteTable(conn, "human_gene_data", h)
  dbWriteTable(conn, "rat_gene_data", r)
  dbWriteTable(conn, "mouse_gene_data", m)
  dbWriteTable(conn, "genebank", genes)
  
}

# 2. read counttables and colDatas into the database 
{
  Counts <- str_sort(list.files(path = "./data/CountTables", pattern = "*.csv")) # names of the count tables

  names = as.vector(substr(Counts, 0, nchar(Counts) - 4))
  
  for (i in 1:length(names)){
    name = names[i]
    value = read.csv(paste0("./data/CountTables/", Counts[i]), header = TRUE)
    dbWriteTable(conn, name, value)
  }

  CDs <- str_sort(list.files(path = "./data/ColDatas", pattern = "*.csv")) # names of colDatas 
  
  names <- substr(CDs, 0, nchar(CDs) - 4)# assign name 
  
  # read in colData
  for (i in 1:length(names)) {
    name = names[i]
    value = read.csv(paste0("./data/ColDatas/", CDs[i]), header = TRUE, row.names = 1)
    
    dbWriteTable(conn, name, value, row.names = TRUE, overwrite = TRUE)
  }
}

deg_df_names = c("subtype_deg_df", "mouse_deg_df", "rat_deg_df", "ipsc_deg_df", "db_deg_df", "cts_deg_df")

# 3. read deg_dfs from each study into the database; generate a combined deg_df for each species  
{
  # write deg_dfs into db 
 
  # you need to specify the names of the dataset (experiments) and store them in the list below
  Dataset = c("mouse/Subtype DRG (Barry)/", "mouse/Mouse DRG Bulk/", "rat/Rat DRG Bulk/",
              "human/Human iPSC HSN1/", "human/Human skin Diabetes/", "human/Human skin Carpal Tunnel (CTS)/")
  
  human_deg_df = data.frame()
  mouse_deg_df = data.frame()
  rat_deg_df = data.frame()
  
  # for each experiment
  for (j in 1:length(Dataset)) {
    name = Dataset[j] # name of experiment
    species = unlist(strsplit(name, split = "/"))[[1]]
    
    # DEs are a list of datasets in this experiment 
    DEs = str_sort(list.files(path = paste("./data/DegData/", name, sep=""), pattern = "*.csv"))
    DE_names <- substr(DEs, 0, nchar(DEs) - 4)# assign name to each population
    deg_list = list()
    
    # read in those datasets into one df for each experiment
    for (i in 1:length(DE_names)) {
      df_name = DE_names[i]
      value = read.csv(paste("./data/DegData/", name, DEs[i], sep=""), header = TRUE, row.names = 1)
      # dbWriteTable(conn, df_name, value) 
      deg_list = append(deg_list, list(value))
    }
    
    dataset_combine_deg = generate_degdf(deg_list, DE_names, unlist(strsplit(name, split = "/"))[[2]])
    dbWriteTable(conn, deg_df_names[j], dataset_combine_deg)
    # add to species df
    if (species == "human"){
      human_deg_df = rbind(human_deg_df, dataset_combine_deg)
    }
    if (species == "mouse"){
      mouse_deg_df = rbind(mouse_deg_df, dataset_combine_deg)
    }
    if (species == "rat"){
      rat_deg_df = rbind(rat_deg_df, dataset_combine_deg)
    }
    
  }
  
  # write species deg_df into table
  dbWriteTable(conn, "mouse_all_deg_df", mouse_deg_df)
  dbWriteTable(conn, "human_deg_df", human_deg_df)
  dbWriteTable(conn, "rat_deg_df", rat_deg_df)

}

# 3. generate combined deg_df for each species, using combined df above 
{

  #' This function generates a combined deg_df for multiple studies 
  #' @param deg_df_list a list object containing the list of deg_dfs 
  #' @param pop_list a list of names of those dfs 
  #' @param species the species of the list 
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
 
}
  
 
{
  # if you have any additional data that you want to write into the table, use the code below 
  Datas <- str_sort(list.files(path = "./data/AdditionalData", pattern = "*.csv")) 
  
  names = as.vector(substr(Datas, 0, nchar(Datas) - 4))
  
  for (i in 1:length(names)){
    name = names[i]
    value = read.csv(paste0("./data/AdditionalData/", Datas[i]), header = TRUE)
    value = value[, !names(value) %in% c("description")]
    dbWriteTable(conn, name, value, overwrite = TRUE)
  }

}

# Now, you can use the dir of db you created to run the app 
# To do this, you will need to change the dir of the .db at the top to your db_dir, then run the app. 


datasets = c("HIV", "SNI_adult", "SNT_L5vsSHAM")
for (exp in datasets) {
  df = allde[allde$experiment_id == exp,]
  df = df[df$gene %in% rat_gene_data$ensembl_gene_id,]
  df$symbol = rat_gene_data[rat_gene_data$ensembl_gene_id %in% df$gene,]$rgd_symbol
  dbWriteTable(conn, paste(exp, "_deg_df", sep = ""), df, overwrite = TRUE)
}

df = allde[allde$experiment_id == "bone_cancer",]
df = df[df$gene %in% mouse_gene_data$ensembl_gene_id,]
df$symbol = mouse_gene_data[mouse_gene_data$ensembl_gene_id %in% df$gene,]$mgi_symbol
dbWriteTable(conn, paste("bone_cancer", "_deg_df", sep = ""), df, overwrite = TRUE)


######################## Option 2: Create R object #############################################################

#2. Create Seq_Data object and load in essential ones 
{
  Seq_Data <- vector(mode = "list", 3)
  names(Seq_Data) = c("CountTable", "ColData", "DegData")
}

# 3. load in count data 
{
  Counts <- str_sort(list.files(path = "./data/CountTables", pattern = "*.csv")) # names of the count table
  # add them into the object 
  Seq_Data[[1]] <- vector(mode = "list", length(Counts))
  names(Seq_Data[[1]]) <- substr(Counts, 0, nchar(Counts) - 4) # assign name 
  
  #Load in data
  for (i in 1:length(Seq_Data[[1]])) {
    Seq_Data[[1]][[i]] = read.csv(paste0("./data/CountTables/", Counts[i]), header = TRUE, row.names = 1)
  }
}

# 4. load in colData
{
  CDs <- str_sort(list.files(path = "./data/ColDatas", pattern = "*.csv")) # names of colDatas 
  
  Seq_Data[[2]] <- vector(mode = "list", length(CDs))
  names(Seq_Data[[2]]) <- substr(CDs, 0, nchar(CDs) - 4)# assign name 
  
  
  # read in colData
  for (i in 1:length(Seq_Data[[2]])) {
    Seq_Data[[2]][[i]] = read.csv(paste0("./data/ColDatas/", CDs[i]), header = TRUE, row.names = 1)
  }
  
}

#5. load in DegData 

{
  Seq_Data[[3]] <- vector(mode = "list", 3)
  names(Seq_Data[[3]]) = c("HumanDEG", "MouseDEG", "RatDEG")
  
  
  # first, read in human deg data
  Seq_Data[[3]][[1]] <- vector(mode = "list", length(hDEs))
  hDEs = str_sort(list.files(path = "./data/DegData/human", pattern = "*.csv"))
  names <- substr(hDEs, 0, nchar(hDEs) - 4)# assign name 
  names(Seq_Data[[3]][[1]]) = substr(hDEs, 0, nchar(hDEs) - 4)# assign name 
  
  for (i in 1:length(names)) {
    name = names[i]
    value = read.csv(paste0("./data/DegData/human/", hDEs[i]), header = TRUE, row.names = 1)
    Seq_Data[[3]][[1]][[i]] = value
  }
  
  
  # mouse deg_dfs
  mDEs = str_sort(list.files(path = "./data/DegData/mouse", pattern = "*.csv"))
  names <- substr(mDEs, 0, nchar(mDEs) - 4)# assign name 
  Seq_Data[[3]][[2]] <- vector(mode = "list", length(mDEs))
  names(Seq_Data[[3]][[2]]) <- substr(mDEs, 0, nchar(mDEs) - 4)# assign name 
  
  # read in mousedeg data into table
  for (i in 1:length(names)) {
    name = names[i]
    value = read.csv(paste0("./data/DegData/mouse/", mDEs[i]), header = TRUE, row.names = 1)
    Seq_Data[[3]][[2]][[i]] = value
  }
  
  # rat
  rDEs = str_sort(list.files(path = "./data/DegData/rat", pattern = "*.csv"))
  names <- substr(rDEs, 0, nchar(rDEs) - 4)# assign name 
  Seq_Data[[3]][[3]] <- vector(mode = "list", length(rDEs))
  names(Seq_Data[[3]][[3]]) <- substr(rDEs, 0, nchar(rDEs) - 4)# assign name 
  
  
  # read in human deg data into table
  for (i in 1:length(names)) {
    name = names[i]
    value = read.csv(paste0("./data/DegData/rat/", rDEs[i]), header = TRUE, row.names = 1)
    Seq_Data[[3]][[3]][[i]] = value
  }
  
}


#  6. create and write combined datasets mouse_all_deg_df, human_deg_df, rat_deg_df
{
  
  # generate combined deg_dfs 
  generate_degdf <- function(deg_df_list, pop_list, species) {
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
    if (!is.null(species)) {
      final_df$Species = rep(species, nrow(final_df))
    }
    final_df = mutate(final_df, sig=ifelse(final_df$padj<0.05, "SIG", "NS"))
    return(final_df)
  }
  
  # for all ones 
  human_deg_df = generate_degdf(Seq_Data[[3]][[1]], pop_list = names(Seq_Data[[3]][[1]]), species = "human")
  mouse_deg_df = generate_degdf(Seq_Data[[3]][[2]], pop_list = names(Seq_Data[[3]][[2]]), species = "mouse")
  rat_deg_df = generate_degdf(Seq_Data[[3]][[3]], pop_list = names(Seq_Data[[3]][[3]]), species = "rat")
  
  
  # deg_df for individual datasets? 
  
  # then, combined into combined
  mouse_deg_df = generate_degdf(list(b10d2, balb), c("b10d2", "balb"), "Mouse DRG")
  rat_deg_df = generate_degdf(list(ratdrg), c("DRG"), "Rat DRG")
  ipsc_deg_df = generate_degdf(list(young, old), c("young", "old"), "iPSC HSN1")
  db_deg_df = generate_degdf(list(DB, DB_female, DB_male), c("Diabetes", "Diabetes_male","Diabetes_female"), "DPN (Skin)")
  cts_deg_df = generate_degdf(list(HS), c("skin"), "CTS (Skin)")
  
  
  
}

