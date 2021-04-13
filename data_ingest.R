#load packages
library(readr); library(Seurat); 
library(dplyr); library(stringr)

#set globals
cache <-  "G:/My Drive/data/ES_singleCell"

#bring in sc data
ES_ingest <- function(cache) {
  df <- read_tsv(file = paste0(cache, "/inducible_cell_line_data/aggregated_1964.txt"))
  return(df)
}

#create metadata df from information in column names of raw count data
extract_metadata <- function(df) {
  cells <- colnames(df)[2:length(colnames(df))]
  meta_df <- data.frame(cell_name = cells) %>% 
    mutate(cell_line = str_extract(cell_name, "^[:alpha:][:digit:]{3}"))
}

df_seurat <- CreateSeuratObject(df) # why does this take an eternity?