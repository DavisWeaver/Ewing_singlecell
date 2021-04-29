#load packages
library(readr); library(Seurat); 
library(dplyr); library(stringr)

#set globals
cache <-  "G:/My Drive/data/ES_singleCell"

#bring in sc data from our seed paper
ES_ingest <- function(cache) {
  exp_df <- read_tsv(file = paste0(cache, "/inducible_cell_line_data/aggregated_1964.txt"))
  meta_df <- read_tsv(file = paste0(cache, "/inducible_cell_line_data/sample_description.txt")) 
  
  #clean up the metadata - here we limit our cells to ES cells from cell lines that
  #are either low expressing of EWS_fli1 or high expressing. We also select MSC cell line cells and myoblasts.
  meta_df <- meta_df %>% janitor::clean_names() %>% 
    filter(group == "MSC" | group == "MYOBLASTS" | 
             group == "ASP14" & subgroup %in% c("7", "22")) %>% 
    mutate(subgroup = ifelse(subgroup == "7", "EWS_low", subgroup), 
           subgroup = ifelse(subgroup == "22", "EWS_high", subgroup), 
           subgroup = ifelse(group != "ASP14", group, subgroup)) %>% 
    select(-name) %>% 
    rename(cell = sample)
  exp_df <- exp_df[,c("GENE",meta_df$cell)] %>% 
    rename(gene = GENE)
  
  return(list(meta_df, exp_df))
}

neuro_ingest <- function() {
  
  exp_df <- read_tsv("https://cells-test.gi.ucsc.edu/early-brain/exprMatrix.tsv.gz")
  meta_df <- read_tsv("https://cells-test.gi.ucsc.edu/early-brain/meta.tsv")
  meta_df <- meta_df %>% janitor::clean_names() %>% 
    filter(cell_type %in% c("Neuroepithelial", "Radial Glial")) %>% 
    mutate(group = "brain_study", 
           subgroup = cell_type) %>% 
    select(cell, group, subgroup)
  
  exp_df <- exp_df[,c("gene",meta_df$cell)] #keep only the columns in meta_df
  return(list(meta_df, exp_df))
  
}


ES_list <- ES_ingest(cache)
neuro_list <- neuro_ingest()

meta_df <- bind_rows(ES_list[[1]], neuro_list[[1]])
exp_df <- left_join(ES_list[[2]], neuro_list[[2]])

rm(ES_list, neuro_list)
cleaned_data <- list(meta_df, exp_df)
save(cleaned_data, file = paste0(cache, "/seurat_ready.Rda"))

