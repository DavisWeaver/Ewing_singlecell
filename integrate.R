#load packages
library(readr); library(Seurat); 
library(dplyr); library(stringr)

##set globals
cache <-  "G:/My Drive/data/ES_singleCell"

load(file = paste0(cache, "/seurat_ready.Rda"))

clean_counts <- function() {
  counts <- cleaned_data[[2]] %>% 
    distinct(gene, .keep_all = TRUE) %>% 
    filter(!is.na(gene)) %>% 
    as.data.frame()
  rownames(counts) <- counts$gene
  counts <- select(counts, -gene) %>% 
    as.matrix()
  return(counts)
}
counts <- clean_counts()
meta <- cleaned_data[[1]]

## THis takes a while 
seur_obj <- CreateSeuratObject(counts = counts, project = "ES_origin", 
                               assay = "RNA"
                               meta.data = meta)

seur_list <- SplitObject(seur_obj, split.by = "group")
