#load packages
library(readr); library(Seurat); 
library(dplyr); library(stringr)
library(naniar)

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
  #remove all non-numbers from cell names
  colnames(counts) <- str_remove_all(colnames(counts), "[[:alpha:]]+")
  colnames(counts) <- str_remove_all(colnames(counts), "_")
  #assume all missing data had zero read counts
  counts[is.na(counts)] <- 0
  return(counts)
}

clean_meta <- function() {
  meta <- cleaned_data[[1]]
  meta <- meta %>% 
    mutate(cell = str_remove_all(cell, "[[:alpha:]]+"), 
           cell = str_remove_all(cell, "_"))
  return(meta)
}
counts <- clean_counts()
meta <- clean_meta()


## THis takes a while 
seur_obj <- CreateSeuratObject(counts = counts, project = "ES_origin", 
                               assay = "RNA",
                               meta.data = meta)
##add mitochondrial gene pct
seur_obj[["percent.mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seur_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seur_obj <- 
  NormalizeData(seur_obj, normalization.method = "LogNormalize", 
                scale.factor = 10000)
