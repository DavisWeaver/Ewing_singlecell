library(readr); library(Seurat); 
library(dplyr); library(stringr)

#set globals
cache <-  "..."

#bring in sc data from our seed paper
ES_ingest <- function(cache) {
  exp_df <- read_tsv(file = paste0(cache, "/aggregated_1964.txt"))
  meta_df <- read_tsv(file = paste0(cache, "/sample_description.txt")) 
  
  #clean up the metadata - here we limit our cells to ES cells from cell lines that
  #are either low expressing of EWS_fli1 or high expressing. We also select MSC cell line cells and myoblasts.
  meta_df <- meta_df %>% janitor::clean_names() %>% 
    filter(group == "MSC") %>% 
    mutate(subgroup = ifelse(subgroup == "22", "EWS_high", subgroup)) %>% 
    select(-name) %>% 
    rename(cell = sample)
  exp_df <- exp_df[,c("GENE",meta_df$cell)] %>% 
    rename(gene = GENE)
  
  #only meta data cells that are present in the counts matrix
  meta_df <- meta_df %>% filter(cell %in% colnames(exp_df))
  return(list(meta_df, exp_df))
}

ES_ingest_raw <- function(cache) {
  test <- read_tsv(paste0(cache, "/GSE130025_RAW/", "GSM3730172_A472U295.mapped.counts.txt.gz"))
}

neuro_ingest <- function() {
  
  exp_df <- read_tsv("https://cells-test.gi.ucsc.edu/early-brain/exprMatrix.tsv.gz")
  meta_df <- read_tsv("https://cells-test.gi.ucsc.edu/early-brain/meta.tsv")
  meta_df <- meta_df %>% janitor::clean_names() %>% 
    filter(cell_type %in% c("Neuroepithelial")) %>% 
    mutate(group = "brain_study", 
           subgroup = cell_type) %>% 
    select(cell, group, subgroup) %>% 
    sample_n(size = 1000)
  
  exp_df <- exp_df[,c("gene",meta_df$cell)] #keep only the columns in meta_df
  
  #only meta data cells that are present in the counts matrix
  meta_df <- meta_df %>% filter(cell %in% colnames(exp_df))
  return(list(meta_df, exp_df))
  
}


ES_list <- ES_ingest(cache)
neuro_list <- neuro_ingest()

meta_df <- bind_rows(ES_list[[1]], neuro_list[[1]])
exp_df <- left_join(neuro_list[[2]],ES_list[[2]])

rm(ES_list, neuro_list)
cleaned_data <- list(meta_df, exp_df)
save(cleaned_data, file = paste0(cache, "/seurat_ready.Rda"))

#load packages
library(readr); library(Seurat); 
library(dplyr); library(stringr)
library(naniar)

##set globals
cache <-  "/Users/roop/Desktop/SYBB459/"

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

#Normalization
seur_obj <- 
  NormalizeData(seur_obj, normalization.method = "LogNormalize", 
                scale.factor = 10000)

#Visualize feature-feature relationships
plot1 <- FeatureScatter(seur_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seur_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Finding Variable Features
seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seur_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seur_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(seur_obj)
seur_obj <- ScaleData(seur_obj, features = all.genes)

#Perform linear dimensional reduction (PCA) and visualize using DimHeatmap
seur_obj <- RunPCA(seur_obj, features = VariableFeatures(object = seur_obj))
DimHeatmap(seur_obj, dims = 1:15, cells = 500, balanced = TRUE)
DimPlot(seur_obj, reduction = "pca")

#Determine dimensionality of the dataset
seur_obj <- JackStraw(seur_obj, num.replicate = 100)
seur_obj <- ScoreJackStraw(seur_obj, dims = 1:20)
JackStrawPlot(seur_obj, dims = 1:15)
ElbowPlot(seur_obj)

#Cluster Cells
seur_obj <- FindNeighbors(seur_obj, dims = 1:10)
seur_obj <- FindClusters(seur_obj, resolution = 0.015)
head(Idents(seur_obj), 5)

#Run non-linear dimensionality reduction
seur_obj <- RunUMAP(seur_obj, dims = 1:10)
DimPlot(seur_obj, reduction = "umap")

#Finding differentially expressed features (cluster biomarkers)
seur_obj.markers <- FindAllMarkers(seur_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- seur_obj.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(seur_obj, features = top5$gene) + NoLegend()

#Featureplot with neural markers
FeaturePlot(seur_obj, features = c("CCND1", "LHX2", "PAX6", "SOX2", "GABA", "GAP43"))

#Featureplot with Positive MSC markers
FeaturePlot(seur_obj, features = c("CD73", "CD90", "CD105", "CD29", "CD44", "CD54", "CD106", "CD166", "CD349", "STRO-1", "TNAP"))

#Featureplot with Negative MSC markers
FeaturePlot(seur_obj, features = c("CD14", "CD34", "CD45", "CD19", "HLA-DR"))

VlnPlot(seur_obj, features = c("LHX2"))
VlnPlot(seur_obj, features = c("CD44"))

#Assigning cell type identity to clusters and plotting UMAP
new.cluster.ids <- c("EWS", "NEC", "MSC")
names(new.cluster.ids) <- levels(seur_obj)
seur_obj <- RenameIdents(seur_obj, new.cluster.ids)
DimPlot(seur_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(seur_obj, file = "FinalPlot.rds")
