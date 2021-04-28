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
    mutate(cell_line = str_extract(cell_name, "^[:alpha:][:digit:]{3}")) #This is regex to grab cell line names
}

df_seurat <- CreateSeuratObject(df) # why does this take an eternity?

#Pre-Processing Seurat Workflow
df_seurat[["percent.mt"]] <- PercentageFeatureSet(df_seurat, pattern = "^MT-") #Creating a percent.mt column for showing percentage of mitochondrial genes
VlnPlot(df_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #Visualize QC metrics as a violin plot
df_seurat <- subset(df_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) #subsetting the dataset and removing bad quality cells

#Normalizing the dataset
df_seurat <- NormalizeData(df_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
df_seurat <- FindVariableFeatures(df_seurat, selection.method = "vst", nfeatures = 2000)

#Integration of Datasets
df_seurat.anchors <- FindIntegrationAnchors(object.list = df_seurat.list, anchor.features = features)

#Scaling the data
all.genes <- rownames(df_seurat)
df_seurat <- ScaleData(df_seurat, features = all.genes)

#Perform linear dimensional reduction (PCA) and visualize using DimHeatmap
df_seurat <- RunPCA(df_seurat, features = VariableFeatures(object = df_seurat))
DimHeatmap(df_seurat, dims = 1:15, cells = 500, balanced = TRUE)

#Cluster Cells
df_seurat <- FindNeighbors(df_seurat, dims = 1:10)
df_seurat <- FindClusters(df_seurat, resolution = 0.5)

#Run non-linear dimensionality reduction
df_seurat <- RunUMAP(df_seurat, dims = 1:10)

#Finding differentially expressed features (cluster biomarkers)
df_seurat.markers <- FindAllMarkers(df_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Assigning cell type identity to clusters and plotting UMAP
new.cluster.ids <- c("EWS", "MSC", "NEC")
names(new.cluster.ids) <- levels(df_seurat)
df_seurat <- RenameIdents(df_seurat, new.cluster.ids)
DimPlot(df_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(df_seurat, file = "X.rds")

