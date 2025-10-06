library(Seurat)
library(dplyr)
adata=readRDS(file='~/HvQvM900k.rds')

adata[["RNA"]] <- split(adata[["RNA"]], f = adata$species)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata)
adata <- ScaleData(adata)
adata <- RunPCA(adata)

adata <- IntegrateLayers(
  object = adata, method=FastMNNIntegration,
  orig.reduction = "pca", new.reduction = "integrated.mnn",
verbose = TRUE,group.by = 'species'
)

adata <- FindNeighbors(adata, reduction = "integrated.mnn", dims = 1:30)
adata <- FindClusters(adata, resolution = 5, algorithm=1,cluster.name = "mnn_clusters")
adata <- RunUMAP(adata, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

p1 <- DimPlot(
  adata,
  reduction = "umap.mnn",
  group.by = c("species", "mnn_clusters"),
  combine = FALSE,raster=FALSE
)
p1
adata <- JoinLayers(adata)
all.markers <- FindAllMarkers(adata, min.pct = 0.15)
write.csv(all.markers, "~/900kMNNallmarkers.csv", row.names=TRUE)
saveRDS(adata, file='~/mnnMerge900k.rds')
