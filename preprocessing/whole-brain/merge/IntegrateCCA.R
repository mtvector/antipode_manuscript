library(Seurat)
library(dplyr)
adata=readRDS(file='~/HvQvM900k.rds')

adata[["RNA"]] <- split(adata[["RNA"]], f = adata$species)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata)
adata <- ScaleData(adata)
adata <- RunPCA(adata)

adata <- IntegrateLayers(
  object = adata, method =CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
verbose = TRUE,group.by = 'species'
)

adata <- FindNeighbors(adata, reduction = "integrated.cca", dims = 1:30)
adata <- FindClusters(adata, resolution = 5, algorithm=1,cluster.name = "cca_clusters")
adata <- RunUMAP(adata, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
saveRDS(adata, file='~/ccaMerge900k.rds')

p1 <- DimPlot(
  adata,
  reduction = "umap.cca",
  group.by = c("species", "cca_clusters"),
  combine = FALSE,raster=FALSE
)
p1
adata <- JoinLayers(adata)
all.markers <- FindAllMarkers(adata, min.pct = 0.15)
write.csv(all.markers, "~/900kCCAallmarkers.csv", row.names=TRUE)
saveRDS(adata, file='~/ccaMerge900k.rds')
