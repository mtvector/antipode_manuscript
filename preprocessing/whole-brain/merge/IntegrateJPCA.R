library(Seurat)
library(dplyr)
adata=readRDS(file='~/HvQvM900k.rds')

adata[["RNA"]] <- split(adata[["RNA"]], f = adata$species)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata)
adata <- ScaleData(adata)
adata <- RunPCA(adata)

adata <- IntegrateLayers(
  object = adata, method =JointPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.jpca",
verbose = TRUE,group.by = 'species'
)

adata <- FindNeighbors(adata, reduction = "integrated.jpca", dims = 1:30)
adata <- FindClusters(adata, resolution = 5, algorithm=1,cluster.name = "jpca_clusters")
adata <- RunUMAP(adata, reduction = "integrated.jpca", dims = 1:30, reduction.name = "umap.jpca")
saveRDS(adata, file='~/jpcaMerge900k.rds')
p1 <- DimPlot(
  adata,
  reduction = "umap.jpca",
  group.by = c("species", "jpca_clusters"),
  combine = FALSE,raster=FALSE
)
p1
adata <- JoinLayers(adata)
all.markers <- FindAllMarkers(adata, min.pct = 0.15)
write.csv(all.markers, "~/900kJPCAallmarkers.csv", row.names=TRUE)
saveRDS(adata, file='~/jpcaMerge900k.rds')
