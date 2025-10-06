library(Seurat)
library(dplyr)
adata=readRDS(file='~/Merged_mOP.rds')

adata[["RNA"]] <- split(adata[["RNA"]], f = adata$species)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata)
adata <- ScaleData(adata)
adata <- RunPCA(adata,npcs=150)

adata <- IntegrateLayers(
  object = adata, method =RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
verbose = TRUE,group.by = 'species'
)
saveRDS(adata, file='~/rpcaMergeMOp.rds')

adata <- FindNeighbors(adata, reduction = "integrated.rpca", dims = 1:150)
adata <- FindClusters(adata, resolution = 5, algorithm=1,cluster.name = "rpca_clusters")
adata <- RunUMAP(adata, reduction = "integrated.rpca", dims = 1:150, reduction.name = "umap.rpca")
saveRDS(adata, file='~/rpcaMergeMOp.rds')
p1 <- DimPlot(
  adata,
  reduction = "umap.rpca",
  group.by = c("species", "rpca_clusters"),
  combine = FALSE,raster=FALSE
)
p1
adata <- JoinLayers(adata)
all.markers <- FindAllMarkers(adata, min.pct = 0.15)
write.csv(all.markers, "~/MOpRPCAallmarkers.csv", row.names=TRUE)
saveRDS(adata, file='~/rpcaMergeMOp.rds')
