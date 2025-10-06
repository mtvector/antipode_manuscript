library(Seurat)
library(dplyr)
adata=readRDS(file='~/HvQvM900k.rds')

adata[["RNA"]] <- split(adata[["RNA"]], f = adata$species)

adata <- NormalizeData(adata)
adata <- FindVariableFeatures(adata)
adata <- ScaleData(adata)
adata <- RunPCA(adata)

adata <- IntegrateLayers(
  object = adata, method =HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harm",
verbose = TRUE,group.by = 'species'
)
saveRDS(adata, file='~/harmMerge900k.rds')

adata <- FindNeighbors(adata, reduction = "integrated.harm", dims = 1:30)
adata <- FindClusters(adata, resolution = 5, algorithm=1,cluster.name = "harm_clusters")
adata <- RunUMAP(adata, reduction = "integrated.harm", dims = 1:30, reduction.name = "umap.harm")

p1 <- DimPlot(
  adata,
  reduction = "umap.harm",
  group.by = c("species", "harm_clusters"),
  combine = FALSE,raster=FALSE
)
p1
adata <- JoinLayers(adata)
all.markers <- FindAllMarkers(adata, min.pct = 0.15)
write.csv(all.markers, "~/900kHarmallmarkers.csv", row.names=TRUE)
saveRDS(adata, file='~/harmMerge900k.rds')
