## In this script we filter out cells which are not of target cell type,
## low quality and so on.
source("code/func.r")
source("code/assets.r")
library(Seurat)

counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

cellanno$sample <- paste(cellanno$condition, cellanno$Mouse.ID, cellanno$Experimental.Batch)

seu <- Seurat::CreateSeuratObject(counts, names.delim = ",")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu)
seu <- Seurat::ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:50)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:50)
seu <- FindClusters(seu, resolution = 1)

rownames(cellanno) <- colnames(counts)
seu <- AddMetaData(seu, cellanno)
q <- UMAPPlot(seu, group.by = "condition")
ggsave(file.path("results", "figs", "umap-all-samples.png"), q, width = 6, height = 6)


q <- UMAPPlot(seu, group.by = "Experimental.Batch")
ggsave(file.path("results", "figs", "umap-all-samples-batches.png"), q, width = 6, height = 6)
q <- UMAPPlot(seu, group.by = "Mouse.ID")
ggsave(file.path("results", "figs", "umap-all-samples-mouse.png"), q, width = 6, height = 6)


marks <- c(
  "pecam1",
  "apoa1",
  "glul",
  "acly",
  "asl",
  "cyp2e1",
  "cyp2f2",
  "ass1",
  "alb",
  "mup3",
  "pck1",
  "g6pc"
)
FeaturePlot(seu, marks)
