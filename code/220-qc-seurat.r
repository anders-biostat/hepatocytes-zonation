## In this script we filter out cells which are not of target cell type,
## low quality and so on.
source("code/func.r")
source("code/assets.r")
library(Seurat)

counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

cellanno$sample <- paste(cellanno$condition, cellanno$Mouse.ID, cellanno$Experimental.Batch)

makeSeu <- function(counts, cellanno, dims = 1:50, res = 1) {
  seu <- Seurat::CreateSeuratObject(counts, names.delim = ",",
    min.cells = 3, min.features = 200)
  rownames(cellanno) <- colnames(counts)
  seu <- AddMetaData(seu, cellanno)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <- FindVariableFeatures(seu)
  seu <- Seurat::ScaleData(seu)
  seu <- RunPCA(seu, npcs = max(dims))
  seu <- RunUMAP(seu, reduction = "pca", dims = dims)
  seu <- FindNeighbors(seu, reduction = "pca", dims = dims)
  seu <- FindClusters(seu, resolution = res)
  seu
}

seu <- makeSeu(counts, cellanno)

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
FeaturePlot(seuHEP, marks)

heps <- cellanno$Cell.type == "HEP"
seuHEP <-  makeSeu(counts[, heps], cellanno[heps,])

x <- seuHEP@meta.data
table(x$sample, x$seurat_clusters)

q <- UMAPPlot(seuHEP, group.by = "condition")
ggsave(file.path("results", "figs", "umap-HEP-samples.png"), q, width = 6, height = 6)
q <- UMAPPlot(seuHEP, group.by = "Experimental.Batch")
ggsave(file.path("results", "figs", "umap-HEP-samples-batches.png"), q, width = 6, height = 6)
q <- UMAPPlot(seuHEP, group.by = "Mouse.ID")
ggsave(file.path("results", "figs", "umap-HEP-samples-mouse.png"), q, width = 6, height = 6)


lsec <- cellanno$Cell.type != "HEP"
seuLSEC <-  makeSeu(counts[, lsec], cellanno[lsec,])

q <- UMAPPlot(seuLSEC, group.by = "condition")
ggsave(file.path("results", "figs", "umap-LSEC-samples.png"), q, width = 6, height = 6)
q <- UMAPPlot(seuLSEC, group.by = "Experimental.Batch")
ggsave(file.path("results", "figs", "umap-LSEC-samples-batches.png"), q, width = 6, height = 6)
q <- UMAPPlot(seuLSEC, group.by = "Mouse.ID")
ggsave(file.path("results", "figs", "umap-LSEC-samples-mouse.png"), q, width = 6, height = 6)

## seuHEP <- FindNeighbors(seuHEP, dims = 1:10)
seuHEP <- FindClusters(seuHEP, resolution = 0.5)
seuHEP.markers <- FindAllMarkers(seuHEP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
