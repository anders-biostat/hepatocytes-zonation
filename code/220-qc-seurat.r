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
ggsave(figpattern("umap-all-samples.png"), q, width = 6, height = 6)
q <- UMAPPlot(seu, group.by = "Experimental.Batch")
ggsave(figpattern("umap-all-samples-batches.png"), q, width = 6, height = 6)
q <- UMAPPlot(seu, group.by = "Mouse.ID")
ggsave(figpattern("umap-all-samples-mouse.png"), q, width = 6, height = 6)

## subset by cell type and create a list of seurats by tissue
tissues <- c("HEP", "LSEC")
seus <- list()
for (tissue in tissues) {
  cells <- cellanno$Cell.type == tissue
  seus[[tissue]] <- makeSeu(counts[, cells], cellanno[cells, ])
}

## plot umaps for various groupings
for (tissue in tissues) {
  for (grouping in list(
    list("condition",          "condition"),
    list("Experimental.Batch", "batches"),
    list("Mouse.ID",           "mouse"))) {
    q <- UMAPPlot(seus[[tissue]], group.by = grouping[[1]])
    ggsave(figpattern("umap-{tissue}-{group}.png",
      tissue = tissue, group = grouping[[2]]), q, width = 6, height = 6)
  }
}

## compare heps and lsec
x <- FindMarkers(seu, ident.1 = "LSEC", ident.2 = "HEP", group.by = "Cell.type")
write.csv(x, tablepattern("lsec-vs-heps-seurat-find-markers.csv"))


marks <- c(
  "pecam1",
  "stab2",
  "lyve1",
  "kdr",
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
q <- FeaturePlot(seu, marks)
ggsave(figpattern("umap-markers.png"), q, width = 12, height = 12, dpi = 100)
