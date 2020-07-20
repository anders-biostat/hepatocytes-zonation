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
    list("seurat_clusters",    "cluster"),
    list("Experimental.Batch", "batches"),
    list("Mouse.ID",           "mouse"))) {
    q <- UMAPPlot(seus[[tissue]], group.by = grouping[[1]])
    ggsave(figpattern("umap-{tissue}-{group}.png",
      tissue = tissue, group = grouping[[2]]), q, width = 6, height = 6)
  }
}

markers <- loadConfig()$markers
## all genes are found in count matrix
setdiff(unlist(markers),  rownames(counts))

## plot expression of marker gene for all and for every sorting
purrr::iwalk(list(
  "all" = seu,
  "HEP" = seus$HEP,
  "LSEC" = seus$LSEC), function(seuobj, celltype) {
    q <- FeaturePlot(seuobj, unlist(markers))
    ggsave(figpattern("umap-markers-{name}.png", name = celltype),
      q, width = 13, height = 13, dpi = 200)
  })

q <- UMAPPlot(seu, group.by = "Mouse.ID", split.by = "Cell.type")
ggsave(figpattern("umap-all-by-celltype.png"), q, width = 10, height = 5, dpi = 200)

s <- FindAllMarkers(seus$LSEC)

## plot how samples contribute to clusters and how sampels are spread among clusters
with(seu@meta.data,
  table(seurat_clusters, sample) %>% prop.table(2) %>%
    pheatmap(filename = figpattern("all-prop-of-clusters-in-samples.png"),
      width = 8, height = 8, units = "in"))
dev.off()

with(seu@meta.data,
  table(seurat_clusters, sample) %>% prop.table(1) %>%
    pheatmap(filename = figpattern("all-prop-of-samples-in-clusters.png"),
      width = 8, height = 8, units = "in"))
dev.off()

saveRDS(list(
  all = seu,
  HEP = seus$HEP,
  LSEC = seus$LSEC
), file = file.path(rdsDir, "seurats-qc.rds"))
