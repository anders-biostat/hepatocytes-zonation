## cluster 6 expresses markers of stellate cells (dcn), so we remove them.
## There are other smaller clusters, but they express stab2, kdr, cdh5, so
## I do not exclude them.
source("code/func.r")
source("code/assets.r")
library(Seurat)

figdir <- file.path("results", "figs", "qc")
counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

s <- readRDS(file.path(rdsDir, "seurats-qc.rds"))
lsec <- s$LSEC

table(lsec@meta.data$seurat_clusters)
##
## 0   1   2   3   4   5   6   7   8   9  10  11  12
## 580 397 319 309 301 254 122  55  55  20  19  18  18

## we plot two ko genes and stellate cell markers, color by cluster 6
plotKO <- function(lsec) {
  genes <- c("prelp", "ecm1", "gdf2", "rspo3", "dcn", "wls")
  totals <- colSums(lsec@assays$RNA@counts)
  expr <- lsec@assays$RNA@counts[genes,]
  i <- lsec$Genotype == "DoubleKO"
  pairs(
    jitter(as.matrix(t(expr) / totals)[i, ]),
    pch = 16, cex = .8,
    col = (d$LSEC$seurat_clusters == 6)[i] + 1)
}

png(figpattern("ko-gene-expr.png"), width = 7, height = 7, units = "in", res = 200)
plotKO(lsec)
dev.off()

inuse <- colnames(lsec)[lsec@meta.data$seurat_clusters != 6]
saveRDS(inuse, file.path(rdsDir, "selected-lsec.rds"))

mtpercent <- getGenePercent(counts, totals)
d <- cbind(cellanno, mt = mtpercent, total = totals)
q <- ggplot(data = d %>% filter(Cell.type == "LSEC")) +
  geom_point(aes(y = mt, x = total, colour = barcode %in% inuse), size = .5) +
  scale_x_log10() +
  scale_colour_discrete(name = "selected") +
  facet_wrap(~ sample, scales = "free_x")
ggsave(figpattern("selected-lsec-mito-percent.png"), q, width = 10, height = 9, dpi = 200)
