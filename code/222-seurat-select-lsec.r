## cluster 6 expresses markers of stellate cells (dcn), so we remove them.
## There are other smaller clusters, but they express stab2, kdr, cdh5, so
## I do not exclude them.
source("code/func.r")
source("code/assets.r")
library(Seurat)

counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

s <- readRDS(file.path(rdsDir, "seurats-qc.rds"))
lsec <- s$LSEC

table(lsec@meta.data$seurat_clusters)
##
## 0   1   2   3   4   5   6   7   8   9  10  11  12
## 580 397 319 309 301 254 122  55  55  20  19  18  18

inuse <- colnames(lsec)[lsec@meta.data$seurat_clusters != 6]
saveRDS(inuse, file.path(rdsDir, "selected-lsec.rds"))
