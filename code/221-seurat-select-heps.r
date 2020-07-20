## cluster 4 expresses markers of stellate cells, Kupfer cells and LSEC
## we remove only it and save cell names for future subsettings
source("code/func.r")
source("code/assets.r")
library(Seurat)

counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

s <- readRDS(file.path(rdsDir, "seurats-qc.rds"))
heps <- s$HEP

inuse <- colnames(heps)[heps@meta.data$seurat_clusters != 4]
saveRDS(inuse, file.path(rdsDir, "selected-heps.rds"))
