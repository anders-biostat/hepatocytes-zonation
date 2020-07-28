## cluster 4 expresses markers of stellate cells, Kupfer cells and LSEC
## we remove only it and save cell names for future subsettings
source("code/func.r")
source("code/assets.r")
library(Seurat)
library(ggplot2)
figdir <- file.path("results", "figs", "qc")

counts <- loadFiles("counts")
totals <- colSums(counts)
cellanno <- loadFiles("cellanno")

s <- readRDS(file.path(rdsDir, "seurats-qc.rds"))
heps <- s$HEP
inuse <- colnames(heps)[heps@meta.data$seurat_clusters != 4]

mtpercent <- getGenePercent(counts, totals)
muppercent <- getGenePercent(counts, totals, "^mup")
d <- cbind(cellanno, mt = mtpercent, total = totals, mup = muppercent)
q <- ggplot(data = d %>% filter(Cell.type == "HEP")) +
  geom_point(aes(y = mt, x = total, colour = barcode %in% inuse), size = .5) +
  scale_x_log10() +
  scale_colour_discrete(name = "selected") +
 facet_wrap(~ sample, scales = "free_x")
ggsave(figpattern("selected-heps-mito-percent.png"), q, width = 10, height = 9, dpi = 200)


q <- ggplot(data = d %>% filter(Cell.type == "HEP")) +
  geom_point(aes(y = mup, x = total, colour = barcode %in% inuse), size = .5) +
  scale_x_log10() +
  scale_colour_discrete(name = "selected") +
  facet_wrap(~ sample, scales = "free_x")
ggsave(figpattern("selected-heps-mup-percent.png"), q, width = 10, height = 9, dpi = 200)

saveRDS(inuse, file.path(rdsDir, "selected-heps.rds"))
