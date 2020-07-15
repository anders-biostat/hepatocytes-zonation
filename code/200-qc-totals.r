library(ggplot2)
library(dplyr)
library(tidyr)

source("code/assets.r")
source("code/func.r")
figpath <- function(x) file.path("results", "figs", x)

counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

cellanno$batch <- factor(cellanno$Experimental.Batch)
cellanno$mouse <- factor(cellanno$Mouse.ID)

totals <- colSums(counts)

qraw <- qplot(y = totals, log = 'y', colour = cellanno$batch) +
  scale_colour_discrete(name = "batch")
ggsave(filename = figpath("total-counts.png"), qraw, dpi = 200, width = 9, height = 6)

qdens <- ggplot(
  data = cbind(
    cellanno, total = totals)) +
  stat_density(aes(x = log10(total),
    colour = mouse,
    linetype = batch), geom = "line") +
  facet_wrap(~condition)

ggsave(filename = figpath("total-counts-density.png"),
  qdens, dpi = 200, width = 9, height = 9)

m <- model.matrix(~ sample + 0,
  data = data.frame(
      sample = paste(cellanno$condition, cellanno$mouse, cellanno$batch)))

## Sum up pseudobulks and plot hclust for moderate expression
pbulks <- counts %*% m
colnames(pbulks) <- gsub("sample", "", colnames(pbulks))
normedBulks <- t(pbulks) /colSums(pbulks)
expressed <- rowMeans(pbulks) > 10

png(figpath("qc-dendrogram-bulks.png"), res = 200, width = 7, height = 7, units = "in")
par(mar = c(0,0,0,10))
plot(as.dendrogram(
  hclust(dist(normedBulks[, expressed]))), horiz = TRUE)
dev.off()
