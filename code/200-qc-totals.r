## for the whole data set (HEP + LSEC) plot pseudobulks clusters and correlations

library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)

source("code/assets.r")
source("code/func.r")
figpath <- function(x) file.path("results", "figs", "qc", x)

counts <- loadFiles("counts")
cellanno <- loadFiles("cellanno")

cellanno$batch <- factor(cellanno$Experimental.Batch)
cellanno$mouse <- factor(cellanno$Mouse.ID)

totals <- colSums(counts)

qbox <- cbind(cellanno, total = totals) %>%
  ggplot(data = .) +
  geom_boxplot(aes(y = total, x = paste(condition, mouse, batch))) +
  scale_y_log10() +
  coord_flip()
ggsave(filename = figpath("total-counts-box.png"), qbox, dpi = 200, width = 9, height = 6)

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

## Sum up pseudobulks/normalize by total and plot hclust for moderate expression
pbulks <- counts %*% m
colnames(pbulks) <- gsub("sample", "", colnames(pbulks))
normedBulks <- t(pbulks) / colSums(pbulks)
expressed <- rowMeans(pbulks) > 10

x <- cor(as.matrix(t(normedBulks[, expressed])))
pheatmap::pheatmap(x, filename = file.path("results", "figs", "qc-correlation-all.png"),
  width = 9, height = 9)

png(figpath("qc-dendrogram-bulks.png"), res = 200, width = 7, height = 7, units = "in")
par(mar = c(0,0,0,10))
plot(as.dendrogram(
  hclust(dist(normedBulks[, expressed]))), horiz = TRUE)
dev.off()

cols <- gsub(".+\\|", "", rownames(normedBulks))
cols <- gsub(" .+", "", cols)
pca <- prcomp(normedBulks[, expressed])
qpca <- qplot(x = pca$x[,1],
  y = pca$x[,2],
  label = rownames(normedBulks),
  geom = "text",
  colour = cols,
  size = I(1.2)) +
  coord_cartesian(xlim = c(-.2,.4)) +
  xlab("pc1") + ylab("pc2")
ggsave(filename = figpath("pca-pbulks.png"),
  qpca, dpi = 200, width = 5, height = 5)
