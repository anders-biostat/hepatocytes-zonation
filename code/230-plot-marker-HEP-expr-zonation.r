## Zonation of hepatocytes
## - exploratory plots
## - zone assignment

library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation")

{cat("reading...")
  counts <- loadFiles("counts-heps")
  cellanno <- loadFiles("cellanno-heps")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)

  setdiff(unlist(markers), rownames(counts))
  ## character(0)

  # table with [type: portal/central], [gene name]
  markerAnnotation <- bind_rows(
    map(markers,  ~ data.frame(gene = .x, row.names = .x)),
    .id = "type")

  cat("OK\n")
}

## create data frames with expression fracs and eta
conditions <- c("HEP|DoubleKO", "HEP|Wildtype")
etaDF <- map(set_names(conditions),
  function(condition) {
    cells <- cellanno$condition == condition
    makeEtaDF(cells, cellanno, totals, fracs, markers)
  })
etaDFwide <- map(etaDF, ~ pivot_wider(.x, names_from = "gene", values_from = "frac"))
cellannoByCond <- split(cellanno, cellanno$condition)

## plot correlations for marker genes
for (condition in conditions) {
  m <- as.matrix(etaDFwide[[condition]][, -(1:3)])
  corr <- cor(m)
  pheatmap(corr,
    filename = figpattern("correlation-{cond}.png", cond = condition),
    annotation_row = markerAnnotation["type"],
    width = 8, height = 8)
  pheatmap(corr,
    filename = figpattern("correlation-{cond}-no-clust.png", cond = condition),
    annotation_row = markerAnnotation["type"],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    width = 8, height = 8)
}

## plot heatmaps with cells ordered by eta
## we scale sqrt of fracs by quantile
## change breaks to deal with outliers
for (condition in conditions) {
  cells <- order(etaDFwide[[condition]]$eta)
  m <- as.matrix(etaDFwide[[condition]][, -(1:3)])
  m <- t(scaleByQuantile(sqrt(m), .99, 2))[, cells]
  colnames(m) <- etaDFwide[[condition]]$cell
  zones <- data.frame(eta = etaDFwide[[condition]]$eta)[cells,, drop = FALSE]
  zoneNum <- 5
  zones <- zones %>%
    mutate(etaq = rank(eta) / length(eta), zone = cut(etaq, 0:zoneNum / zoneNum))
  rownames(zones) <- colnames(m)
  pheatmap(m,
    filename = figpattern("marker-heatmap-cells-{condition}.png"),
    annotation_row = markerAnnotation["type"],
    width = 12, height = 8,
    show_colnames = FALSE,
    breaks = seq(0, quantile(m, .99), length.out = 100),
    cluster_cols = FALSE)
  pheatmap(m,
    filename = figpattern("marker-heatmap-cells-{condition}-no-clust.png"),
    annotation_row = markerAnnotation["type"],
    show_colnames = FALSE,
    width = 12, height = 8,
    breaks = seq(0, quantile(m, .99), length.out = 100),
    cluster_cols = FALSE, cluster_rows = FALSE)
  zoneColor <- colorRampPalette(c("blue", "red"))(zoneNum)
  names(zoneColor) <- levels(zones$zone)
  pheatmap(m,
    filename = figpattern("marker-heatmap-cells-{condition}-clust.png"),
    annotation_row = markerAnnotation["type"],
    annotation_col = zones["zone"],
    annotation_colors = list(zone = zoneColor),
    show_colnames = FALSE,
    width = 12, height = 8,
    breaks = seq(0, quantile(m, .99), length.out = 100))
}

## pca vs eta
for (condition in conditions) {
  m <- as.matrix(etaDFwide[[condition]][, -(1:3)])
  scaledM <- scaleByQuantile(sqrt(m), .99)
  pr <- prcomp(scaledM)
  png(figpattern("pca-fracs-vs-eta-{cond}.png", cond = condition),
    width = 12, height = 12, res = 200, units = "in")
  par(mfrow = c(4, 4))
  batches <- filter(cellanno, condition == !!condition)$Experimental.Batch
  for (i in 1:16)
    plot(etaDFwide[[condition]]$eta, pr$x[, i],
      xlab = "eta",
      ylab = sprintf("PC%d", i),
      col = batches)
  dev.off()
}

## plot marker expression vs eta
for(i in conditions) {
  d <- etaDF[[i]]
  ## plot expression along eta and smoothing
  q <- d %>%
    ggplot(data = .) +
    geom_smooth(aes(x = rank(eta), y = sqrt(frac)), method = "loess") +
    facet_wrap(~gene, scales = "free_y")
  ggsave(filename = figpattern("marker-frac-smooth-eta-{i}.png", i = i),
    q, width = 12, height = 8, dpi = 200)

  q <- d %>%
    ggplot(data = .) +
    geom_smooth(aes(x = rank(eta), y = sqrt(frac)), method = "loess") +
    geom_point(aes(x = rank(eta), y = sqrt(frac)), shape = ".") +
    facet_wrap(~gene, scales = "free_y")
  ggsave(filename = figpattern("marker-frac-smooth-points-eta-{i}.png", i = i),
    q, width = 12, height = 8, dpi = 200)

  q <- d %>%
    ggplot(data = .) +
    geom_point(aes(x = rank(eta), y = frac * total), shape = ".") +
    ylab("count") +
    facet_wrap(~gene, scales = "free")

  ggsave(filename = figpattern(sprintf("marker-count-vs-eta-{i}.png", i = i)),
    q, width = 12, height = 8, dpi = 200)

}

colours <- set_names(c("dodgerblue", "black"), conditions)

plotCompareExpressionOverEta <- function(etaDF) {
  x <- bind_rows(etaDF, .id = "condition")
  q <- qplot(
    data = x,
    x = eta,
    colour =  condition,
    y = sqrt(frac),
    alpha = I(.5),
    size = I(.1)) +
    scale_colour_manual(values = colours) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    facet_wrap(~gene, scales = "free_y")
  q
}

q <- plotCompareExpressionOverEta(etaDF)
ggsave(filename = figpattern("compare-hep-marker-vs-eta.png"),
  plot = q, width = 16, height = 12, dpi = 200)


plotCompareExpressionOverEtaRank <- function(etaDF, alpha = .5) {
  x <- bind_rows(etaDF, .id = "condition")
  x <- x %>% group_by(condition) %>%
    mutate(etaRank = percent_rank(eta))

  q <- qplot(
    data = x,
    x = etaRank,
    colour =  condition,
    y = sqrt(frac),
    alpha = I(alpha),
    size = I(.1)) +
    scale_colour_manual(values = colours) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    facet_wrap(~gene, scales = "free_y")
  q
}

q <- plotCompareExpressionOverEtaRank(etaDF)
ggsave(filename = figpattern("compare-hep-marker-vs-eta-rank.png"),
  plot = q, width = 16, height = 12, dpi = 200)

qsmooth <- plotCompareExpressionOverEtaRank(etaDF, alpha = .1) +
  geom_smooth(aes(x = etaRank, y = sqrt(frac)), method = "loess")
ggsave(filename = figpattern("compare-hep-marker-vs-eta-rank-smooth.png"),
  plot = qsmooth, width = 16, height = 12, dpi = 200)

q <- q + facet_wrap(~gene, scales = "fixed")
ggsave(filename = figpattern("compare-hep-marker-vs-eta-rank-fixed.png"),
  plot = q, width = 16, height = 12, dpi = 200)

averageGenesOverZones <- function(d, zoneNum) {
  x <- d %>%
    mutate(etaq = rank(eta) / length(eta), zone = cut(etaq, 0:zoneNum / zoneNum))
  x <- x %>%
    select(gene, zone, frac) %>%
    group_by(gene, zone) %>%
    summarise(frac = mean(frac)) %>%
    ungroup()
  x <- x %>% pivot_wider(names_from = gene, values_from = frac)

  m <- as.matrix(x[, -1])
  rownames(m) <- x$zone
  t(m)
}

## plot heatmaps for average expression at quantiles of eta
for (condition in conditions) {
  m <- averageGenesOverZones(etaDF[[condition]], zoneNum = 6)
  fn <- figpattern("marker-heatmap-zones-{condition}.png")
  pheatmap(
    filename = fn,
    annotation_row = markerAnnotation[, "type", drop = FALSE],
    m,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    width = 10, height = 8, res = 200)
}
