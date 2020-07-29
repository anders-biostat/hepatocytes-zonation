library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

library(foreach)
library(doParallel)
registerDoParallel(8)


figdir <- file.path("results", "figs", "functional")
dir.create(figdir, showWarnings = FALSE)

{cat("reading...")
  counts <- loadFiles("counts-heps")
  cellanno <- loadFiles("cellanno-heps")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)
  funcgenes <- loadFiles("functional-genes")
  cat("OK\n")
}

conditions <- c("HEP|DoubleKO", "HEP|Wildtype")
etaDF <- map(set_names(conditions),
  function(condition) {
    cells <- cellanno$condition == condition
    makeEtaDF(cells, cellanno, totals, fracs, markers, othergenes = unlist(funcgenes))
  })
etaDFwide <- map(etaDF, ~ pivot_wider(.x, names_from = "gene", values_from = "frac"))
cellannoByCond <- split(cellanno, cellanno$condition)


geneZonePlot <- function(x) {
  x <- x %>% mutate(etaq = cume_dist(eta))
  q <- ggplot(data = x) +
    geom_smooth(aes(x = etaq, y = sqrt(frac)), method = "loess") +
    geom_point(aes(x = etaq, y = sqrt(frac)), shape = ".") +
    facet_wrap(~gene, scales = "free_y")
  q
}

plotCompareExpressionOverEtaRank <- function(etaDF, alpha = .5) {
  x <- bind_rows(etaDF, .id = "condition")
  x <- x %>%
    mutate(etaq = cume_dist(eta))
  colours <- set_names(c("dodgerblue", "black"), conditions)
  q <- qplot(
    data = x,
    x = etaq,
    colour =  condition,
    y = sqrt(frac),
    alpha = I(alpha),
    size = I(.1)) +
    scale_colour_manual(values = colours) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    facet_wrap(~gene, scales = "free_y")
  q + geom_smooth(aes(x = etaq, y = sqrt(frac)), method = "loess")
}

underscorize <- function(x) gsub("[ ]", "_", x)
genesets <- c(funcgenes[[1]], funcgenes[[2]])
## all except wnt from Shalev
genesets <- genesets[map_int(genesets, length) < 100]

foreach(geneset = names(genesets)) %dopar% {
  x <- map(etaDF, ~ filter(.x, gene %in% genesets[[geneset]]))
  ## q <- geneZonePlot(x)
  q <- plotCompareExpressionOverEtaRank(x)
  N <- length(genesets[[geneset]])
  ggsave(filename = figpattern("{underscorize(geneset)}-frac-eta.png"),
         q,
         width = min(30, 4 + sqrt(N) * 2.5),
         height = min(30, 4 + (sqrt(N)) * 2),
         dpi = 200)
}

averageGenesOverZones <- function(d, zoneNum) {
  x <- d %>%
    mutate(etaq = cume_dist(eta), zone = cut(etaq, 0:zoneNum / zoneNum))
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

genesets <- c(funcgenes[[1]], funcgenes[[2]])
for (condition in conditions) {
  foreach(geneset = names(genesets)[map_int(genesets, length) > 1]) %dopar% {
    g <- genesets[[geneset]]
    x <- averageGenesOverZones(etaDF[[condition]] %>% filter(gene %in% g), zoneNum = 6)
    x <- x[rowSums(x) > 0,]
    fn <- figpattern("heatmap-zones-{underscorize(geneset)}-{condition}.png")
    pheatmap(
      filename = fn,
      x,
      scale = "row",
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      width = 10, height = 10, res = 200)
  }
}
