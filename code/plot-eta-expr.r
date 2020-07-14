library(Matrix)
library(purrr)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")
figdir <- file.path("results", "figs")
figpath <- function(x) file.path(figdir, x)


{cat("reading...")
  counts <- loadFiles("counts")
  meta <- loadFiles("meta")
  cellanno <- loadFiles("cellanno")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")
  cat("OK\n")
}

## genes we do not have:
setdiff(unlist(markers$itzkevitz), rownames(counts))
## [1] "cml2" "cyb5" "c1s"  "dak"
## use only genes we have
markers$itzkevitz <- lapply(markers$itzkevitz, function(a) intersect(a, rownames(counts)))


## select condition
sets <- list(
  "HEP-WT"  = cellanno$Cell.type == "HEP"  & cellanno$Genotype != "DoubleKO",
  "HEP-DKO" = cellanno$Cell.type == "HEP"  & cellanno$Genotype == "DoubleKO",
  "LSEC-WT" = cellanno$Cell.type == "LSEC" & cellanno$Genotype != "DoubleKO",
  "LSEC-DKO" = cellanno$Cell.type == "LSEC" & cellanno$Genotype == "DoubleKO"
)

## plot marker expression vs eta
for(i in names(sets)[1:2]) {
  d <- makeEtaDF(sets[[i]], cellanno, totals, fracs, markers$itzkevitz)
  ## plot expression along eta and smoothing
  q <- d %>%
    ggplot(data= .) +
    geom_smooth(aes(x = rank(eta), y = frac), method = "loess") +
    facet_wrap(~gene, scales = "free_y")
  ggsave(filename = figpath(sprintf("marker-frac-smooth-eta-%s.png", i)),
    q, width = 12, height = 8, dpi = 200)

  q <- d %>%
    ggplot(data= .) +
    geom_smooth(aes(x = rank(eta), y = frac), method = "loess") +
    geom_point(aes(x = rank(eta), y = frac), shape = ".") +
    facet_wrap(~gene, scales = "free_y")
  ggsave(filename = figpath(sprintf("marker-frac-smooth-points-eta-%s.png", i)),
    q, width = 12, height = 8, dpi = 200)

  q <- d %>%
    ggplot(data= .) +
    geom_point(aes(x = rank(eta), y = frac*total), shape = ".") +
    facet_wrap(~gene, scales = "free")

  ggsave(filename = figpath(sprintf("marker-count-vs-eta-%s.png", i)),
    q, width = 12, height = 8, dpi = 200)

}


averageGenesOverZones <- function(d, zoneNum) {
  x <- d %>%
    mutate(etaq = rank(eta)/length(eta), zone = cut(etaq, 0:zoneNum/zoneNum))
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
for(i in names(sets)[1:2]) {
  d <- makeEtaDF(sets[[i]], cellanno, totals, fracs, markers$itzkevitz)
  m <- averageGenesOverZones(d, zoneNum = 6)

  pheatmap(
    m,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE)

  fn <- figpath(sprintf("marker-heatmap-zones-%s.png", i))

  geneType <- data.frame(
    marker = c("portal", "central")[1 + rownames(m) %in% markers$itzkevitz$central])
  rownames(geneType) <- rownames(m)

  pheatmap(
    filename = fn,
    annotation_row = geneType,
    m,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    width = 10, height = 8, res = 200)

  x <- d %>% pivot_wider(names_from = gene, values_from = frac)
  x <- x[,-(1:3)]
  x <- as.matrix(x)
  x <- cor(x)

  pheatmap(x,
    filename = figpath(sprintf("marker-corr-%s.png", i)),
    breaks = seq(-.2,.2, length.out = 100),
    cluster_cols = FALSE,
    cluster_rows = FALSE)

}


plotCompareExpressionOverEta <- function() {

  ds <- lapply(names(sets)[1:2], function(i)
    d <- makeEtaDF(sets[[i]], cellanno, totals, fracs, markers$itzkevitz))
  names(ds) <- names(sets)[1:2]

  x <- bind_rows(ds, .id = "condition")
  q <- qplot(
    data = x,
    x = eta,
    colour =  condition,
    y = frac,
    alpha = .5,
    size = I(.2)) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    facet_wrap(~gene, scales = "free_y")
  q
}

q <- plotCompareExpressionOverEta()
ggsave(filename = figpath("compare-hep-marker-vs-eta.png"),
  plot = q, width = 12, height = 10, dpi = 200)
