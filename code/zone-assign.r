library(Matrix)
library(purrr)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
figdir <- file.path("results", "figs")
figpath <- function(x) file.path(figdir, x)


{cat("reading...")
  rdsDir <- file.path("results", "rds")
  countsFile <- file.path(rdsDir, "allcouts.rds")
  counts <- readRDS(countsFile)
  meta <- readRDS(file.path(rdsDir, "meta.rds"))
  cellanno <- readRDS(file.path(rdsDir, "cellanno.rds"))
  rownames(counts) <- tolower(rownames(counts))
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)

  markers <- getZoneMarkers()
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

for(i in names(sets)[1:2]) {
  d <- makeEtaDF(sets[[i]], cellanno, counts, fracs, markers$itzkevitz)
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
  x <- x %>% group_by(gene, zone) %>%
    summarise(frac = mean(frac)) %>%
    ungroup()
  x <- x %>% pivot_wider(names_from = gene, values_from = frac)

  m <- as.matrix(x[, -1])
  rownames(m) <- x$zone
  t(m)
}


for(i in names(sets)[1:2]) {
d <- makeEtaDF(sets[[i]], cellanno, counts, fracs, markers$itzkevitz)
m <- averageGenesOverZones(d, zoneNum = 6)

pheatmap(
  m,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE)

fn <- figpath(sprintf("marker-heatmap-zones-%s.png", i))

geneType <- data.frame(
  marker = c("central", "portal")[1 + rownames(m) %in% markers$itzkevitz$central])
rownames(geneType) <- rownames(m)

pheatmap(
  filename = fn,
  annotation_row = geneType,
  m,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  width = 10, height = 8, res = 200)

}
