library(Matrix)
library(ggplot2)
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
  markers$itzkevitz <- lapply(markers$itzkevitz, function(a) intersect(a, rownames(counts)))
  cat("OK\n")
}

landmarks <- add6LandMarks(list())
cellanno$condition <- with(cellanno, paste(Cell.type, Genotype, sep = "|"))

## all are present!
intersect(unlist(landmarks), rownames(counts))

d <- makeEtaDF(
  cellanno$condition == "HEP|Wildtype",
  cellanno, totals, fracs, markers$itzkevitz,
  othergenes = list(unlist(landmarks)))

q <- d %>%
  filter(gene %in% unlist(landmarks)) %>%
  ggplot(data = .) +
  geom_point(aes(x = eta, y = frac), size = .2) +
  facet_wrap(~gene, scales = "free_y")

ggsave(figpath("landmarks-6-in-HEPs-WT-vs-eta.png"),
  q, width = 8, height= 6, dpi = 200)

d <- makeEtaDF(
  cellanno$condition == "HEP|DoubleKO",
  cellanno, totals, fracs, markers$itzkevitz,
  othergenes = list(unlist(landmarks)))

q <- d %>%
  filter(gene %in% unlist(landmarks)) %>%
  ggplot(data = .) +
  geom_point(aes(x = eta, y = frac), size = .2) +
  facet_wrap(~gene, scales = "free_y")

ggsave(figpath("landmarks-6-in-HEPs-DoubleKO-vs-eta.png"),
  q, width = 8, height= 6, dpi = 200)
