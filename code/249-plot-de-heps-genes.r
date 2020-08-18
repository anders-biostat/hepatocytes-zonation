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


figdir <- file.path("results", "figs", "de")
dir.create(figdir, showWarnings = FALSE)

{cat("reading...")
  counts <- loadFiles("counts-heps")
  cellanno <- loadFiles("cellanno-heps")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)
  cat("OK\n")
}

de <- readRDS("results/rds/pbulks-deseq.rds")
degenes <- rownames(de$heps$result[1:100,])

conditions <- c("HEP|DoubleKO", "HEP|Wildtype")
etaDF <- map(set_names(conditions),
  function(condition) {
    cells <- cellanno$condition == condition
    makeEtaDF(cells, cellanno, totals, fracs, markers, othergenes = degenes)
  })
etaDFwide <- map(etaDF, ~ pivot_wider(.x, names_from = "gene", values_from = "frac"))
cellannoByCond <- split(cellanno, cellanno$condition)

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

x <- map(etaDF, ~ filter(.x, gene %in% degenes))
q <- plotCompareExpressionOverEtaRank(x)
ggsave(filename = figpattern("hep-de-frac-eta.png"),
  q,
  width = 20,
  height = 20,
  dpi = 100)
