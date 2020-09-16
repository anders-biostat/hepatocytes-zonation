library(mgcv)
library(locfit)
library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation-tests")
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
conditions <- unique(cellanno$condition)

samples <- unique(cellanno$sample)

normedMarkers <- map(c(1, .99), function(quant) {
  x <- map(set_names(samples), function(i) {
    scaleByQuantile(fracs[unlist(markers),cellanno$sample == i],
      quant, 2)
  })
  x <- do.call(cbind, x)
  x[x > 1] <- 1
  x <- t(x)
  x <- sparse2df(x)
  x
})


sums <- map(normedMarkers, function(x) {
  sums <- imap(markers, ~ data.frame(type = .y, value = rowSums(x[,.x])))
  sums <- map(sums, function(x) {
    x$cell <- rownames(x)
    x
  })
  sums <- do.call(rbind, sums)
  sums
})
sums[[1]]$quant <- 1
sums[[2]]$quant <- .99
sums <- bind_rows(sums) 
sums$quant <- factor(sums$quant)

q <- ggplot() +
stat_ecdf(data = sums, aes(x = value, colour = type, linetype = quant))

ggsave(filename = file.path("results", "figs", "qc", "scales-1-vs-.99.png"),
  q, width = 5, height = 5)

odds <- sums %>%
  spread(type, value) %>%
  mutate(odd = portal / central) %>%
  select(quant, odd, cell) %>%
  spread(quant, odd)

q <- qplot(
  x = rank(odds[,2]),
  y = rank(odds[,3]),
  xlab = ".99 quant",
  ylab = "1.0 quant"
)

ggsave(filename = file.path("results", "figs", "qc", "scales-1-vs-.99-odds.png"),
  q, width = 5, height = 5)
