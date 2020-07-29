library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation")
markers <- loadFiles("markers")
markers$itzkevitz <- add6LandMarks(markers$itzkevitz)
names(markers)[1] <- "heps"

read <- function(cellType) {
  cat("reading...", cellType)
  counts <- loadFiles(glue::glue("counts-{cellType}"))
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  cellanno <- loadFiles(glue::glue("cellanno-{cellType}"))
  cellanno$sample <- paste(cellanno$condition, cellanno$Mouse.ID, cellanno$Experimental.Batch)
  cat("OK\n")
  list(cellanno = cellanno, counts = counts, fracs = fracs, totals = totals)
}

celltypes <- c("heps", "lsec") %>% set_names
dat <- lapply(celltypes, read)

z <- dat[[cellType]]
cellType <- "heps"
idx <- split(z$cellanno$barcode, z$cellanno$sample)
x <- map(idx, ~ makeEtaDF(.x, z$cellanno, z$totals, z$fracs, markers[[cellType]]))
x <- bind_rows(x, .id = "sample")
x <- inner_join(x, z$cellanno, by = c("cell" = "barcode", "sample" = "sample"))
x <- x %>% group_by(sample) %>%
  mutate(etaq = cume_dist(eta)) %>%
  ungroup

q <- ggplot(data = x %>% filter(gene == "pck1")) +
  geom_point(
    aes(x = etaq, y = frac, colour = condition, group = sample)
  ) +
  facet_wrap(~gene+sample, scales = "free_x")
q

q <- ggplot(data = x %>% filter(gene == "gc")) +
  geom_point(
    aes(x = etaq, y = frac, colour = condition, group = sample)
  ) +
  facet_wrap(~gene+sample, scales = "free_x")
q
