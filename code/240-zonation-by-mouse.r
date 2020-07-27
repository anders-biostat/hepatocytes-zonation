library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation")

read <- function(cellType) {
  cat("reading...", cellType)
  counts <- loadFiles(glue::glue("counts-{cellType}"))
  cellanno <- loadFiles(glue::glue("cellanno-{cellType}"))
  cellanno$sample <- paste(cellanno$condition, cellanno$Mouse.ID, cellanno$Experimental.Batch)
  cat("OK\n")
  list(cellanno = cellanno, counts = counts)
}

celltypes <- c("heps", "lsec") %>% set_names
dat <- lapply(celltypes, read)
