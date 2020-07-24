library(DESeq2)
library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation")

makeDE <- function(counts, cellanno) {
  m <- model.matrix(~ 0 + sample, data = cellanno)
  pbulks <- counts %*% m
  pbulks <- pbulks[rowMeans(pbulks) > .5, ]
  ## "sampleLSEC|Wildtype 3 9" -> "LSEC|Wildtype 3 9" 
  colnames(pbulks) <- gsub("sample", "", colnames(pbulks))
  bulkanno <- unique(cellanno[c("Genotype", "sample")])
  rownames(bulkanno) <- bulkanno$sample
  pbulks <- pbulks[, rownames(bulkanno)]

  dds <- DESeq2::DESeqDataSetFromMatrix(pbulks, colData = bulkanno, design = ~ Genotype)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("Genotype", "Wildtype", "DoubleKO"))

  res <- res %>% as.data.frame %>% arrange(pvalue)
  list(result = res, pbulks = pbulks, anno = bulkanno)
}

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
results <- lapply(dat, do.call, what = makeDE)

saveRDS(results, file.path(rdsDir, "pbulks-deseq.rds"))

