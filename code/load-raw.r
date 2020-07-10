## This file reads raw counts into a single table and creates meta table for cells.
## There were  seurat objects provided with clustering information,
## no source code for that, but it looks like some cells were filtered out,
## cell identities can not be mapped back to batches in obvious way, only guessing
## from gene expression, cell number and condition.
library(purrr)
library(Matrix)
meta <- read.table(
  "from_Ki/Ki experimental design.txt", sep = "\t",
  header = TRUE)
meta$ID <- sub(".+?_", "", meta$ID, perl=TRUE)

countDir <- "from_Ki/Raw data/"

## there two files with _2 suffix in filename. Ki's said, he made a mistake
## and one needs to use the _2 files, as in meta.
x <- list.files("from_Ki/Raw data/", recursive = TRUE) %>%
  strsplit(.Platform$file.sep)
x <- Filter(function(a) grepl("csv$", a[[2]]), x)
ids <- sub(".+?_", "", map_chr(x, 2), perl=TRUE) %>% sub(".csv", "", x = .)
x <- x[ids %in% meta$ID]

counts <- purrr::map(x, ~ read.csv(file.path(countDir, .x[1], .x[2])))

allGenes <- Reduce(union, map(counts, "X"))

shapeCounts <- function(x) {
  rownames(x) <- x$X
  x$X <- NULL
  res <- Matrix(0, nrow = length(allGenes), ncol = ncol(x), sparse = TRUE)
  rownames(res) <- allGenes
  colnames(res) <- colnames(x)
  res[rownames(x),] <- as.matrix(x)
  res
}

sparseCounts <- map(counts, shapeCounts)
allcounts <- do.call(cbind, sparseCounts)
saveRDS(allcounts, file.path("results", "rds", "allcouts.rds"))


meta$ID <- sub("_.+", "", meta$ID, perl=TRUE)
barcodes <- gsub("Batch", "", colnames(allcounts))
cellanno <- data.frame(
  barcode = barcodes,
  ID = gsub("_Cell.+", "", barcodes),
  stringsAsFactors = FALSE)
cellanno <- cellanno %>% left_join(meta)
saveRDS(cellanno, file.path("results", "rds", "cellanno.rds"))

saveRDS(meta, file.path("results", "rds", "meta.rds"))
