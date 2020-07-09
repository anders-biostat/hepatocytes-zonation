library(purrr)
library(Matrix)
meta <- read.table(
  "from_Ki/Ki experimental design.txt", sep = "\t",
  header = TRUE)
saveRDS(meta, file.path("results", "rds", "meta.rds"))

countDir <- "from_Ki/Raw data/"

x <- list.files("from_Ki/Raw data/", recursive = TRUE) %>%
  strsplit(.Platform$file.sep)
x <- Filter(function(a) grepl("csv$", a[[2]]), x)

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


barcodes <- gsub("Batch", "", colnames(counts))
cellanno <- data.frame(
  barcode = barcodes,
  ID = gsub("_Cell.+", "", barcodes),
  stringsAsFactors = FALSE)
cellanno <- cellanno %>% left_join(meta)
saveRDS(cellanno, file.path("results", "rds", "cellanno.rds"))
