
rdsDir <- file.path("results", "rds")
figdir <- file.path("results", "figs")
tabledir <- file.path("results", "tables")

figpattern <- function(pattern, ...) {
  x <- glue::glue(pattern, ...)
  file.path(figdir, x)
}

tablepattern <- function(pattern, ...) {
    x <- glue::glue(pattern, ...)
    file.path(tabledir, x)
}

getCounts <- function() {
  countsFile <- file.path(rdsDir, "allcouts.rds")
  counts <- readRDS(countsFile)
  rownames(counts) <- tolower(rownames(counts))
  counts
}

getCellAnno <- function() {
  cellanno <- readRDS(file.path(rdsDir, "cellanno.rds"))
  cellanno$condition <- with(cellanno, paste(Cell.type, Genotype, sep = "|"))
  cellanno
}

loads <- list(
  counts = getCounts,
  meta = function() readRDS(file.path(rdsDir, "meta.rds")),
  cellanno = getCellAnno,
  markers = function() getZoneMarkers()
)

loadFiles <- function(...) {
  l <- list(...)
  x <- lapply(l, function(x) loads[[x]]())
  names(x) <- l
  if (length(l) == 1) return(x[[1]])
  x
}
