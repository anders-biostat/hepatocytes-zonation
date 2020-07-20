
rdsDir <- file.path("results", "rds")
figdir <- file.path("results", "figs")
tabledir <- file.path("results", "tables")

loadConfig <- function()
  yaml::read_yaml("code/config.yaml")

figpattern <- function(pattern, ...) {
  x <- glue::glue(pattern, ..., .envir = parent.frame())
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

getBarcodes <- function(celltype) {
  readRDS(file.path(rdsDir,
    glue::glue("selected-{celltype}.rds", celltype = celltype)))
}

getCountsCelltype <- function(celltype) {
  counts <- getCounts()
  counts[, getBarcodes(celltype)]
}

getCellAnnoCelltype <- function(celltype) {
  cellanno <- getCellAnno()
  cellanno[match(getBarcodes(celltype), cellanno$barcode), ]
}

loads <- list(
  counts = getCounts,
  "counts-heps" = function() getCountsCelltype("heps"),
  "counts-lsec" = function() getCountsCelltype("lsec"),
  meta = function() readRDS(file.path(rdsDir, "meta.rds")),
  cellanno = getCellAnno,
  "cellanno-heps" = function() getCellAnnoCelltype("heps"),
  "cellanno-lsec" = function() getCellAnnoCelltype("lsec"),
  markers = function() getZoneMarkers()
)

loadFiles <- function(...) {
  l <- list(...)
  x <- lapply(l, function(x) loads[[x]]())
  names(x) <- l
  if (length(l) == 1) return(x[[1]])
  x
}
