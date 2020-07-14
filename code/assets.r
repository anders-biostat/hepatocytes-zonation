
rdsDir <- file.path("results", "rds")


getCounts <- function() {
  countsFile <- file.path(rdsDir, "allcouts.rds")
  counts <- readRDS(countsFile)
  rownames(counts) <- tolower(rownames(counts))
  counts
}

loads <- list(
  counts = getCounts,
  meta = function() readRDS(file.path(rdsDir, "meta.rds")),
  cellanno = function() readRDS(file.path(rdsDir, "cellanno.rds")),
  markers = function() getZoneMarkers()
)

loadFiles <- function(...) {
  l <- list(...)
  x <- lapply(l, function(x) loads[[x]]())
  names(x) <- l
  if (length(l) == 1) return(x[[1]])
  x
}
