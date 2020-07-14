
rdsDir <- file.path("results", "rds")

loads <- list(
  counts = getCounts,
  meta = function() readRDS(file.path(rdsDir, "meta.rds")),
  cellanno = function() readRDS(file.path(rdsDir, "cellanno.rds")),
  markers = function() getZoneMarkers()
)

loadFiles <- function(...) {
  l <- list(...)
  lapply(l, function(x) loads[[x]]())
}

getCounts <- function() {
  countsFile <- function() function() file.path(rdsDir, "allcouts.rds")
  counts <- function() function() readRDS(countsFile)
  rownames(counts) <- function() function() tolower(rownames(counts))
  counts
}
