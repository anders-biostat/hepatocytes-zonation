
library(rmatio)

getItzkevitzParams <- function() {
  zoneParamsFile <- "from_Ki/Itzkovitz zonation metlab/matlab_code//Zonation_params.mat"
  zonationParams <- read.mat(zoneParamsFile)
  zonationParams
}



getMarkers <- function() {
  zonationParams <- getItzkevitzParams()
  markers <- list(
    central = zonationParams$genes_cv,
    portal = zonationParams$genes_pn)
  markers <- lapply(markers, unlist)
  markers
}

## Here we reproduce Itykovitz 2017/2018 algo
etaDensity <- function(eta, zone, gammaParams) {
  shape <- gammaParams[zone, 1]
  scale <- gammaParams[zone, 2]
  dgamma(eta, shape = shape, scale = scale)
}

calcEta <- function(x) {
  res <- x$portal / (x$portal + x$central)
  normEta(res)
}

normEta <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

calcEtaLikelihoods <- function(eta) {
  zones <- 1:8
  do.call(cbind, map(zones, ~ etaDensity(eta, .x)))
}

getPriors <- function(n) {
  cellnums <- 2 * seq(1, n)
  cellnums / sum(cellnums)
}

calcPosteriors <- function(zoneLiks, priors) {
  x <- t(t(zoneLiks) * priors)
  x / rowSums(x)
}
