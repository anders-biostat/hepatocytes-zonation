
## Itzkovitz publication related
etaDensity <- function(eta, zone, gammaParams) {
  shape <- gammaParams[zone, 1]
  scale <- gammaParams[zone, 2]
  dgamma(eta, shape = shape, scale = scale)
}

normEta <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

getPriors <- function(n) {
  cellnums <- 2 * seq(1, n)
  cellnums / sum(cellnums)
}

calcPosteriors <- function(zoneLiks, priors) {
  x <- t(t(zoneLiks) * priors)
  x / rowSums(x)
}

normByMax <- function(x)
  x / apply(x, 1, max)

calcEtaLikelihoods <- function(eta) {
  zones <- 1:8
  do.call(cbind, map(zones, ~ etaDensity(eta, .x)))
}

## 2018 publ. for reference
## Itzkovitz use fracs in matlab file, not doc'ed in publ.
shalevZoneAssign <- function(mexpr, markers, fracs) {
  ## 1) normalize by gene frac max across heps
  nmexpr <- map(mexpr, normByMax)
  ## 2) compute eta score
  eta <- calcEta(map(nmexpr, colSums))
  eta <- calcEtaAlt(map(map(markers, ~ counts[.x, ]), colSums))
  ## 2) posterior probability of zones P cells x zones
  zoneLiks <- calcEtaLikelihoods(eta)
  priors <- getPriors(8)
  posteriors <- calcPosteriors(zoneLiks, priors)
  ## 3) column-normalize P  -> W
  weights <- posteriors / colSums(posteriors, na.rm=TRUE)
  weights[is.na(weights)] <- 0
  ## 4) multiply expression x W
  exprByZone <-  fracs %*% weights
}
