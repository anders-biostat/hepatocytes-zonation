
library(rmatio)

getItzkevitzParams <- function() {
  zoneParamsFile <- "from_Ki/Itzkovitz zonation metlab/matlab_code//Zonation_params.mat"
  zonationParams <- read.mat(zoneParamsFile)
  zonationParams
}



getItzkeitzMarkers <- function() {
  zonationParams <- getItzkevitzParams()
  markers <- list(
    central = zonationParams$genes_cv,
    portal = zonationParams$genes_pn)
  markers <- lapply(markers, unlist)
  markers <- lapply(markers, tolower)
  markers
}

getLSECMarkers <- function() {
  x <- list(
    portal = read.table("from_Ki/Gene_of_interest_Ki_zonated_periportal")$V2,
    central = read.table("from_Ki/Gene_of_interest_Ki_zonated_pericentral")$V2
    )
  lapply(x, as.character)
}

getZoneMarkers <- function() {
  list(itzkevitz = getItzkeitzMarkers(),
    lsec = getLSECMarkers())
}

## Here we reproduce Itykovitz 2017/2018 algo
etaDensity <- function(eta, zone, gammaParams) {
  shape <- gammaParams[zone, 1]
  scale <- gammaParams[zone, 2]
  dgamma(eta, shape = shape, scale = scale)
}

calcEta <- function(x, markers) {
  x <- map(markers, ~ x[.x, ])
  x <- map(x, normByMax)
  x <- map(x, colSums)
  res <- x$portal / (x$portal + x$central)
}

calcEtaAlt <- function(x, totals) {
  res <- (x$portal + .5) / (0.5 + x$central)
  res
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

normByMax <- function(x)
  x / apply(x, 1, max)

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

sparse2df <- function(x) as.data.frame(as.matrix(x))

makeEtaDF <- function(idx, cellanno, counts, fracs, markers, etamethod = calcEta) {
  cellanno <- cellanno[idx, ]
  counts <- counts[, idx]
  fracs <- fracs[, idx]
  totals <- totals[idx]

  ## just ratio of portal/(central + portal)
  eta <- etamethod(fracs, markers)
  d <- do.call(rbind, map(markers, ~ fracs[.x,]))
  d <- sparse2df(t(d))
  d$eta <- eta
  d$total <- totals
  d$cell <- colnames(counts)
  d <- d %>% gather(gene, frac, -eta, -total, -cell)
  d$gene <- factor(d$gene, levels = unlist(markers))
  d
}
