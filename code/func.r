
library(rmatio)

getItzkevitzParams <- function() {
  zoneParamsFile <- "from_Ki/Itzkovitz zonation metlab/matlab_code//Zonation_params.mat"
  zonationParams <- read.mat(zoneParamsFile)
  zonationParams
}


add6LandMarks <- function(x) {
  x$central <- c(x$central, "glul", "cyp2e1")
  x$portal <- c(x$portal, "ass1", "asl", "alb", "cyp2f2")
  lapply(x, unique)
}

getItzkeitzMarkers <- function() {
  zonationParams <- getItzkevitzParams()
  markers <- list(
    central = zonationParams$genes_cv,
    portal = zonationParams$genes_pn)
  markers <- lapply(markers, unlist)
  markers <- lapply(markers, tolower)
  markers <- lapply(markers, transformGeneAltNames)
  markers
}

## substitute values given by Ki, return vector, input vector
transformGeneAltNames <- function(markers) {
  ## [1] "cml2" "cyb5" "c1s"  "dak"
  ## use only genes we have
  ## Ki found alternative names
  i <- tolower(c("Nat8f2", "Cyb5a", "C1s1", "Tkfc"))
  names(i) <- c("cml2","cyb5","c1s", "dak")
  subGene <- function(g)
    ifelse(g %in% names(i), i[g], g)
  purrr::map_chr(markers, subGene)
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

makeEtaDF <- function(idx, cellanno, totals, fracs, markers, othergenes=NULL,
                      etamethod = calcEta) {
  cellanno <- cellanno[idx, ]
  fracs <- fracs[, idx]
  totals <- totals[idx]
  ## just ratio of portal/(central + portal)
  eta <- etamethod(fracs, markers)
  allgenes <- unique(unlist(c(markers, othergenes)))
  d <- fracs[allgenes, ]
  d <- sparse2df(t(d))
  d$eta <- eta
  d$total <- totals
  d$cell <- colnames(fracs)
  d <- d %>% gather(gene, frac, -eta, -total, -cell)
  d$gene <- factor(d$gene, levels = unlist(allgenes))
  d
}
