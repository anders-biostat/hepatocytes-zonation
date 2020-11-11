library(rmatio)
library(purrr)
library(dplyr)
library(tidyr)

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
    central = read.table("from_Ki/Gene_of_interest_Ki_zonated_pericentral")$V2,
    portal = read.table("from_Ki/Gene_of_interest_Ki_zonated_periportal")$V2
    )
  x <- lapply(x, tolower)
  lapply(x, as.character)
}

getZoneMarkers <- function() {
  list(itzkevitz = getItzkeitzMarkers(),
    lsec = getLSECMarkers())
}

## eta calculations
scaleByQuantile <- function(mat, prob, MARGIN = 2) {
  qs <- apply(mat, MARGIN, function(x) quantile(x, prob))
  if (MARGIN == 1)
    return(mat / qs)
  if (MARGIN == 2)
    return(t(t(mat) / qs))
}

calcEta <- function(x, markers, prob = 1) {
  x <- map(markers, ~ x[.x, ])
  x <- map(x, scaleByQuantile, MARGIN = 1, prob = prob)
  x <- map(x, colSums)
  res <- x$portal / (x$portal + x$central)
  res
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
  d$odds <- eta / (1 - eta)
  d$total <- totals
  d$cell <- colnames(fracs)
  d <- d %>% gather(gene, frac, -eta, -total, -cell, -odds)
  d$gene <- factor(d$gene, levels = unlist(allgenes))
  d
}


getGenePercent <- function(counts, totals, pattern = "mt-") {
  mtgenes <- grepl(pattern, rownames(counts))
  mtpercent <- colSums(counts[mtgenes, ]) / totals
  mtpercent
}

createSplineMat <- function(position, degree) {
  X <- splines::bs(position, df = degree,
    Boundary.knots = c(0,1), intercept = TRUE)
  colnames(X) <- paste0("X", colnames(X))
  X
}

assignPosition <- function(counts, cellanno, totals, fracs, markers) {
  conditions <- unique(cellanno$condition)
  etadf <- map(set_names(conditions),
    function(condition) {
      cells <- cellanno$condition == condition
      makeEtaDF(cells, cellanno, totals, fracs, markers)
    })
  etadf <- do.call(rbind, etadf)
  res <- etadf %>% inner_join(cellanno, by = c(cell = "barcode")) %>%
    group_by(sample) %>%
    mutate(etaq = cume_dist(eta)) %>%
    ungroup %>%
    select(eta, etaq, cell) %>%
    unique %>%
    inner_join(cellanno, by = c(cell = "barcode"))
  res$total <- totals
  res$mouse <- factor(res$Mouse.ID)
  res$condition <- factor(res$condition)
  res
}

capitalize <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
