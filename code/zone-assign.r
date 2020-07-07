library(Matrix)
library(purrr)
library(dplyr)
library(tidyr)

{cat("reading...")
  rdsDir <- file.path("results", "rds")
  countsFile <- file.path(rdsDir, "allcouts.rds")
  counts <- readRDS(countsFile)
  meta <- readRDS(file.path(rdsDir, "meta.rds"))
  cellanno <- readRDS(file.path(rdsDir, "cellanno.rds"))
  cat("OK\n")
}

## we work only with hepatocytes
heps <- Negate(is.na)(
  cellanno$Cell.type == "HEP"
  & cellanno$Genotype == "Wildtype")
cellanno <- cellanno[heps,]
counts <- counts[, heps]

## 1) normalize by cell totals
## 2) posterior probability of zones P cells x zones
## 3) column-normalize P (????) -> W
## 4) multiply expression x W

## get parameters from Itzkovitz data
rownames(counts) <- tolower(rownames(counts))
totals <- colSums(counts)
fracs <- t(t(counts) / totals)

source("code/func.r")
gammaParams <- getItzkevitzParams()$Gamma_params
markers <- getMarkers()
## use only genes we have
markers <- lapply(markers, function(a) intersect(a, tolower(rownames(counts))))

etaDensity <- function(eta, zone) {
  shape <- gammaParams[zone, 1]
  scale <- gammaParams[zone, 2]
  dgamma(eta, shape = shape, scale = scale)
}

calcEta <- function(x) {
  mexpr <- map(markers, ~ colSums(x[.x,]))
  mexpr$portal / (mexpr$portal + mexpr$central)
}

normEta <- function(x)
  (x - min(x)) / (max(x) - min(x))

calcEtaLikelihoods <- function(eta) {
  zones <- 1:8
  do.call(cbind, map(zones, ~ etaDensity(eta, .x)))
}

getPdfs <- function(eta, numZones) {
  pdfMat <- matrix(rep(0, numZones * length(eta)), ncol = numZones, nrow = length(eta))
  for (i in 1:numZones){
    pdf <- dgamma(eta, shape = gammaParams[i,1], scale = gammaParams[i,2])
    pdf[pdf == 0] <- min(pdf[pdf > 0])
    pdfMat[,i] <- pdf
  }
  pdfMat
}
