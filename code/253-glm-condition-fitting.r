library(Matrix)
library(MASS)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation-tests")
dir.create(figdir, showWarnings = FALSE)
SPLINE_DEGREE <- 4


highexprp <- function(counts, cellanno, minproportion = 0.05) {
  propexpressed <- apply(counts > 0, 1, mean)
  highgenes <- propexpressed > minproportion
  highgenes
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


fitModels <- function(counts, cellanno, markers,
                      modelform, cores = 6, minproportion = 0.05) {
  highgenes <- highexprp(counts, cellanno, minproportion)
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  cellpositions <- assignPosition(counts, cellanno, totals, fracs, markers)
  X <- createSplineMat(cellpositions$etaq, SPLINE_DEGREE)
  splineVars <- colnames(X)
  mat <- cbind(cellpositions, X)
  form <- modelform(splineVars)
  models <- parallel::mclapply(
    which(highgenes),
    ## "gulo",
    function(.x) {
    try({
      m <- glm.nb(form, data = cbind(gene = counts[.x,], mat))
      list(coef = coef(m), vcov = vcov(m))
    })}, mc.cores = cores)
  models <- models[map(models, class) != "try-error"]
  list(models = models, mat = mat, splineVars = splineVars)
}


## HEPATOCYTES

{cat("reading...")
  d <- readCellType("heps")
  counts <- d$counts; cellanno <- d$cellanno
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)
  cat("OK\n")
}

## model with offset by total UMI and spline coefs for every mouse
mouseform <- function(splineVars) {as.formula(
  sprintf("gene ~ -1 + offset(log(total)) + (%s):condition",
    paste0(splineVars, collapse = "+")))
}

mousemodels <- fitModels(counts, cellanno, markers, modelform = mouseform, cores = 7,
  minproportion = .1)

betas <- map(mousemodels$models, "coef")
vars <- map(mousemodels$models, "vcov")

saveRDS(
  list(betas = betas, vars = vars, mat = mousemodels$mat),
  file.path(rdsDir, "hepa-glm-condition-nb-hvg.rds"))


## LSEC
{cat("reading...")
  d <- readCellType("lsec")
  counts <- d$counts; cellanno <- d$cellanno
  markers <- loadFiles("markers")$lsec
  cat("OK\n")
}

mousemodels <- fitModels(counts, cellanno, markers, modelform = mouseform, cores = 7,
  minproportion = 0.1)

betas <- map(mousemodels$models, "coef")
vars <- map(mousemodels$models, "vcov")

saveRDS(
  list(betas = betas, vars = vars, mat = mousemodels$mat),
  file.path(rdsDir, "lsec-glm-condition-nb-hvg.rds"))
