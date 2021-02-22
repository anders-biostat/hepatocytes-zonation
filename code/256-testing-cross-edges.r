## We test edge expression as following:
##
## Let beta_1, beta_2, beta_3, beta_4 are coefs for a single mouse.
## beta * X_0.8 - beta * X_0.2 = beta * (X_0.8 - X_0.2) = C * beta
## beta = mean(beta_i) mathematically, but we implement beta as a long vector of
## concatenated per mouse estimations. This results in concatenation of
## contrasts into matrix C, which we need to divide by number of mice to have
## mean, not sum, in the end.
## Test statistics is s = (beta * C)^T (C V C^T)^{-1} (C * beta),
## C * beta ~ N(0, C^T V C)
## Variance is also defined as a large matrix for all mice coefs together.
##
## We add additional variance due to population variability:
## var_add = var(beta_1) over mice, which we add to diagonal elements of V.
##
##

library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
source("code/func.r")
source("code/assets.r")
resdir <- "results/tables/"

res <- list(
  hep  = readRDS(file.path(rdsDir, "hepa-glm-mouse-nb-hvg.rds")),
  lsec = readRDS(file.path(rdsDir, "lsec-glm-mouse-nb-hvg.rds")))

## we will compare only two positions, near CV, 0.2, and PV, 0.8.
positions <- c("0.2" = .2, "0.8" = .8)
SPLINE_DEGREE <- 4


createContrastsEdges <- function(positions, X, mice, splineVars) {
  ## take two positions and calculate the difference
  ## it gives us coef only for one mouse
  x <- predict(X, positions)
  x <- x["0.8", ] - x["0.2", ]
  ## now assign it to the specified group
  cntr <- numeric(length(splineVars))
  names(cntr) <- splineVars
  mice <- as.character(mice)
  for (m in mice) {
    cntr[grep(paste0("mouse", m), splineVars)] <- x / length(mice)
  }
  as.matrix(cntr, ncol = 1)
}


cntrVar <- function(cntr, v) {
  ## diag only to make it vector(1)
  diag(t(cntr) %*% v %*% (cntr))
}

populationVar <- function(betas, cntr, mice) {
  betamat <- matrix(betas, ncol = SPLINE_DEGREE, byrow = FALSE)
  vars <- apply(betamat[mice, ], 2, var)
  cntrmat <- matrix(cntr, ncol = SPLINE_DEGREE)[mice,]
  sum(cntrmat^2 %*% vars)
}

testModel <- function(beta, v, cntr, X, mice) {
  cntrvar <- cntrVar(cntr, v) + populationVar(beta, cntr, mice)
  list(l2fc = as.vector(beta %*% cntr / log(2)),
       sd = sqrt(cntrvar) / log(2),
       pval = as.numeric(2 * pnorm(-abs(beta %*% cntr), sd = sqrt(cntrvar))))
}

celltype <- "hep"

testEdges <- function(x) {
  { # prepare cell type variables
    splineVars <- grep("X[0-9]", colnames(x$vars[[1]]), value = TRUE)
    mouseCondition <- unique(x$mat[c("mouse", "condition")])
    phenotypes <- list(
      Wildtype = mouseCondition$mouse[grepl("Wild", mouseCondition$condition)],
      DoubleKO = mouseCondition$mouse[grepl("DoubleKO", mouseCondition$condition)]
    )
    phenotypes <- map(phenotypes, as.character)
    X <- createSplineMat(x$mat$etaq, degree = 4)
  }

  tests <- list()
  for (pheno in names(phenotypes)) {
    cntr <- createContrastsEdges(positions, X, phenotypes[[pheno]], splineVars)
    tests[[pheno]] <- lapply(set_names(names(x$betas)), function(g) {
      testModel(
        x$betas[[g]],
        x$vars[[g]],
        cntr,
        phenotypes[[pheno]]
      )})
  }

  tests <- purrr::map(tests, ~ bind_rows(.x, .id = "gene"))
  tests <- dplyr::bind_rows(tests, .id = "genotype")
  tests
}

tests <- list()
tests[["hep"]] <- testEdges(res[["hep"]])
tests[["lsec"]] <- testEdges(res[["lsec"]])
tests <- bind_rows(tests, .id = "genotype")

