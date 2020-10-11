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

positions <- c("0.2" = .2, "0.5" = 0.5, "0.8" = .8)

createContrasts <- function(positions, splineVars, X, wildmice, komice) {
  h <- matrix(
    nrow = length(positions),
    ncol = length(splineVars))
  rownames(h) <- names(positions)
  colnames(h) <- splineVars
  x <- predict(X, positions)
  ## Double - Wildtype
  colmouse <- gsub(".+mouse", "", colnames(h))
  for(m in wildmice) {
    h[,colmouse == m] <- -x / length(wildmice)
  }
  for(m in komice) {
    h[,colmouse == m] <- x / length(komice)
  }
  h
}

crossgroupVar <- function(beta, positions, X, wildmice, komice) {
  x <- predict(X, positions)
  colmouse <- gsub(".+mouse", "", names(beta))
  betamat <- matrix(beta, ncol = ncol(x))
  rownames(betamat) <- colmouse[seq_len(nrow(betamat))]
  ## divide by n^2 since it is var of mean for mice groups
  vwild <- apply(betamat[wildmice,], 2, var)/ (length(wildmice))^2
  vko <- apply(betamat[komice,], 2, var) / (length(komice))^2
  (x*x)%*%(vko + vwild)
}

## compute variance from vcov and contrast matrix
cntrVar <- function(cntr, v) {
  diag(cntr %*% v %*% t(cntr))
}

## wald test of linear hypothesis
testModel <- function(beta, v, cntr, positions, X, wildmice, komice) {
  cntrvar <- cntrVar(cntr, v) + crossgroupVar(beta, positions, X, wildmice, komice)
  list(l2fc = cntr %*% beta / log(2),
    sd = sqrt(cntrvar) / log(2),
    pval = 2 * pnorm(-abs(cntr %*% beta), sd = sqrt(cntrvar)))
}

tests <- lapply(res, function(x) {
  ## get mouse annotation
  splineVars <- grep("X[0-9]", colnames(x$vars[[1]]), value = TRUE)
  mouseCondition <- unique(x$mat[c("mouse", "condition")])
  wildmice <- mouseCondition$mouse[grepl("Wild", mouseCondition$condition)]
  komice <- mouseCondition$mouse[!mouseCondition$mouse %in% wildmice]
  X <- createSplineMat(x$mat$etaq, degree = 4)
  cntr <- createContrasts(positions, splineVars, X, wildmice, komice)
  map(set_names(names(x$betas)),
    ~ testModel(x$betas[[.x]], x$vars[[.x]], cntr,
      positions, X, wildmice, komice))
})

restab <- lapply(tests, function(x) {
  x <- lapply(x, function(a) {
    a <- data.frame(a)
    a$position <- rownames(a)
    rownames(a) <- NULL
    a
  })
  x <- bind_rows(x, .id = "gene")
  x
})

for(celltype in names(restab)) {
  write.csv(restab[[celltype]],
    file.path(resdir, sprintf("%s-glm-de-test.csv", celltype)))
}
