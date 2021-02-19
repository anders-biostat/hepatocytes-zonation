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

positions <- c("0.2" = .2, "0.8" = .8)

x <- res$hep
splineVars <- grep("X[0-9]", colnames(x$vars[[1]]), value = TRUE)
mouseCondition <- unique(x$mat[c("mouse", "condition")])
wildmice <- mouseCondition$mouse[grepl("Wild", mouseCondition$condition)]
komice <- mouseCondition$mouse[!mouseCondition$mouse %in% wildmice]
X <- createSplineMat(x$mat$etaq, degree = 4)


createContrastsEdges <- function(positions, X, mice) {
  ## take two positions and calculate the difference
  ## it gives us coef only for one mouse
  x <- predict(X, positions)
  x <- x["0.2", ] - x["0.8", ]
  ## now assign it to the specified group
  cntr <- numeric(length(splineVars))
  names(cntr) <- splineVars
  mice <- as.character(mice)
  for (m in mice) {
    cntr[grep(paste0("mouse", m), splineVars)] <- x / length(mice)
  }
  as.matrix(cntr, ncol = 1)
}
## for understanding, what kind of result should be
## > . + > createContrastsEdges(positions, X, wildmice)
## X1:mouse1 X1:mouse2 X1:mouse3 X1:mouse4 X1:mouse5 X1:mouse6 X1:mouse7 X2:mouse1
## 0.504     0.504     0.504     0.000     0.000     0.000     0.000     0.288
## X2:mouse2 X2:mouse3 X2:mouse4 X2:mouse5 X2:mouse6 X2:mouse7 X3:mouse1 X3:mouse2
## 0.288     0.288     0.000     0.000     0.000     0.000    -0.288    -0.288
## X3:mouse3 X3:mouse4 X3:mouse5 X3:mouse6 X3:mouse7 X4:mouse1 X4:mouse2 X4:mouse3
## -0.288     0.000     0.000     0.000     0.000    -0.504    -0.504    -0.504
## X4:mouse4 X4:mouse5 X4:mouse6 X4:mouse7
## 0.000     0.000     0.000     0.000

cntr <- createContrastsEdges(positions, X, wildmice)

cntrVar <- function(cntr, v) {
  diag(t(cntr) %*% v %*% (cntr))
}

v <-  res$hep$vars[[1]]
cntrVar(cntr, v)
