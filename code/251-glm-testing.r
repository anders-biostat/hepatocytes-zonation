library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "zonation-tests")
dir.create(figdir, showWarnings = FALSE)

{cat("reading...")
  counts <- loadFiles("counts-heps")
  cellanno <- loadFiles("cellanno-heps")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)
  conditions <- unique(cellanno$condition)
  samples <- unique(cellanno$sample)
  cat("OK\n")
}

{## select highly expressed genes
  propexpressed <- apply(
    counts[, cellanno$condition == "HEP|Wildtype"] > 0, 1, mean)
  plot(ecdf(propexpressed))
  highgenes <- rank(-propexpressed) < 5000
}

createHEPdata <- function() {
  conditions <- c("HEP|DoubleKO", "HEP|Wildtype")
  etadf <- map(set_names(conditions),
    function(condition) {
      cells <- cellanno$condition == condition
      makeEtaDF(cells, cellanno, totals, fracs, markers)
    })
  etawide <- map(etadf, ~ pivot_wider(.x, names_from = "gene", values_from = "frac"))
  etawide <- bind_rows(etawide, .id = "condition")
  etawide <- etawide %>%
    inner_join(cellanno, by = c(cell = "barcode", condition = "condition"))

  ## etaq is quantile among the condition
  etawide <- etawide %>%
    group_by(sample) %>%
    mutate(etaq = cume_dist(eta)) %>%
    ungroup

  stopifnot(nrow(etawide) == ncol(counts))

  etawide <- etawide[match(colnames(counts), etawide$cell),]

  etawide$condition <- relevel(factor(etawide$condition), ref = "HEP|Wildtype")
  etawide
}

dat <- createHEPdata()

X <- splines::bs(dat$etaq, df = 4, Boundary.knots = c(0,1), intercept = TRUE)
splineVars <- colnames(X) <- paste0("X", colnames(X))
mat <- as.data.frame(X)
mat$condition <- dat$condition
mat$total <- dat$total
Xpred <- mat
Xpred$total <- 1e4

modelform <- as.formula(
  sprintf("gene ~ -1 + offset(log(total)) + (%s):condition",
    paste0(splineVars, collapse = "+")))

g <-which(highgenes)[1] 
m <- glm.nb(modelform, data = cbind(gene = counts[g,], mat))
m2 <- glm.nb(gene ~ - 1 + offset(log(total)) + . - total - condition,
  data = cbind(gene = counts[g,], mat))


## wald test of linear hypothesis
positions <- c("0.2" = .2, "0.8" = .8)

parvariance <- vcov(m)
h <- matrix(
  nrow = length(positions),
  ncol = 2 * length(splineVars))
x <- predict(X, positions)
## Double - Wildtype
h[,seq_along(splineVars)*2-1] <- -x
h[,seq_along(splineVars)*2] <- x

beta <- coef(m)
h %*% beta / log(2)
pnorm(h%*%beta, sd = sqrt(diag(h %*% parvariance %*% t(h))))
