## here we explore modeling of gene expression along the assigned position
## within portal-central dimension.
library(mgcv)
library(locfit)
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
  cat("OK\n")
}
conditions <- unique(cellanno$condition)

samples <- unique(cellanno$sample)

## we select only hepatocytes to create a model with spline along the position
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



createMat <- function(etawide, gene1) {
  dat <- data.frame(
    cell = etawide$cell,
    genecount = counts[gene1, etawide$cell],
    total = etawide$total,
    condition = etawide$condition,
    mouse = factor(etawide$Mouse.ID),
    etaq = etawide$etaq
  )
  dat$gene <- dat$genecount / dat$total
  dat$condition <- relevel(dat$condition, ref = "HEP|Wildtype")

  Y <- t(counts[,dat$cell])
  X <- splines::bs(dat$etaq,
                   intercept = TRUE,
                   Boundary.knots = c(0,1),
                   df = 4)
  mat <- as.data.frame(X)
  splineVars <- paste0("X", names(mat))
  colnames(mat) <- splineVars
  mat <- cbind(dat, mat)
  list(mat = mat, splineVars = splineVars)
}

nbformula <- function(splineVars, useOffset = TRUE) {
  modelform <- genecount ~ -1
                                        # model with offset for total number
  if (useOffset)
    modelform <- update(modelform, ~ . + offset(log(total)))
  form <- paste("~ ",
                paste(sprintf("%s:condition", splineVars),
                      collapse = "+"), " + .")
  modelform <- update(modelform, as.formula(form))
  modelform
}

binomialFormula <- function(splineVars) {
  modelform <-   cbind(genecount, total - genecount) ~ - 1
  form <- paste("~ ",
                paste(sprintf("%s:condition", splineVars),
                      collapse = "+"), " + .")
  modelform <- update(modelform, as.formula(form))
  modelform
}

## compare the models for a couple of genes, which were identified as de by deseq
## spaghetti code follows

gene1 <- "aqp9"
gene1 <- "cyp2a5"

testGenes <- c("aqp9", "cyp2a5")


{# neg binomial
  for (g in testGenes) {
    a <- createMat(etawide, g)
    mat <- a$mat
    splineVars <- a$splineVars
    rm(a)
    m <- glm.nb(nbformula(splineVars), data = mat)

    png(file.path(figdir, sprintf("nb-model-%s.png", g)),
        res = 200, width = 8, height = 8, units = "in")
    par(mfrow = c(2,2))
    plot(m)
    dev.off()

    png(file.path(figdir, sprintf("nb-model-rstandard-%s.png", g)),
                  res = 200, width = 8, height = 8,
                  units = "in")
        plot(rstandard(m), x = dat$etaq, xlab = "eta quantile", ylab = "resid. st.")
        dev.off()

    }
}


{ ## qbinomial fit
  for (g in testGenes) {
    a <- createMat(etawide, g)
    mat <- a$mat
    splineVars <- a$splineVars
    rm(a)

    mqb <- glm(binomialFormula(splineVars), data = mat,
               family = quasibinomial)
    png(
      file.path(figdir, sprintf("quasibinomial-model-%s.png", g)),
      res = 200, width = 8, height = 8,
      units = "in")
    par(mfrow = c(2,2))
    plot(mqb)
    dev.off()
    ## since it is count data, we see negative residuals for zero counts, and there are a lot.
    png(file.path(figdir, sprintf("qb-model-rstandard-%s.png", g)),
        res = 200, width = 8, height = 8, units = "in")
    plot(rstandard(mqb), x = dat$etaq, xlab = "eta quantile", ylab = "resid. st.",
         col = factor(dat$genecount == 0))
    dev.off()
  }
}

{## what if we fit only for a single condition
  modelform <- update(modelform, genecount ~ -1 + offset(log(total)))
  form <- paste("~ ", paste(sprintf("%s", splineVars), collapse = "+"), " + .")
  modelform <- update(modelform, as.formula(form))
  m <- glm.nb(modelform, data = mat[mat$condition == "HEP|Wildtype", ])
  png(file.path(figdir, "nb-model-wildtype.png"), res = 200, width = 8, height = 8,
      units = "in")
  par(mfrow = c(2,2))
  plot(m)
  dev.off()
}

