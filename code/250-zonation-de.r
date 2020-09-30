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


gene1 <- "cyp2a5"

dat <- data.frame(
  cell = etawide$cell,
  genecount = counts[gene1, etawide$cell],
  total = etawide$total,
  condition = etawide$condition,
  mouse = factor(etawide$Mouse.ID),
  etaq = etawide$etaq
)
dat$gene <- dat$genecount / dat$total
## dat <- dat %>% filter(etaq != 0, etaq != 1)
dat$condition <- relevel(dat$condition, ref = "HEP|Wildtype")


binomialFormula <- function(splineVars) {
  modelform <-   cbind(genecount, total - genecount) ~ - 1
  form <- paste("~ ",
                paste(sprintf("%s:condition", splineVars),
                      collapse = "+"), " + .")
  modelform <- update(modelform, as.formula(form))
  modelform
}

nbformula <- function(splineVars) {
  # model with offset for total number
  modelform <-  genecount ~ -1  + offset(log(total))
  form <- paste("~ ",
                paste(sprintf("%s:condition", splineVars),
                      collapse = "+"), " + .")
  modelform <- update(modelform, as.formula(form))
  modelform
}

Y <- t(counts[,dat$cell])
X <- splines::bs(dat$etaq,
                 ## knots = c(.2, .4, .6, .8),
                 intercept = TRUE,
                 Boundary.knots = c(0,1),
                 df = 4)
mat <- as.data.frame(X)
splineVars <- paste0("X", names(mat))
colnames(mat) <- splineVars
mat <- cbind(dat, mat)

m <- glm.nb(nbformula(splineVars), data = mat)

png(file.path(figdir, "nb-model.png"), res = 200, width = 8, height = 8,
    units = "in")
par(mfrow = c(2,2))
plot(m)
dev.off()

mqb <- glm(binomialFormula(splineVars), data = mat,
         family = quasibinomial)
png(file.path(figdir, "quasibinomial-model.png"), res = 200, width = 8, height = 8,
    units = "in")
par(mfrow = c(2,2))
plot(mqb)
dev.off()
