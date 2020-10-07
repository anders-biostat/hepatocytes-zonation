## TODO refactor to a list of 0.2, 0.8
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

{cat("reading...")
  counts <- loadFiles("counts-heps")
  cellanno <- loadFiles("cellanno-heps")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)
  conditions <- unique(cellanno$condition)
  samples <- unique(cellanno$sample)
  hepde <- readRDS(file.path(rdsDir, "pbulks-deseq.rds"))$heps$result
  cat("OK\n")
}

{## select highly expressed genes
  propexpressed <- apply(
    counts[, cellanno$condition == "HEP|Wildtype"] > 0, 1, mean)
  plot(ecdf(propexpressed))
  highgenes <- rank(-propexpressed) < 2000
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
  names(etawide)[names(etawide) == "Mouse.ID"] <- "mouse"
  etawide
}

dat <- createHEPdata()

X <- splines::bs(dat$etaq, degree = 4, Boundary.knots = c(0,1), intercept = TRUE)
splineVars <- colnames(X) <- paste0("X", colnames(X))
mat <- as.data.frame(X)
mat$condition <- dat$condition
mat$total <- dat$total
mat$mouse <- factor(dat$mouse)
Xpred <- mat
Xpred$total <- 1e4

## fit a spline for every mouse
modelform <- as.formula(
  sprintf("gene ~ -1 + offset(log(total)) + (%s):mouse",
    paste0(splineVars, collapse = "+")))

qs <- map(rownames(hepde)[1:10], function(g) {
  try({
  m <- glm.nb(modelform, data = cbind(gene = counts[g,], mat))
  pred <- predict(m, newdata = Xpred)
  ggplot(data = cbind(dat, expr = pred) %>% arrange(etaq)) +
    geom_line(aes(x = etaq, y = expr, colour = condition, group = mouse)) +
    ggtitle(g)
  })
})
q <- plot_grid(plotlist = Filter(function(x) !is(x, "try-error"), qs))
ggsave(filename = file.path(figdir, "fit-per-mouse.png"),
  q, width = 10, height = 10, dpi = 200)


## take a gene with prominent difference and counts
g <- "gulo"
## model pro mouse and model pro condition
modelform <- as.formula(
  sprintf("gene ~ -1 + offset(log(total)) + (%s):mouse",
    paste0(splineVars, collapse = "+")))
mmouse <- glm.nb(modelform, data = cbind(gene = counts[g,], mat))

modelform <- as.formula(
  sprintf("gene ~ -1 + offset(log(total)) + (%s):condition",
    paste0(splineVars, collapse = "+")))
mcondition <- glm.nb(modelform, data = cbind(gene = counts[g,], mat))

## compare with random effects for every mouse and fixed effects for condition
## quite time consuming...
modelform <- as.formula(
  paste("gene ~ -1 + offset(log(total)) +  ",
    paste0(sprintf("%s:condition +(0+%s|mouse)", splineVars, splineVars), collapse = "+")))
mmixed <- lme4::glmer.nb(modelform, data = cbind(gene = counts[g,], mat))



## compare variances with simulations
l <- list(mixed = mmixed, condition = mcondition, mouse = mmouse)

## Indeed, using a simpler model with conditions only may result in more
## false positives, since the variance is smaller and does not account for
## biological variation.
# fixed effects pro condition
##                            Estimate Std. Error   z value      Pr(>|z|)
## X1:conditionHEP|Wildtype  -7.109897  0.1111517 -63.96572  0.000000e+00
## X1:conditionHEP|DoubleKO  -7.903542  0.1495198 -52.85949  0.000000e+00
## X2:conditionHEP|Wildtype  -9.402495  0.3612774 -26.02569 2.535786e-149
## X2:conditionHEP|DoubleKO -10.097444  0.5201289 -19.41335  5.951912e-84
## X3:conditionHEP|Wildtype  -8.025290  0.6422258 -12.49606  7.844596e-36
## X3:conditionHEP|DoubleKO -10.634433  1.0045191 -10.58659  3.438856e-26
## X4:conditionHEP|Wildtype -13.778540  0.6148102 -22.41105 3.071288e-111
## X4:conditionHEP|DoubleKO -12.961319  1.0825578 -11.97287  4.929594e-33
## X5:conditionHEP|Wildtype -11.414183  0.2856873 -39.95342  0.000000e+00
## X5:conditionHEP|DoubleKO -12.940257  0.6496367 -19.91922  2.772610e-88

# mixed effect model
## > a1
##                            Estimate Std. Error   z value      Pr(>|z|)
## X1:conditionHEP|Wildtype  -7.103366  0.1256707 -56.52366  0.000000e+00
## X1:conditionHEP|DoubleKO  -7.926666  0.1613386 -49.13064  0.000000e+00
## conditionHEP|Wildtype:X2  -9.348529  0.4174393 -22.39494 4.408787e-111
## conditionHEP|DoubleKO:X2 -10.133175  0.5549899 -18.25830  1.777282e-74
## conditionHEP|Wildtype:X3  -7.981858  0.6627147 -12.04418  2.081273e-33
## conditionHEP|DoubleKO:X3 -10.723040  1.0239517 -10.47221  1.159020e-25
## conditionHEP|Wildtype:X4 -13.776912  0.6080382 -22.65797 1.164234e-113
## conditionHEP|DoubleKO:X4 -12.856112  1.0833392 -11.86712  1.754105e-32
## conditionHEP|Wildtype:X5 -11.413403  0.2861650 -39.88399  0.000000e+00
## conditionHEP|DoubleKO:X5 -12.977334  0.6527144 -19.88210  5.814484e-88


q <- qplot(
  y = predict(mmixed, newdata = Xpred),
  x = dat$etaq,
  ylab = "prediction",
  xlab = "etaq",
  size = I(.4),
  group = dat$mouse,
  colour = dat$condition) +
  scale_colour_discrete(name = "condition")
ggsave(filename = file.path(figdir, sprintf("%s-mixed-prediction.png", g)), q,
  width = 5, height = 4, dpi = 200)

modelform <- as.formula(
  paste("sqrt(gene/total) ~ -1 +   ",
    paste0(sprintf("%s:condition +(0+%s|mouse)", splineVars, splineVars), collapse = "+")))
mmixedfrac <- lme4::lmer(modelform, data = cbind(gene = counts[g,], mat))
