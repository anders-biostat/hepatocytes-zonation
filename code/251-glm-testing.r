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

## vector
models <- purrr::map(which(highgenes),
  ~ glm.nb(modelform, data = cbind(gene = counts[.x,], mat)))
betas <- map(models, coef)
vars <- map(models, vcov)
saveRDS(list(betas, vars), file.path(rdsDir, "hepa-glm-nb-hvg.rds"))

teststats <- map(betas, ~ h %*% .x)
teststats <- do.call(cbind, teststats)
colnames(teststats) <- names(highgenes)[highgenes]

png(file.path(figdir, "contrast-distribution-hepa.png"),
  width = 4, height = 4, units = "in", res = 200)
plot(density(teststats[1,]), main = "", ylim = c(0,2))
lines(density(teststats[2,]), col = 2)
legend("topright", legend = c("0.2", "0.8"), col = 1:2, lty = 1)
dev.off()

res <- data.frame(
  gene = colnames(teststats), stringsAsFactors = FALSE
)
res <- cbind(res, t(teststats))
colnames(res)[2:3] <- c(names(positions))

sdContrasts <- function(v)
  sqrt(diag(h %*% v %*% t(h)))
x <- map(vars, sdContrasts)
x <- do.call(rbind, x)
colnames(x) <- sprintf("sd-%s", names(positions))

res <- cbind(res, x)
res[["pval-0.2"]] <- pnorm(res[,2], sd = res[,4], lower.tail = TRUE)
res[["pval-0.8"]] <- pnorm(res[,3], sd = res[,5], lower.tail = TRUE)
res$wtexpr <- apply(
  counts[res$gene, dat$condition == "HEP|Wildtype"], 1, mean)

saveRDS(res, file.path(rdsDir, "glm-hepa-res.rds"))

x <- res %>%
  filter(`pval-0.2` < 1e-4, abs(`0.2`) > .7 ) %>%
  arrange(`pval-0.2` )
g <- x$gene[2]
q <- qplot(
  x = dat$etaq,
  y = counts[g,]/dat$total,
  color = dat$condition,
  size = I(.6)
) +
  stat_smooth() +
  scale_y_continuous(trans = "sqrt") +
  xlab("eta quantile") +
  ylab("fraction") +
  ggtitle(g) +
  theme(legend.position = c(.8,.8)) +
  scale_colour_discrete(name="condition") 
ggsave(
  filename = file.path(figdir, sprintf("0.2-wt-high-%s.png", g)), q,
  width = 6, height = 4, dpi = 200
)

x <- res %>%
  filter(`pval-0.8` < 1e-4, abs(`0.8`) > .7) %>%
  arrange(`pval-0.8` )
g <- x$gene[1]
g <- "gm10718"
q <- qplot(
  x = dat$etaq,
  y = counts[g,]/dat$total,
  color = dat$condition,
  size = I(.6)
) +
  stat_smooth() +
  scale_y_continuous(trans = "sqrt") +
  xlab("eta quantile") +
  ylab("fraction") +
  ggtitle(g) +
  scale_colour_discrete(name="condition") +
  theme(legend.position = c(.8,.8))
ggsave(
  filename = file.path(figdir, sprintf("0.8-wt-high-%s.png", g)), q,
  width = 6, height = 4, dpi = 200
)


a <- readRDS(file.path(rdsDir, "pbulks-deseq.rds"))$heps$result
res <- cbind(res, a[res$gene,])

x <- res %>%
  filter(
    `0.2` < -.1 ,
    pvalue > .1 
    ) %>%
  arrange(`pval-0.2`)

g <- x$gene[1]
qs <- map(1:16, function(i) {
  g <- x$gene[i]
q <- qplot(
  x = dat$etaq,
  y = counts[g,]/dat$total,
  color = dat$condition,
  alpha = I(.4),
  shape = I("."),
  size = I(.2)
) +
  stat_smooth() +
  scale_y_continuous(trans = "sqrt") +
  xlab("eta quantile") +
  ylab("fraction") +
  ggtitle(g) +
  scale_colour_discrete(name="condition") +
  theme(legend.position = c(.8,.8))
})
q <- plot_grid(plotlist = qs)
ggsave(
  filename = file.path(figdir, sprintf("missed-deseq-%s.png", g)), q,
  width = 13, height = 10, dpi = 200
)
