library(Matrix)
library(purrr)
library(dplyr)
library(tidyr)
source("code/func.r")


{cat("reading...")
  rdsDir <- file.path("results", "rds")
  countsFile <- file.path(rdsDir, "allcouts.rds")
  counts <- readRDS(countsFile)
  meta <- readRDS(file.path(rdsDir, "meta.rds"))
  cellanno <- readRDS(file.path(rdsDir, "cellanno.rds"))
  ## we work only with hepatocytes
  heps <- cellanno$Cell.type == "HEP" & cellanno$Genotype == "DoubleKO"
  heps[is.na(heps)] <- FALSE
  cellanno <- cellanno[heps,]
  counts <- counts[, heps]
  rownames(counts) <- tolower(rownames(counts))
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)

  cat("OK\n")
}


## get parameters from Itzkovitz data
gammaParams <- getItzkevitzParams()$Gamma_params
markers <- getMarkers()
## genes we do not have:
setdiff(unlist(markers), rownames(fracs))
## [1] "cml2" "cyb5" "c1s"  "dak"
## use only genes we have
markers <- lapply(markers, function(a) intersect(a, tolower(rownames(counts))))
## Itzkovitz use fracs in matlab file, not doc'ed in publ.
mexpr <- map(markers, ~ fracs[.x,])

## 2018 publ.
## 1) normalize by gene frac max across heps
normByMax <- function(x)
  x / apply(x, 1, max)
nmexpr <- map(mexpr, normByMax)
## 2) compute eta score
eta <- calcEta(map(nmexpr, colSums))
## 2) posterior probability of zones P cells x zones
zoneLiks <- calcEtaLikelihoods(eta)
priors <- getPriors(8)
priors <- rep(1/8, 8)
posteriors <- calcPosteriors(zoneLiks, priors)
## 3) column-normalize P  -> W
weights <- posteriors / colSums(posteriors, na.rm=TRUE)
weights[is.na(weights)] <- 0
## 4) multiply expression x W
exprByZone <-  fracs %*% weights

x <- exprByZone / apply(exprByZone, 1, max)

genes_interest_hep <- c("Glul","Cyp2e1","Lgr5","Slc22a3","Slc13a3","Slc1a2",
  "Slc1a2","Cib2", "Npr2","Cyp2c37","Lhpp","Cyp2a5","Slc1a4","Cyp2c38",
  "Pcp4l1", "Arg1","Mup3","Pck1","Ass1","C2","Cyp2f2","G6pc","Uox","Igf1",
  "Pigr", "Acly","Gls2","Mfsd2a","Celsr1","Bdh2")
test <- c("Gpc1","Gas2","Sdsl","Ugt2b38","Hsd17b6","Cryl1","Serpina12","Mmd2",
  "Clec2h","Setd7","Aldh1b1")

pheatmap(
  x[tolower(genes_interest_hep),],
  cluster_rows =  F,
  cluster_cols = F,
  show_colnames = F,
  ## scale = "row",
  color = rev(viridis::viridis(300)))

pheatmap(x[unlist(markers),],
  cluster_rows =  F,
  cluster_cols = F,
  show_colnames = F,
  scale = "row",
  color = (viridis::viridis(300)))

x <- as.data.frame(as.matrix(t(counts[unlist(markers),])))
x <- x %>% tibble::rownames_to_column("cell") %>%
  gather(gene, count, -cell) %>%
  inner_join(data.frame(cell = colnames(fracs), eta = eta))
levels(x$gene) <-  unique(unlist(markers))

  ggplot(data = x) +
  geom_point(aes(x = eta, y = count)) +
  facet_wrap(~gene, scales = "free_y",)
