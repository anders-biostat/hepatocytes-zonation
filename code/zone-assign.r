library(Matrix)
library(purrr)
library(dplyr)
library(tidyr)
library(pheatmap)
source("code/func.r")
figdir <- file.path("results", "figs")
figpath <- function(x) file.path(figdir, x)


{cat("reading...")
  rdsDir <- file.path("results", "rds")
  countsFile <- file.path(rdsDir, "allcouts.rds")
  counts <- readRDS(countsFile)
  meta <- readRDS(file.path(rdsDir, "meta.rds"))
  cellanno <- readRDS(file.path(rdsDir, "cellanno.rds"))
  rownames(counts) <- tolower(rownames(counts))
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)

  markers <- getZoneMarkers()
  cat("OK\n")
}

## genes we do not have:
setdiff(unlist(markers$itzkevitz), rownames(counts))
## [1] "cml2" "cyb5" "c1s"  "dak"
## use only genes we have
markers$itzkevitz <- lapply(markers$itzkevitz, function(a) intersect(a, rownames(counts)))

makeEtaDF <- function(idx, cellanno, counts, fracs, markers) {
  cellanno <- cellanno[idx, ]
  counts <- counts[, idx]
  fracs <- fracs[, idx]
  totals <- totals[idx]
  mexpr <- map(markers, ~ fracs[.x, ])

  ## just ratio of portal/(central + portal)
  eta <- calcEta(map(mexpr,  colSums))
  d <- do.call(rbind, map(markers, ~ fracs[.x,]))
  d <- sparse2df(t(d))
  d$eta <- eta
  d$total <- totals
  d <- d %>% gather(gene, frac, -eta, -total)
  d$gene <- factor(d$gene, levels = unlist(markers))
  d
}


## select condition
sets <- list(
  "HEP-WT" = cellanno$Cell.type == "HEP" & cellanno$Genotype != "DoubleKO",
  "HEP-DKO" = cellanno$Cell.type == "HEP" & cellanno$Genotype == "DoubleKO",
  "LSEC-WT" = cellanno$Cell.type == "LSEC" & cellanno$Genotype != "DoubleKO",
  "LSEC-WT" = cellanno$Cell.type == "LSEC" & cellanno$Genotype != "DoubleKO"
  )

for(i in names(sets)[1:2]) {
  d <- makeEtaDF(sets[[i]], cellanno, counts, fracs, markers$itzkevitz)
  ## plot expression along eta and smoothing
  q <- d %>%
    ggplot(data= .) +
    geom_smooth(aes(x = log(eta), y = frac), method = "loess") +
    ## geom_point(aes(x = log(eta), y = frac), shape = ".") +
    facet_wrap(~gene, scales = "free_y")

  ggsave(filename = figpath(sprintf("marker-frac-smooth-eta-%s.png", i)),
    q, width = 12, height = 8, dpi = 200)

  q <- d %>%
    ggplot(data= .) +
    geom_smooth(aes(x = log(eta), y = frac), method = "loess") +
    geom_point(aes(x = log(eta), y = frac), shape = ".") +
    facet_wrap(~gene, scales = "free_y")
  ggsave(filename = figpath(sprintf("marker-frac-smooth-points-eta-%s.png", i)),
    q, width = 12, height = 8, dpi = 200)

}


## left overs of ki's script
genes_interest_hep <- c("Glul","Cyp2e1","Lgr5","Slc22a3","Slc13a3","Slc1a2",
  "Slc1a2","Cib2", "Npr2","Cyp2c37","Lhpp","Cyp2a5","Slc1a4","Cyp2c38",
  "Pcp4l1", "Arg1","Mup3","Pck1","Ass1","C2","Cyp2f2","G6pc","Uox","Igf1",
  "Pigr", "Acly","Gls2","Mfsd2a","Celsr1","Bdh2")
test <- c("Gpc1","Gas2","Sdsl","Ugt2b38","Hsd17b6","Cryl1","Serpina12","Mmd2",
  "Clec2h","Setd7","Aldh1b1")



