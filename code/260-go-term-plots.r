library(purrr)
library(dplyr)
library(org.Mm.eg.db)
library(yaml)
library(clusterProfiler)
source("code/func.r")
source("code/assets.r")
## https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Interpreting_Leading_Edge
godir <- file.path("results", "tables", "go")
dir.create(godir, showWarnings = FALSE)
gopath <- function(x) file.path(godir, x)
figdir <- file.path("results", "figs", "gse")

{cat("reading...")
  selectedGO <- yaml::read_yaml("from_Ki/go-terms.yaml")
  celltypes <- set_names(c("hep", "lsec"))
  go <- map(celltypes, ~names(unlist(selectedGO[[.x]])))
  geres <- readRDS(file.path(rdsDir, "geres.rds"))
  gseares <- readRDS(file.path(rdsDir, "gseares.rds"))
  deseq_ge <- readRDS(file.path(rdsDir, "gea-deseq.rds"))
  deseq_gsea <- readRDS(file.path(rdsDir, "gsea-deseq.rds"))
  cat("OK\n")
}

figfile <- function(what, celltype, position = NULL)
  file.path(figdir, paste0(paste0(c(what, celltype, position), collapse = "-"), ".png"))

plotResult <- function(dat, what) {
  for (celltype in names(dat)) {
    for (position in names(dat[[celltype]])) {
      if (!is.null(dat[[celltype]][[position]])) {
        q <- dotplot(dat[[celltype]][[position]], showCategory = 15)
        ggsave(figfile(what, celltype, position), q, width = 11, height = 6, dpi = 200)
      }
    }
  }
}

plotResult(geres, "overrepresentation")
plotResult(gseares, "gsea")

plotResultNoPosition <- function(dat, what) {
  for (celltype in names(dat)) {
      if (!is.null(dat[[celltype]])) {
        q <- dotplot(dat[[celltype]], showCategory = 15)
        ggsave(figfile(what, celltype), q, width = 11, height = 6, dpi = 200)
      }
  }
}

plotResultNoPosition(deseq_ge, "overrepresentation-on-pseudobulk")
plotResultNoPosition(deseq_gsea, "gsea-on-pseudobulk")
