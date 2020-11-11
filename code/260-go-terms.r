library(purrr)
library(dplyr)
library(org.Mm.eg.db)
library(yaml)
library(clusterProfiler)

PVAL_THRES <- .01
ABS_LFC_THRES <- 0

{cat("reading...")
  selectedGO <- yaml::read_yaml("from_Ki/go-terms.yaml")

  source("code/func.r")
  source("code/assets.r")

  celltypes <- set_names(c("hep", "lsec"))
  glmtests <- map(celltypes,
    ~ read.csv(file.path("results", "tables", sprintf("%s-glm-de-test.csv", .x)),
      stringsAsFactors = FALSE))
  cat("OK\n")
}

{cat("gene set analysis...")
  ## for gene set analysis split the gene list by thresholds for p-values and lfc
  prepareGeneLists <- function(x, positions=.2, pvalthres = 1e-2, abslfc = 0) {
    x <- x %>% arrange(-l2fc) %>%
      filter(position %in% positions)
    x$gene <- capitalize(x$gene)

    gene.df <- bitr(x$gene,
      fromType = "SYMBOL",
      toType = c("ENTREZID"),
      OrgDb = org.Mm.eg.db)

    gene.df <- gene.df %>% inner_join(x, by = c("SYMBOL"="gene"))
    i <- abs(gene.df$l2fc) > abslfc & gene.df$pval < pvalthres
    topgenes <- gene.df$SYMBOL[i]
    list(gene.df = gene.df, topgenes = topgenes)
  }

  gseForPosition <- function(x, position) {
    a <- prepareGeneLists(x, position, PVAL_THRES, ABS_LFC_THRES)

    res <- enrichGO(
      gene = a$topgenes,
      universe      = a$gene.df$SYMBOL,
      keyType = "SYMBOL",
      OrgDb         = org.Mm.eg.db,
      ont           = "ALL",
      pAdjustMethod = "none",
      pvalueCutoff = 1,
      qvalueCutoff = 1)
    res
  }

  celltype <- "hep"

  geres <- map(celltypes, function(celltype) {
    x <- glmtests[[celltype]]
    positions <- unique(x$position)
    xres <- map(set_names(positions), ~gseForPosition(x, .x))
    xres
  })

  saveRDS(geres, file.path(rdsDir, "geres.rds"))
  cat("OK\n")
}

{cat("gsea .....")
  gseaForPosition <- function(x, position) {
    a <- prepareGeneLists(x, position)
    geneList <- set_names(a$gene.df$l2fc, a$gene.df$SYMBOL)
    gseres <- gseGO(
      geneList,
      OrgDb  = org.Mm.eg.db,
      keyType = "SYMBOL",
      ont  = "ALL",
      pAdjustMethod = "none",
      pvalueCutoff = 1)
    gseres
  }

  gseares <- map(celltypes, function(celltype) {
    x <- glmtests[[celltype]]
    positions <- unique(x$position)
    xres <- map(set_names(positions), ~gseForPosition(x, .x))
    xres
  })

  saveRDS(gseares, file.path(rdsDir, "gseares.rds"))
  cat("OK\n")
}

go <- names(unlist(selectedGO[[celltype]]))
