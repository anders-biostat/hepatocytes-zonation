library(Matrix)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
source("code/func.r")
source("code/assets.r")
resdir <- "results/tables/"

figdir <- file.path("results", "figs", "zonation-tests")


testsres <- map(c("hep", "lsec"),
  ~ read.csv(file.path(resdir, sprintf("%s-glm-de-test.csv", .x))))

res <- list(
  hep  = readRDS(file.path(rdsDir, "hepa-glm-mouse-nb-hvg.rds")),
  lsec = readRDS(file.path(rdsDir, "lsec-glm-mouse-nb-hvg.rds")))

{## HEP
  {cat("reading...")
    d <- readCellType("heps")
    counts <- d$counts; cellanno <- d$cellanno
    mat <- res$hep$mat
    cat("OK\n")
  }

  x <- testsres[[1]]
  genes <- x %>%
    mutate(padj = p.adjust(pval)) %>%
    filter(abs(l2fc) > .5, padj < 0.05) %>%
    arrange(padj) %>%
    select(gene) %>% unique %>% head(100) %>% unlist

  x <- cbind(cell = mat$cell, sparse2df(t(counts[genes, mat$cell])))
  x <- x %>% gather(gene, value, -cell) %>%
    inner_join(mat)

  q <- ggplot(data = x) +
    geom_point(aes(x = etaq, y = value/total, colour = condition),
      shape = ".") +
    stat_smooth(aes(x = etaq, y = value/total, colour = condition)) +
    facet_wrap(~gene, scales = "free_y") +
    scale_y_continuous(trans = "sqrt") +
    ylab("frac")
  ggsave(filename = file.path(figdir, "hep-test-glm.png"),
    width = 20, height = 20, dpi = 100)

}

{## LSEC

  {cat("reading...")
    d <- readCellType("lsec")
    counts <- d$counts; cellanno <- d$cellanno
    mat <- res$lsec$mat
    cat("OK\n")
  }

  x <- testsres[[2]]
  genes <- x %>% arrange(pval) %>% select(gene) %>% unique %>% head(50) %>% unlist

  x <- cbind(cell = mat$cell, sparse2df(t(counts[genes, mat$cell])))
  x <- x %>% gather(gene, value, -cell) %>%
    inner_join(mat)

  q <- ggplot(data = x) +
    geom_point(aes(x = etaq, y = value/total, colour = condition),
      shape = ".") +
    stat_smooth(aes(x = etaq, y = value/total, colour = condition)) +
    scale_y_continuous(trans = "sqrt") +
    facet_wrap(~gene, scales = "free_y") +
    ylab("frac")

  ggsave(filename = file.path(figdir, "lsec-test-glm.png"),
    width = 20, height = 20, dpi = 100)
}
