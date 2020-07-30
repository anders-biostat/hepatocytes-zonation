## cleaning list of genes sent by Ki
library(Matrix)
library(purrr)
library(dplyr)
library(tidyr)
library(readxl)
source("code/func.r")
source("code/assets.r")

figdir <- file.path("results", "figs", "functional-genes")

geneFile <- "from_Ki/24072020_Liver funcational gene lists.xlsx"

xl2lists <- function(x) {
  x <- as.list(x)
  map(x, ~ Filter(Negate(is.na), .x))
}

{ cat("reading...")
  counts <- loadFiles("counts-heps")
  cellanno <- loadFiles("cellanno-heps")
  totals <- colSums(counts)
  fracs <- t(t(counts) / totals)
  markers <- loadFiles("markers")$itzkevitz
  markers <- add6LandMarks(markers)
  cat("OK\n")
}

shalevList <- read_excel(geneFile, sheet = 1)
shalevList <- map(xl2lists(shalevList), tolower)

## Some of cells include many genes, Ki explained, that it is a mistake and
## one should use only the first gene.
## "hsd3b4;hsd3b5" -> "hsd3b4"
shalevList <- map(shalevList, ~gsub(";.+", "", .x))

hepatoList <- map(xl2lists(read_excel(geneFile, sheet = 2, skip = 2)), tolower)
hepatoList <- map(hepatoList, ~gsub(";.+", "", .x))

allgenes <- unique(c(unlist(shalevList), unlist(hepatoList)))
a <- setdiff(allgenes, rownames(counts))

capitalize <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}

library(org.Mm.eg.db)
corrections <- select(org.Mm.eg.db,
       columns = c("SYMBOL", "ALIAS"),
       keytype = "ALIAS",
       keys = capitalize(a))

corrections <- corrections %>% filter(!is.na(SYMBOL))

success <- tolower(corrections$SYMBOL) %in% rownames(counts)
## correction map
alias2symbol <- tolower(corrections$SYMBOL)
names(alias2symbol) <- tolower(corrections$ALIAS)
alias2symbol <- alias2symbol[success]

correctGenes <- function(x, corrmap) {
  unlist(map_if(x, ~ .x %in% names(corrmap), ~ corrmap[.x]))
}

shalevList <- map(shalevList, correctGenes, alias2symbol)
hepatoList <- map(hepatoList, correctGenes, alias2symbol)

## read corrections from Ki
kilist <- unlist(yaml::read_yaml(file.path("from_Ki", "gene-aliases-ki.yaml")))
shalevList <- map(shalevList, correctGenes, kilist)
hepatoList <- map(hepatoList, correctGenes, kilist)

allgenes <- unique(c(unlist(shalevList), unlist(hepatoList)))
a <- setdiff(allgenes, rownames(counts))
## still 32 genes not found, try to recover using biomart

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
symbolattrs <- c("external_gene_name", "hgnc_symbol", "mgi_symbol", "uniprot_gn_symbol",
  "external_synonym")

martsymbols <- map(symbolattrs[5], ~getBM(
  attributes = symbolattrs,
  filters = .x,
  values = a,
  mart = ensembl))

x <- martsymbols[[1]]
x <- unique(x[c("mgi_symbol", "external_synonym")])
x <- map_df(x, tolower)
i <- x$mgi_symbol %in% rownames(counts)

alias2symbol <- set_names(x$mgi_symbol[i], x$external_synonym[i])
## there are some ambiguous pairs (asp, cs1)
## 1700052n19rik: armt1
## 2700094k13rik: selenoh
## 4930420k17rik: tmem243
## 4933433p14rik: gskip
## 8430408g22rik: depp1
## 8430410k20rik: msantd4
## asp: ropn1l
## asp: a
## asp: tmprss11d
## c230081a13rik: peak1
## cs1: itprid2
## cs1: slamf7
## lat1: slc7a5
## mrp4: abcc4
## pox1: prodh2
length(setdiff(a, names(alias2symbol)))
## 20 genes still missing
alias2symbol <- alias2symbol[! alias2symbol %in% c("a", "tmprss11d", "itprid2")]

shalevList <- map(shalevList, correctGenes, alias2symbol)
hepatoList <- map(hepatoList, correctGenes, alias2symbol)

shalevList <- map(shalevList, ~ .x[.x %in% rownames(counts)])
hepatoList <- map(hepatoList, ~ .x[.x %in% rownames(counts)])
saveRDS(list(shalev = shalevList, ki = hepatoList),
  file.path(rdsDir, "functional-genes.rds"))
