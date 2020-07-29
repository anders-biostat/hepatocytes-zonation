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
  OK()
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
allgenes <- unique(c(unlist(shalevList), unlist(hepatoList)))
a <- setdiff(allgenes, rownames(counts))

x <- map(a, ~ agrep(.x, x = rownames(counts), max.distance = 1, value = TRUE))
names(x) <- a
map_int(x, length) %>% sort 

shalevList <- map(shalevList, ~ .x[.x %in% rownames(counts)])
hepatoList <- map(hepatoList, ~ .x[.x %in% rownames(counts)])
saveRDS(list(shalev = shalevList, ki = hepatoList),
        file.path(rdsDir, "functional-genes.rds"))

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
symbolattrs <- c("external_gene_name", "hgnc_symbol", "mgi_symbol", "uniprot_gn_symbol")

martsymbols <- map(symbolattrs[1], ~getBM(
  attributes = symbolattrs,
  filters = .x,
  values = corrections[!success,]$SYMBOL,
  mart = ensembl))
