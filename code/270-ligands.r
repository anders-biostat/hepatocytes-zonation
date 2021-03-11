## Here we try to select genes, which are expressed in HEP
## and zonated in LSEC:
## the idea is to have receptor possibly able to accept a signal,
## where the angiocrine signal from LSEC is zonated.

library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RSQLite)
library(DBI)
library(Matrix)

source("code/func.r")
source("code/assets.r")
resdir <- "results/tables/"

res <- list(
  hep  = readRDS(file.path(rdsDir, "hepa-glm-mouse-nb-hvg.rds")),
  lsec = readRDS(file.path(rdsDir, "lsec-glm-mouse-nb-hvg.rds")))

edgesTests <- read.csv(
  file.path(resdir, "edges-test-hep-lsec.csv"))

getMeanExpressionHepsWT <- function() {
  heps <- readCellType("heps")
  wt.heps <- heps$cellanno$Genotype == "DoubleKO"
  x <- heps$counts[wt.heps, ]
  x <- rowMeans(t(t(x)/colSums(x)))
  x
}
meanGeneFracsHEP <- getMeanExpressionHepsWT()

con <- dbConnect(RSQLite::SQLite(), "ext/cellphone.db")
dbtables <-  dbListTables(con)
cptables <- map(set_names(dbtables), dbReadTable, conn = con)
## [1] "complex_composition_table" "complex_table"
## [3] "gene_table"                "interaction_table"
## [5] "multidata_table"           "protein_table"

annotation <- cptables$protein_table %>%
  inner_join(cptables$gene_table, by = c("id_protein" = "protein_id"))

## transform to hgnc
mouse2humanFile <- file.path(rdsDir, "mouse2human.rds")

if (!file.exists(mouse2humanFile)) {
  library("biomaRt")
  mart <- list(
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

  mouseGenes <- capitalize(edgesTests$gene)
  mouse2human <- getLDS(
    filters = "mgi_symbol",
    values = names(meanGeneFracsHEP),
    mart = mart$mouse,
    martL = mart$human,
    attributes = c("mgi_symbol"),
    attributesL = c("hgnc_symbol"),
    uniqueRows=T)
  saveRDS(mouse2human, mouse2humanFile)
}

mouse2human <- readRDS(file.path(rdsDir, "mouse2human.rds"))

## add corresponding mouse genes
annotation <- annotation %>%
  inner_join(mouse2human, by = c("gene_name" = "HGNC.symbol"))

## create a list for every mouse gene


id2gene <- function(id) {
  unique(annotation$MGI.symbol[annotation$id_protein == id])
}

## list of human genes pairs (x, y) -> [(x, y), (y, x)] | mouse ids
createMouseGeneLists <- function(inter) {
  allids <- unique(inter$multidata_1_id, inter$multidata_2_id)
  withprotein <- map(allids, id2gene) %>% map_lgl(~ length(.x) > 0)
  allids <- allids[withprotein]

  x <- map(set_names(allids), function(id) {
    c(inter$multidata_2_id[inter$multidata_1_id == id],
      inter$multidata_1_id[inter$multidata_2_id == id])
  })
  x <- map(x, unlist)
  x <- map(x, ~ map(.x, id2gene))
  x <- map(x, unlist) %>% Filter(function(x) length(x) > 0, x = .)

  ## some human genes are mapped to several mouse genes:
  ## we copy lists for this genes as well
  interactingMGI <- map(names(x), id2gene) %>%
    imap(function(mousegenes, i) {
      set_names(rep(x[i], length(mousegenes)), mousegenes)
    }) %>% do.call(what = c)
  interactingMGI
}

inter <- cptables$interaction_table
interactingMGI <- createMouseGeneLists(inter)

## lfc in zonation in LSEC vs mean WT expression in HEP
genePairs <- map(interactingMGI, function(x) data.frame(second = x)) %>%
  bind_rows(.id = "first") %>% unique %>%
  mutate(first = tolower(first), second = tolower(second))

edges <- edgesTests %>% filter(gene %in% genePairs$first)
dim(edges)
## [1] 68  6
edges <- edges %>%
  inner_join(genePairs, by = c("gene" = "first"))
edges$hepMeanFrac <- meanGeneFracsHEP[edges$second]

result <- edges %>% arrange(pval) %>% filter(celltype == "lsec", pval < 0.1, hepMeanFrac > 1e-7)
result <- result %>%
  rename(lsec.gene = gene, hep.gene = second, PV.vs.CV.l2fc = l2fc)
write.csv(result, file.path(resdir, "interactions-expressed-hep-zonated-lsec.csv"))

## list interacting genes, which are zonated in both
zonatedpairs <- edges %>%
  filter(celltype == "lsec") %>%
  inner_join(edges %>% filter(celltype == "hep"),
             by = c("gene" = "second"), suffix = c(".lsec", ".hep"))
## nothing reasonable
