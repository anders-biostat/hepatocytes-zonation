## here we use the results of the DE from the pseudobulks for tests
## between genotypes and filter the genes which are located in the cellphonedb.
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

pbulkres <- readRDS("./results/rds/pbulks-deseq.rds")

## first get the mapping from biomaRt
library("biomaRt")
createMouse2HumanMapping <- function(mouseGenes) {
  mart <- list(
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

  mouse2human <- getLDS(
    filters = "mgi_symbol",
    values = mouseGenes,
    mart = mart$mouse,
    martL = mart$human,
    attributes = c("mgi_symbol"),
    attributesL = c("hgnc_symbol"),
    uniqueRows=T)
  mouse2human
}

mouse2humanFile <- file.path(rdsDir, "mouse2human-biomart.rds")
mouseGenes <- map(pbulkres, "result") %>% map(rownames) %>% unlist %>% unique
mouse2human <- createMouse2HumanMapping(mouseGenes)
saveRDS(mouse2human, mouse2humanFile)

## add mapping information to the cellphonedb database
mouse2human <- readRDS(mouse2humanFile)
con <- dbConnect(RSQLite::SQLite(), "ext/cellphone.db")
dbtables <-  dbListTables(con)
cptables <- map(set_names(dbtables), dbReadTable, conn = con)
annotation <- cptables$protein_table %>%
  inner_join(cptables$gene_table, by = c("id_protein" = "protein_id"))
annotation <- annotation %>%
  inner_join(mouse2human, by = c("gene_name" = "HGNC.symbol"))

inter <- cptables$interaction_table
interactingMGI <- createMouseGeneLists(inter)

annotation <- annotation %>% inner_join(cptables$multidata_table,
                                        by = c("protein_multidata_id" = "id_multidata"))

genePairs <- map(interactingMGI, function(x) data.frame(second = x)) %>%
  bind_rows(.id = "first") %>% unique

saveRDS(genePairs, file.path(rdsDir, "gene-pairs-cellphone-mouse-biomart.rds"))

## and now finally intersect the DE results and the cellphonedb orthologes
genePairs <- readRDS(file.path(rdsDir, "gene-pairs-cellphone-mouse-biomart.rds"))

depairs <- pbulkres %>%
  map(function(x) {
    x$result %>% tibble::rownames_to_column(var = "gene") %>%
      mutate(gene = capitalize(gene)) %>%
      filter(gene %in% unlist(genePairs))
  })


## add info about the pairs
g <- genePairs
g$first <- genePairs$second
g$second <- genePairs$first
g <- unique(rbind(genePairs, g))


decellphone <- depairs$heps %>%
  inner_join(g, by = c(gene = "first")) %>%
  inner_join(depairs$lsec, by = c(second = "gene"),
             suffix = c(".hep", ".lsec")) %>%
  rename(lsec.gene = "second", hep.gene = "gene")

write.csv(decellphone,
          file = file.path(resdir, "pbulks-wt-ko-cellphonedb-pairs.csv"),
          row.names = FALSE)
