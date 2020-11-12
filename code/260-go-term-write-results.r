## write csv tables from rds results of gse and gsea,
## in addition subset the GO sent by Ki
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


allgeres <- map(geres, function(o) {
  o %>% map(as.data.frame) %>% bind_rows(.id = "position")
})
allgeres <- bind_rows(allgeres, .id = "celltype") %>%
  arrange(qvalue)
allgeres$p.adjust <- NULL
write.csv(allgeres, gopath("ge-go.csv"), row.names = FALSE)

allgseares <- map(gseares, function(o) {
  o %>% map(as.data.frame) %>% bind_rows(.id = "position")
})
allgseares <- bind_rows(allgseares, .id = "celltype") %>%
  arrange(qvalues)
allgseares$p.adjust <- NULL
write.csv(allgseares, gopath("gsea-go.csv"), row.names = FALSE)

allge_deseq <- map(deseq_ge, as.data.frame) %>%
  bind_rows(.id = "celltype") %>%
  arrange(qvalue)
allge_deseq$p.adjust <- NULL
write.csv(allge_deseq, gopath("ge-deseq-go.csv"), row.names = FALSE)

allgsea_deseq <- map(deseq_gsea, as.data.frame) %>%
  bind_rows(.id = "celltype") %>%
  arrange(qvalues)
allgsea_deseq$p.adjust <- NULL
write.csv(allgsea_deseq, gopath("gsea-deseq-go.csv"), row.names = FALSE)


## subset from Ki
gotab <- go %>% imap(~data.frame(celltype = .y, ID = .x)) %>% bind_rows

allgeres %>%
  right_join(gotab) %>%
write.csv(., gopath("ge-go-ROI.csv"), row.names = FALSE)

allgseares %>%
  right_join(gotab) %>%
  write.csv(., gopath("gsea-go-ROI.csv"), row.names = FALSE)

allge_deseq %>%
  right_join(gotab) %>%
  write.csv(., gopath("ge-deseq-go-ROI.csv"), row.names = FALSE)

allgsea_deseq %>%
  right_join(gotab) %>%
  write.csv(., gopath("gsea-deseq-go-ROI.csv"), row.names = FALSE)
