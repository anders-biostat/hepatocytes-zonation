source("code/func.r")
source("code/assets.r")

results <- readRDS(file.path(rdsDir, "pbulks-deseq.rds"))
tabledir <- file.path("results", "tables")

for (i in names(results)) {
  write.csv(results[[i]]$result,
    file.path(tabledir, sprintf("pbulk-deseq-%s.csv", i)))
}

