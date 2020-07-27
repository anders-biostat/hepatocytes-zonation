source("code/func.r")
source("code/assets.r")
library(ggplot2)

figdir <- file.path("results", "figs", "qc")
cellType <- "lsec"
counts <- loadFiles(glue::glue("counts-{cellType}"))
cellanno <- loadFiles(glue::glue("cellanno-{cellType}"))
cellanno$sample <- with(cellanno, paste(condition, Mouse.ID, Experimental.Batch))

kogenes <- c("wls", "rspo3")

a <- sparse2df(t(counts[kogenes,]))
a$barcode <- rownames(a)
a <- a %>% left_join(cellanno)
stopifnot(Negate(anyNA)(a))

q <- a %>%
  pivot_longer(cols = 1:2, names_to = "gene", values_to = "count") %>%
ggplot(data = .) +
  geom_jitter(aes(x = sample, y = count), height = .2) +
  facet_wrap(~ gene, scales = "free_x") +
  coord_flip()
ggsave(figpattern("kogenes-count.png"), q,
  width = 6, height = 6, dpi = 200)
