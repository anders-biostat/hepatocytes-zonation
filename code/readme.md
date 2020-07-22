## Preprocessing

- 100-load-raw-counts.r  
  We get count files for every cell as a csv file and collect all of them into a single matrix. 
  Save rds files for cell annotation.
- 200-qc-totals.r  
  Plot correlation heatmap and clustering of the samples as bulks.
- 220-qc-seurat.r  
  Plot UMAP with marker gene expression and clustering.
- 221-seurat-select-heps.r, 222-seurat-select-lsec.r  
  Subset cells, which express markers of LSEC or hepatocytes,
  exclude Kupfer and stellate cells
  
## Files
- assets.r is for convenient reading the data and preliminary results
- func.r varioud function used in all scripts
- plot-eta-expr.r assign zonation and plot marker expression
- landmarks-6.r cross-checking zonation and 6 landmarkers

## Zone assignment

We use 21 + 33 markers from Shalev 2018 Supplementary.
Itzkovitz uses maximal expression as normalizing factor.
It can work for a single condition, but in case of decreased expression
due to k.o., RNA level can be much lower. In limiting case,
cells with a single read will be "portal" or "central", although,
in the wild type such low level of expression would indicate
an opposite phenotype.

Nevertheless, we followed this way for exploratory analysis with the difference,
that we normalize using 99% quantile instead of maximum.

\Eta is a proportion of portal gene over portal + central.
```
calcEta <- function(x, markers) {
  x <- map(markers, ~ x[.x, ])
  x <- map(x, normByMax)
  x <- map(x, colSums)
  res <- x$portal / (x$portal + x$central)
}
```
