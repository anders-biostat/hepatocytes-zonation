## Files 

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

- 230-plot-marker-HEP-expr-zonation.r 231-plot-marker-LSEC-expr-zonation.r
  Assign zone pragmatically using proportion of the portal genes, scaled by quantile.
- 235-pbulk-de-conditions.r 236-pbulk-de-tables.r
  DESeq2 for bulks in LSEC and HEPs for WT vs KO
- 238-plot-ko-genes.r
  Show expression of KO genes in LSEC

- assets.r is for convenient reading the data and preliminary results
- func.r varioud function used in all scripts
<!-- - landmarks-6.r cross-checking zonation and 6 landmarkers -->

## Zone assignment

We use 21 + 33 markers from Shalev 2018 Supplementary.
Itzkovitz uses maximal expression as normalizing factor.
It can work for a single condition, but in case of decreased expression
due to k.o., RNA level can be much lower. In limiting case,
cells with a single read will be "portal" or "central", although,
in the wild type such low level of expression would indicate
an opposite phenotype.


~~Nevertheless, we followed this way for exploratory analysis with the difference,
that we normalize using 99% quantile instead of maximum.~~
We use maximum, there is no significant difference.


\Eta is a proportion of portal gene over portal + central.
```
calcEta <- function(x, markers) {
  x <- map(markers, ~ x[.x, ])
  x <- map(x, scaleByQuantile)
  x <- map(x, colSums)
  res <- x$portal / (x$portal + x$central)
}
```

## comparison of conditions

We fit glm with nb for every mouse.
Then test contrasts as mean(DKO) - mean(WT) and mean(CV-PV)_{genotype}.

## ligand - receptor

We use cellphoneDB to search for expression of interacting pairs:
we look for the genes, which are zonated in LSEC and 
which have their pair expressed in HEPs.
