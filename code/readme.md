
- assets.r is for convenient reading the data and preliminary results
- load-raw.r is for reading meta and counts, creates rds files
- func.r varioud function used in all scripts
- plot-eta-expr.r assign zonation and plot marker expression
- landmarks-6.r cross-checking zonation and 6 landmarkers

## Zone assignment

We use 23 + 30 markers from Shalev 2018 Supplementary.
Normalize by maximal value of every gene, then sum central and portal up.
\Eta is a proportion of portal gene over portal + central.
```
calcEta <- function(x, markers) {
  x <- map(markers, ~ x[.x, ])
  x <- map(x, normByMax)
  x <- map(x, colSums)
  res <- x$portal / (x$portal + x$central)
}
```
