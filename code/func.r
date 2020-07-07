
library(rmatio)

getItzkevitzParams <- function() {
  zoneParamsFile <- "from_Ki/Itzkovitz zonation metlab/matlab_code//Zonation_params.mat"
  zonationParams <- read.mat(zoneParamsFile)
  zonationParams
}



getMarkers <- function() {
  zonationParams <- getItzkevitzParams()
  markers <- list(
    central = zonationParams$genes_cv,
    portal = zonationParams$genes_pn)
  markers <- lapply(markers, unlist)
  markers
}
