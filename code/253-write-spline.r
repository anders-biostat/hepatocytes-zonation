
source("code/func.r")
source("code/assets.r")
SPLINE_DEGREE <- 4
X <- createSplineMat(seq(0, 1, .01), SPLINE_DEGREE)
write.csv(X, file.path("results", "tables", "spline.csv"), row.names = FALSE)
