
# 00_setup.R — packages and small utilities

required_pkgs <- c(
  "Biostrings","ggplot2","FactoMineR","dynamicTreeCut",
  "e1071","isotree","ineq","entropy","vegan","mclust", "NMF", "caret", "RGCCA",
  "fastICA"
)

optional_pkgs <- c("multiblock", "Rdimtools", "uwot", "dendextend",
                   "TDA", "TDAmapper", "igraph", "xgboost", "SHAPforxgboost")

quiet_attach <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (p in required_pkgs) quiet_attach(p)
for (p in optional_pkgs) {
  if (requireNamespace(p, quietly = TRUE))
    suppressPackageStartupMessages(library(p, character.only = TRUE))
}
