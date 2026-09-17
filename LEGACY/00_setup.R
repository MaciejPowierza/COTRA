
# 00_setup.R — packages and small utilities

required_pkgs <- c(
  "Biostrings","ggplot2","FactoMineR","dynamicTreeCut","multiblock",
  "e1071","isotree","ineq","entropy","vegan","mclust", "NMF", "caret", "RGCCA",
  "Rdimtools", "fastICA"
)

quiet_attach <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (p in required_pkgs) quiet_attach(p)
