# 05b_tda.R — Topological Data Analysis feature extraction
#
# Provides persistent homology (landmark Rips), persistence landscapes,
# persistence statistics, and Mapper graph features.
# Used for per-guide and global structural representations.

library(TDA)
library(TDAmapper)
library(igraph)
library(ineq)
library(entropy)

# ---------------------------------------------------------------------------
# Persistent homology via landmark subsampling
# ---------------------------------------------------------------------------

#' Compute persistent homology using k-means landmarks.
#'
#' @param X Numeric matrix (n x p), feature matrix for one guide or pooled.
#' @param n_landmarks Integer, number of k-means landmarks (default 200).
#' @param maxdimension Integer, max homology dimension (0=components, 1=loops).
#' @param seed Integer, random seed for k-means reproducibility.
#' @return List with: diagram (TDA persistence diagram), landmarks (matrix),
#'   maxscale (numeric, diameter of landmark set).
compute_persistence <- function(X, n_landmarks = 200, maxdimension = 1, seed = 123) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  kk <- min(n_landmarks, n)

  km <- kmeans(X, centers = kk, iter.max = 50, nstart = 3)
  landmarks <- km$centers

  d <- dist(landmarks)
  maxscale <- max(d)

  diag <- ripsDiag(
    X = landmarks,
    maxdimension = maxdimension,
    maxscale = maxscale,
    library = "GUDHI"
  )

  list(diagram = diag[["diagram"]], landmarks = landmarks, maxscale = maxscale)
}

# ---------------------------------------------------------------------------
# Persistence statistics — compact scalar summaries of a persistence diagram
# ---------------------------------------------------------------------------

#' Extract scalar persistence statistics from a persistence diagram.
#'
#' @param pd Persistence diagram (matrix from TDA::ripsDiag), columns: dimension, birth, death.
#' @param maxdimension Integer, max dimension computed.
#' @return Named numeric vector (~8-10 values).
persistence_statistics <- function(pd, maxdimension = 1) {
  # Access columns case-insensitively
  dim_col <- grep("^dimension$", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  death_col <- grep("^death$", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  birth_col <- grep("^birth$", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  stats <- c()

  for (dim in 0:maxdimension) {
    sub <- pd[pd[, dim_col] == dim, , drop = FALSE]
    prefix <- paste0("H", dim, "_")

    if (nrow(sub) == 0) {
      stats[paste0(prefix, "betti")] <- 0
      stats[paste0(prefix, "max_persist")] <- 0
      stats[paste0(prefix, "entropy")] <- 0
      stats[paste0(prefix, "wasserstein")] <- 0
      next
    }

    persistences <- sub[, death_col] - sub[, birth_col]
    # Exclude the infinite H0 feature (death = maxscale, often Inf or very large)
    finite <- is.finite(persistences)
    persistences_finite <- persistences[finite]

    # Betti number: number of features with death > birth (positive persistence)
    stats[paste0(prefix, "betti")] <- sum(persistences > 1e-10)

    # Max persistence
    stats[paste0(prefix, "max_persist")] <- if (length(persistences_finite) > 0) max(persistences_finite) else 0

    # Persistence entropy (Shannon entropy of normalized persistences)
    if (length(persistences_finite) > 0 && sum(persistences_finite) > 0) {
      probs <- persistences_finite / sum(persistences_finite)
      stats[paste0(prefix, "entropy")] <- entropy::entropy(probs, method = "ML")
    } else {
      stats[paste0(prefix, "entropy")] <- 0
    }

    # Wasserstein amplitude (1-Wasserstein distance from diagram to empty diagram)
    stats[paste0(prefix, "wasserstein")] <- if (length(persistences_finite) > 0) sum(abs(persistences_finite)) else 0
  }

  stats
}

# ---------------------------------------------------------------------------
# Persistence landscapes — functional representation vectorized by sampling
# ---------------------------------------------------------------------------

#' Compute persistence landscape features by sampling at fixed t-values.
#'
#' @param pd Persistence diagram (matrix from TDA::ripsDiag).
#' @param dimension Integer, which homology dimension.
#' @param n_samples Integer, number of t-values to sample (default 20).
#' @param KK Integer, which landscape function (1 = top, 2 = second).
#' @return Named numeric vector of length n_samples.
landscape_features <- function(pd, dimension, n_samples = 20, KK = 1) {
  # Access columns case-insensitively without modifying pd
  dim_col <- grep("^dimension$", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  death_col <- grep("^death$", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  birth_col <- grep("^birth$", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  sub <- pd[pd[, dim_col] == dimension, , drop = FALSE]
  persistences <- sub[, death_col] - sub[, birth_col]
  finite <- is.finite(persistences)
  max_death <- max(sub[finite, death_col], na.rm = TRUE)

  if (max_death <= 0 || sum(finite) == 0) {
    return(setNames(rep(0, n_samples), paste0("land_H", dimension, "_k", KK, "_", seq_len(n_samples))))
  }

  tseq <- seq(0, max_death, length.out = n_samples)
  land <- TDA::landscape(pd, dimension = dimension, KK = KK, tseq = tseq)

  setNames(as.numeric(land), paste0("land_H", dimension, "_k", KK, "_", seq_len(n_samples)))
}

# ---------------------------------------------------------------------------
# Mapper graph features — per-locus features from Mapper graph membership
# ---------------------------------------------------------------------------

#' Compute Mapper graph and extract per-locus features.
#'
#' @param X Numeric matrix (n x p), feature matrix.
#' @param edits Numeric vector (n), edit counts for coloring nodes.
#' @param num_intervals Integer, number of intervals for Mapper filter.
#' @param percent_overlap Numeric, overlap fraction between intervals.
#' @param num_bins Integer, bins for clustering within intervals.
#' @return List with: features (data.frame, n rows x ~4 cols), graph (igraph object).
mapper_features <- function(X, edits, num_intervals = 10, percent_overlap = 50,
                             num_bins = 8) {
  X <- as.matrix(X)
  X_scaled <- scale(X)
  X_scaled[!is.finite(X_scaled)] <- 0

  # Use first PCA component as filter (lens)
  pca <- prcomp(X_scaled, rank. = 1)
  filter_vals <- pca$x[, 1]

  # Compute distance matrix
  d <- dist(X_scaled)

  # Run Mapper
  mapper_obj <- mapper1D(
    distance_matrix = d,
    filter_values = filter_vals,
    num_intervals = num_intervals,
    percent_overlap = percent_overlap,
    num_bins_when_clustering = num_bins
  )

  # Build graph
  g <- graph.adjacency(mapper_obj$adjacency, mode = "undirected")

  # Per-locus features from Mapper node membership
  n <- nrow(X)
  node_of <- rep(NA_integer_, n)
  for (v in seq_along(mapper_obj$points_in_vertex)) {
    idxs <- mapper_obj$points_in_vertex[[v]]
    node_of[idxs] <- v
  }

  # Node-level statistics
  node_degree <- degree(g)
  node_size <- sapply(mapper_obj$points_in_vertex, length)
  node_mean_edits <- sapply(mapper_obj$points_in_vertex, function(idx) {
    if (length(idx) == 0) return(NA_real_)
    mean(edits[idx], na.rm = TRUE)
  })

  # Connected component sizes
  comp <- components(g)
  comp_size <- comp$csize

  features <- data.frame(
    mapper_node_degree = node_degree[node_of],
    mapper_node_size = node_size[node_of],
    mapper_node_mean_edits = node_mean_edits[node_of],
    mapper_comp_size = comp_size[comp$membership[node_of]]
  )

  # Replace NA (loci not in any node) with 0
  features[is.na(features)] <- 0

  list(features = features, graph = g, mapper_obj = mapper_obj)
}

# ---------------------------------------------------------------------------
# Master function: compute all TDA features for a set of loci
# ---------------------------------------------------------------------------

#' Compute all TDA features for a feature matrix.
#'
#' @param X Numeric matrix (n x p), feature matrix for one guide or pooled.
#' @param edits Numeric vector (n), edit counts for Mapper coloring.
#' @param n_landmarks Integer, number of landmarks for persistent homology.
#' @param n_landscape_samples Integer, t-values per landscape.
#' @param mapper_params List with num_intervals, percent_overlap, num_bins.
#' @param seed Integer, random seed.
#' @return List with:
#'   $global_stats — named numeric vector (persistence statistics, ~8-10 values)
#'   $landscape_H0 — named numeric vector (landscape features, n_landscape_samples values)
#'   $landscape_H1 — named numeric vector (landscape features, n_landscape_samples values)
#'   $mapper_features — data.frame (n rows x ~4 cols, per-locus)
#'   $persistence — raw persistence diagram and landmarks (for debugging)
compute_tda_features <- function(X, edits, n_landmarks = 200, n_landscape_samples = 20,
                                  mapper_params = list(num_intervals = 10,
                                                        percent_overlap = 50,
                                                        num_bins = 8),
                                  seed = 123) {
  cat("  Computing persistent homology (landmark Rips)...\n")
  pers <- compute_persistence(X, n_landmarks = n_landmarks, maxdimension = 1, seed = seed)
  pd <- pers$diagram

  cat("  Extracting persistence statistics...\n")
  global_stats <- persistence_statistics(pd, maxdimension = 1)

  cat("  Extracting persistence landscapes (H0, H1)...\n")
  land_H0 <- landscape_features(pd, dimension = 0, n_samples = n_landscape_samples, KK = 1)
  land_H1 <- landscape_features(pd, dimension = 1, n_samples = n_landscape_samples, KK = 1)

  cat("  Computing Mapper graph features...\n")
  mp <- mapper_features(
    X, edits,
    num_intervals = mapper_params$num_intervals,
    percent_overlap = mapper_params$percent_overlap,
    num_bins = mapper_params$num_bins
  )

  list(
    global_stats = global_stats,
    landscape_H0 = land_H0,
    landscape_H1 = land_H1,
    mapper_features = mp$features,
    persistence = list(diagram = pd, landmarks = pers$landmarks, maxscale = pers$maxscale),
    mapper_graph = mp$graph,
    mapper_obj = mp$mapper_obj
  )
}

# ---------------------------------------------------------------------------
# Broadcast global TDA stats to all loci (same values for every locus in a guide)
# ---------------------------------------------------------------------------

#' Expand global TDA statistics to one row per locus.
#'
#' @param tda_result Result from compute_tda_features().
#' @param n Integer, number of loci.
#' @return Data.frame (n rows) with global stats + landscape features broadcast.
expand_tda_global <- function(tda_result, n) {
  global_vec <- c(tda_result$global_stats, tda_result$landscape_H0, tda_result$landscape_H1)
  as.data.frame(matrix(rep(global_vec, each = n), nrow = n, byrow = FALSE,
                       dimnames = list(NULL, names(global_vec))))
}
