# 05d_global_structure.R — Global structural representations
#
# Computes structural representations across ALL gRNAs pooled together,
# providing cross-guide context. Includes:
#   - Global TDA statistics (persistence stats on pooled feature space)
#   - Global distance tree (hierarchical clustering across all loci)
#   - Global latent projection (PCA or UMAP on pooled features)
#
# These are "global" in the sense that they pool all loci regardless of guide,
# giving each locus a position within the overall off-target landscape.

# uwot only needed if use_umap = TRUE
if (requireNamespace("uwot", quietly = TRUE)) library(uwot)

# ---------------------------------------------------------------------------
# Master function: compute all global structural features
# ---------------------------------------------------------------------------

#' Compute global structural representations across all loci.
#'
#' @param feat_df Data frame with all loci and their sequence features (from 03_features.R).
#' @param feature_set Character vector, feature column names to use.
#' @param edits Numeric vector (n), total edit counts for Mapper coloring.
#' @param n_landmarks Integer, landmarks for global persistent homology.
#' @param n_landscape_samples Integer, t-values for landscape sampling.
#' @param n_pcs Integer, number of global PCA/UMAP components.
#' @param use_umap Logical, if TRUE use UMAP for global projection, else PCA.
#' @param seed Integer, random seed.
#' @return List with:
#'   $tda_global — data.frame (n rows), global TDA stats + landscapes broadcast
#'   $tree_global — data.frame (n rows), global tree features per locus
#'   $latent_global — data.frame (n rows x n_pcs cols), global projection
compute_global_structure <- function(feat_df, feature_set, edits,
                                      n_landmarks = 200, n_landscape_samples = 20,
                                      n_pcs = 16, use_umap = FALSE, seed = 123,
                                      latent_method = "pca") {
  n <- nrow(feat_df)
  X <- as.matrix(feat_df[, feature_set, drop = FALSE])

  # --- Global TDA ---
  cat("Computing global TDA features...\n")
  tda_res <- compute_tda_features(
    X = X,
    edits = edits,
    n_landmarks = n_landmarks,
    n_landscape_samples = n_landscape_samples,
    mapper_params = list(num_intervals = 15, percent_overlap = 50, num_bins = 8),
    seed = seed
  )
  tda_global <- expand_tda_global(tda_res, n)

  # --- Global distance tree ---
  cat("Computing global distance tree...\n")

  # For global tree, use a synthetic mismatch column based on min mismatch to either guide
  mm_global <- pmin(
    feat_df$mismatch_count_1732 %||% 999,
    feat_df$mismatch_count_4894 %||% 999
  )
  df_for_tree <- feat_df
  df_for_tree$global_mismatch <- mm_global

  # Use the locus with minimum global mismatch as root
  root_idx <- which.min(mm_global)[1]

  tree_res <- build_distance_tree(
    df = df_for_tree,
    mismatch_col = "global_mismatch",
    feature_cols = feature_set,
    ontarget_idx = root_idx,
    cut_heights = c(0.3, 0.5, 0.7)
  )
  tree_global <- tree_res$features

  # --- Global latent projection ---
  cat("Computing global latent projection...\n")

  if (latent_method == "hpca") {
    cat("  Using hPCA...\n")
    dr_result <- dim_reduce(df_feat = feat_df, method = "hpca", k_dr = n_pcs)
    latent <- dr_result$pcs
  } else if (use_umap && n > 100) {
    X_scaled <- scale(X)
    X_scaled[!is.finite(X_scaled)] <- 0
    cat("  Using UMAP...\n")
    latent <- uwot::umap(X_scaled, n_components = n_pcs, seed = seed,
                         n_neighbors = min(30, n - 1), min_dist = 0.1)
  } else {
    X_scaled <- scale(X)
    X_scaled[!is.finite(X_scaled)] <- 0
    cat("  Using PCA...\n")
    pca <- prcomp(X_scaled, rank. = n_pcs)
    latent <- pca$x[, seq_len(min(n_pcs, ncol(pca$x))), drop = FALSE]
  }

  colnames(latent) <- paste0("global_PC", seq_len(ncol(latent)))
  latent_global <- as.data.frame(latent)

  list(
    tda_global = tda_global,
    tree_global = tree_global,
    latent_global = latent_global
  )
}

# Null-coalescing operator (R doesn't have one natively)
`%||%` <- function(a, b) if (is.null(a)) b else a
