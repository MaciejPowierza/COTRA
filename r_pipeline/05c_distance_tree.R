# 05c_distance_tree.R — Distance tree construction and per-locus feature extraction
#
# Builds a hierarchical clustering dendrogram (UPGMA) on sequence/feature distance
# from the on-target root, then extracts per-locus tree features.
# This is a distance tree (not a phylogeny) — off-targets are independently
# generated cut sites, not evolutionarily related. The tree captures the
# hierarchy of sequence similarity to the on-target.

library(dendextend)

# ---------------------------------------------------------------------------
# Identify on-target root loci
# ---------------------------------------------------------------------------

#' Identify on-target loci for a given gRNA.
#'
#' @param df Data frame with mismatch_count column.
#' @param mismatch_col String, column name for mismatch count to this gRNA.
#' @return Integer vector of row indices that are on-target (mismatch == 0).
identify_ontarget <- function(df, mismatch_col) {
  which(df[[mismatch_col]] == 0)
}

# ---------------------------------------------------------------------------
# Compute distance from on-target for each locus
# ---------------------------------------------------------------------------

#' Compute sequence distance from on-target to each off-target.
#'
#' Uses the mismatch_count column (Hamming distance to gRNA+PAM) as the
#' primary distance metric. Optionally augments with flanking k-mer distance.
#'
#' @param df Data frame with mismatch_count and feature columns.
#' @param mismatch_col String, column name for mismatch count.
#' @param feature_cols Character vector, feature column names for k-mer distance.
#' @param ontarget_idx Integer, row index of the on-target (first one if multiple).
#' @return Dist object (pairwise distances among all loci, with on-target as reference).
compute_tree_distances <- function(df, mismatch_col, feature_cols = NULL, ontarget_idx) {
  n <- nrow(df)

  # Primary distance: mismatch count (Hamming to gRNA)
  mm <- df[[mismatch_col]]

  # If we have feature columns, compute Euclidean distance in feature space
  # and combine with mismatch count
  if (!is.null(feature_cols) && length(feature_cols) > 0) {
    feat <- as.matrix(df[, feature_cols, drop = FALSE])
    # Scale features
    feat_scaled <- scale(feat)
    feat_scaled[!is.finite(feat_scaled)] <- 0

    # Euclidean distance in feature space
    feat_dist <- dist(feat_scaled)

    # Weighted combination: mismatch count is the dominant axis
    # Create a distance matrix that combines mismatch count and feature distance
    mm_matrix <- outer(mm, mm, function(a, b) abs(a - b))

    # Combine: 70% mismatch distance, 30% feature distance (normalized)
    feat_dist_matrix <- as.matrix(feat_dist)
    feat_max <- max(feat_dist_matrix)
    if (feat_max > 0) feat_dist_matrix <- feat_dist_matrix / feat_max

    combined <- 0.7 * mm_matrix + 0.3 * feat_dist_matrix
    diag(combined) <- 0
    d <- as.dist(combined)
  } else {
    # Distance based on mismatch count only
    mm_matrix <- outer(mm, mm, function(a, b) abs(a - b))
    diag(mm_matrix) <- 0
    d <- as.dist(mm_matrix)
  }

  d
}

# ---------------------------------------------------------------------------
# Build distance tree and extract per-locus features
# ---------------------------------------------------------------------------

#' Build UPGMA tree and extract per-locus tree features.
#'
#' @param df Data frame with locus information.
#' @param mismatch_col String, column name for mismatch count.
#' @param feature_cols Character vector, feature column names (optional).
#' @param ontarget_idx Integer, row index of on-target root.
#' @param cut_heights Numeric vector, heights at which to cut tree for cluster IDs.
#' @return List with:
#'   $features — data.frame (n rows x ~8 cols)
#'   $dendro — dendrogram object
#'   $cophenetic_dist — named numeric vector of cophenetic distances from on-target
build_distance_tree <- function(df, mismatch_col, feature_cols = NULL,
                                 ontarget_idx = NULL, cut_heights = c(0.3, 0.5, 0.7)) {
  n <- nrow(df)

  if (is.null(ontarget_idx)) {
    ontarget_idx <- identify_ontarget(df, mismatch_col)
    if (length(ontarget_idx) == 0) {
      warning("No on-target found (mismatch_count == 0). Using locus with min mismatch as root.")
      ontarget_idx <- which.min(df[[mismatch_col]])[1]
    }
    ontarget_idx <- ontarget_idx[1]
  }

  cat("  On-target root: row", ontarget_idx, 
      " mismatch=", df[[mismatch_col]][ontarget_idx], "\n")

  # Compute distances
  d <- compute_tree_distances(df, mismatch_col, feature_cols, ontarget_idx)

  # Build UPGMA dendrogram
  cat("  Building UPGMA tree...\n")
  hc <- hclust(d, method = "average")
  dendro <- as.dendrogram(hc)

  # Cophenetic distances (distance in the tree between every pair)
  coph <- as.matrix(cophenetic(hc))

  # Cophenetic distance from on-target to each locus
  ont_coph <- coph[ontarget_idx, ]
  if (is.matrix(ont_coph)) ont_coph <- ont_coph[1, ]  # handle multiple on-targets
  ont_coph <- as.numeric(ont_coph)
  names(ont_coph) <- rownames(df)

  # Tree depth from root for each leaf
  # Depth = number of edges from root to leaf
  get_depths <- function(dend, current_depth = 0, depths = list()) {
    if (is.leaf(dend)) {
      depths[[as.character(attr(dend, "label"))]] <- current_depth
      return(depths)
    }
    for (i in seq_along(dend)) {
      depths <- get_depths(dend[[i]], current_depth + 1, depths)
    }
    depths
  }

  # Ensure dendrogram leaves are labeled with row names (as character)
  dendro_labeled <- dendrapply(dendro, function(node) {
    if (is.leaf(node)) {
      attr(node, "label") <- as.character(attr(node, "label"))
    }
    node
  })

  depths_list <- get_depths(dendro_labeled)
  # Map from original row names to local 1:n indices
  depths <- sapply(seq_len(n), function(i) {
    key <- rownames(df)[i]
    if (key %in% names(depths_list)) depths_list[[key]]
    else NA_integer_
  })

  # Cluster IDs at different cut heights
  cluster_ids <- matrix(NA_integer_, nrow = n, ncol = length(cut_heights))
  for (j in seq_along(cut_heights)) {
    clusters <- cutree(hc, h = cut_heights[j])
    cluster_ids[, j] <- clusters
  }
  colnames(cluster_ids) <- paste0("tree_cluster_h", cut_heights)

  # Branch length from parent (distance between consecutive merges)
  # Approximated by the cophenetic distance to the nearest merge point
  # For simplicity, use the height at which each leaf joins the tree
  merge_heights <- hc$height[hc$order]
  # Map back to original order
  leaf_merge_height <- rep(NA_real_, n)
  for (i in seq_along(hc$order)) {
    if (hc$order[i] <= n) {
      leaf_merge_height[hc$order[i]] <- hc$height[min(i, length(hc$height))]
    }
  }

  # Number of siblings (cluster size at first cut minus 1)
  first_cut <- cutree(hc, h = cut_heights[1])
  siblings <- sapply(seq_len(n), function(i) {
    sum(first_cut == first_cut[i]) - 1
  })

  features <- data.frame(
    tree_depth = depths,
    tree_cophenetic_from_ontarget = ont_coph,
    tree_branch_length = leaf_merge_height,
    tree_n_siblings = siblings,
    cluster_ids,
    stringsAsFactors = FALSE
  )

  # Replace NA with 0
  features[is.na(features)] <- 0

  colnames(features) <- c("tree_depth", "tree_cophenetic_from_ontarget",
                           "tree_branch_length", "tree_n_siblings",
                           paste0("tree_cluster_h", cut_heights))

  list(
    features = features,
    dendro = dendro,
    cophenetic_dist = ont_coph,
    hclust_obj = hc
  )
}
