#!/usr/bin/env Rscript
# 04_visualizations.R — Structural visualizations for both gRNAs
# Generates: persistence diagrams, dendrograms with on-target, PCA biplots with Mapper overlay
# 6 figures total (3 per guide), saved as PNG (300 DPI) + SVG
#
# Usage:
#   Rscript 04_visualizations.R

suppressPackageStartupMessages({
  library(Biostrings)
  library(TDA)
  library(TDAmapper)
  library(igraph)
  library(dendextend)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# ============ CONFIGURATION ============
# EDIT THESE PATHS to match your local setup
# Use absolute paths to avoid working-directory issues

# R library path (where you installed R packages)
.libPaths(c("/workspace/.Rlib", .libPaths()))

# Path to the TSV data file (output of 00_data_prep.py)
DATA_PATH <- "/home/kinga/COTRAs/COTRA_version_3/r_pipeline/data/FINAL_RESULTS_1732_4894_with_ATAC.tsv"

# Directory containing R pipeline scripts (for sourcing helper functions)
R_PIPELINE_DIR <- "r_pipeline"

# Output directory for visualizations
OUT_DIR <- "results/visualizations"
# ========================================

# ── Source pipeline functions ──
setwd(R_PIPELINE_DIR)
source("03_features.R")
source("05b_tda.R")
source("05c_distance_tree.R")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GRNA_SPACERS <- list(
  "1732" = list(spacer = "TGCTCGAGTGGGTCCCCGTG", pam = "AGG",
                full = "TGCTCGAGTGGGTCCCCGTGAGG"),
  "4894" = list(spacer = "GAGGACGAGATGTAAGAGGCTGG", pam = "TGG",
                full = "GAGGACGAGATGTAAGAGGCTGG")
)

# ── Load data ──
cat("Loading data...\n")
df <- read.csv(DATA_PATH, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
cat("  ", nrow(df), "loci\n")

# Identify sample columns using grep (robust against R's column name mangling)
ss <- "1732_1.breakends.noEnds.FILTERED"
se_col <- "Neg_ctrl_1h_DNA1.breakends.noEnds.FILTERED"
s_idx <- grep(paste0("^", gsub("\\.", "\\\\.", ss), "$"), names(df))
e_idx <- grep(paste0("^", gsub("\\.", "\\\\.", se_col), "$"), names(df))
if (length(s_idx) == 0) s_idx <- which(names(df) == ss)
if (length(e_idx) == 0) e_idx <- which(names(df) == se_col)
if (length(s_idx) == 0 || length(e_idx) == 0) {
  s_idx <- min(grep("1732_1\\.breakends", names(df)))
  e_idx <- max(grep("Neg_ctrl.*breakends", names(df)))
}
sample_cols <- names(df)[s_idx:e_idx]

# ── Feature extraction ──
cat("Extracting sequence features...\n")
feat_df <- feature_extraction(df, seq_col = "centered_sequence", draw_plot = FALSE)
feature_set <- get_feature_set(feat_df)
cat("  ", length(feature_set), "features\n")

# ── Helper: save figure ──
save_figure <- function(plot_obj, name, width = 8, height = 6) {
  png_path <- file.path(OUT_DIR, paste0(name, ".png"))
  svg_path <- file.path(OUT_DIR, paste0(name, ".svg"))
  ggsave(png_path, plot_obj, width = width, height = height, dpi = 300)
  ggsave(svg_path, plot_obj, width = width, height = height)
  cat("  Saved:", png_path, "and", svg_path, "\n")
}

# ── Process each guide ──
for (grna_id in c("1732", "4894")) {
  cat("\n========================================\n")
  cat("Processing gRNA", grna_id, "\n")
  cat("========================================\n")

  grna_cfg <- GRNA_SPACERS[[grna_id]]
  mm_col <- paste0("mismatch_count_", grna_id)
  grna_samples <- sample_cols[grepl(paste0("^", grna_id), sample_cols)]

  # Guide membership
  has_edits <- rowSums(feat_df[, grna_samples, drop = FALSE]) > 0
  is_ontarget <- feat_df[[mm_col]] == 0
  guide_idx <- which(has_edits | is_ontarget)
  cat("  Loci for this guide:", length(guide_idx), "\n")

  if (length(guide_idx) < 10) {
    cat("  Too few loci, skipping\n")
    next
  }

  guide_feat <- feat_df[guide_idx, , drop = FALSE]
  X_guide <- as.matrix(guide_feat[, feature_set, drop = FALSE])
  guide_edits <- rowSums(guide_feat[, grna_samples, drop = FALSE])

  # On-target index (within guide subset)
  ont_local <- which(guide_feat[[mm_col]] == 0)
  if (length(ont_local) == 0) ont_local <- which.min(guide_feat[[mm_col]])[1]
  ont_local <- ont_local[1]
  ont_cid <- guide_feat$cluster_id[ont_local]
  cat("  On-target: cluster_id =", ont_cid, "\n")

  # ── PCA ──
  cat("  Computing PCA...\n")
  X_scaled <- scale(X_guide)
  X_scaled[!is.finite(X_scaled)] <- 0
  pca_res <- prcomp(X_scaled, rank. = 2)
  pc_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    edits = guide_edits,
    cluster_id = guide_feat$cluster_id,
    is_ontarget = seq_len(nrow(guide_feat)) == ont_local
  )
  # Top-5 loading features
  loadings <- pca_res$rotation[, 1:2]
  loadings_mag <- sqrt(rowSums(loadings^2))
  top5_idx <- order(loadings_mag, decreasing = TRUE)[1:5]
  top5_names <- rownames(loadings)[top5_idx]
  loadings_df <- data.frame(
    feature = top5_names,
    PC1 = loadings[top5_idx, 1] * 5,
    PC2 = loadings[top5_idx, 2] * 5
  )

  # ── TDA: Persistence diagram ──
  cat("  Computing persistence diagram...\n")
  pers <- compute_persistence(X_guide, n_landmarks = 200, maxdimension = 1, seed = 123)
  pd <- pers$diagram
  pd_df <- data.frame(
    birth = pd[, 2],
    death = pd[, 3],
    dimension = as.factor(pd[, 1])
  )
  pd_df$persistence <- pd_df$death - pd_df$birth
  pd_df$death[is.infinite(pd_df$death)] <- pers$maxscale
  pd_df$persistence[is.infinite(pd_df$persistence)] <- pd_df$death[is.infinite(pd_df$persistence)] - pd_df$birth[is.infinite(pd_df$persistence)]

  # Top-3 most persistent features
  pd_df$label <- ""
  top3 <- order(pd_df$persistence, decreasing = TRUE)[1:3]
  for (i in seq_along(top3)) {
    pd_df$label[top3[i]] <- paste0("H", pd_df$dimension[top3[i]], " persist=", round(pd_df$persistence[top3[i]], 3))
  }

  # ── Persistence diagram figure ──
  cat("  Generating persistence diagram...\n")
  max_val <- max(c(pd_df$birth, pd_df$death), na.rm = TRUE) * 1.1

  p_persist <- ggplot(pd_df, aes(x = birth, y = death, color = dimension, size = persistence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_point(alpha = 0.7) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 10,
                              color = "black", bg.color = "white", bg.r = 0.15) +
    scale_color_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                       labels = c("0" = "H0 (components)", "1" = "H1 (loops)"),
                       name = "Dimension") +
    scale_size_continuous(range = c(1, 6), name = "Persistence") +
    coord_equal(xlim = c(0, max_val), ylim = c(0, max_val)) +
    labs(
      title = paste("Persistence diagram - gRNA", grna_id),
      x = "Birth", y = "Death"
    ) +
    theme_minimal(base_family = "sans") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  save_figure(p_persist, paste0("persistence_diagram_", grna_id), width = 7, height = 7)

  # ── Distance tree ──
  cat("  Computing distance tree...\n")
  tree_res <- build_distance_tree(
    guide_feat, mismatch_col = mm_col,
    feature_cols = feature_set,
    ontarget_idx = ont_local,
    cut_heights = c(0.3, 0.5, 0.7)
  )

  dend <- tree_res$dendro
  edit_colors <- colorRampPalette(c("#0279EE", "#E9ED4C", "#FF9400"))(100)
  ed <- guide_edits
  ed_scaled <- round((ed - min(ed)) / max(1, max(ed) - min(ed)) * 99) + 1

  leaf_order <- order.dendrogram(dend)
  leaf_colors <- edit_colors[ed_scaled[leaf_order]]

  dend_colored <- dend %>%
    set("labels_cex", 0.3) %>%
    set("labels_col", leaf_colors) %>%
    set("leaves_pch", 19) %>%
    set("leaves_col", leaf_colors) %>%
    set("leaves_cex", 0.5)

  ont_label_pos <- which(leaf_order == ont_local)

  # ── Dendrogram figure ──
  cat("  Generating dendrogram...\n")
  png(file.path(OUT_DIR, paste0("dendrogram_", grna_id, ".png")), width = 1200, height = 800, res = 150)
  par(mar = c(5, 4, 4, 8))
  plot(dend_colored,
       main = paste("Distance tree - gRNA", grna_id, "(on-target highlighted)"),
       ylab = "Distance", xlab = "Loci",
       cex.main = 1.2)
  abline(h = c(0.3, 0.5, 0.7), lty = 2, col = "gray60")
  points(ont_label_pos, 0, pch = 8, col = "red", cex = 2, lwd = 2)
  text(ont_label_pos, -0.02, "on-target", col = "red", cex = 0.7, srt = 90, adj = 1)
  legend("topright", legend = c("on-target", "low edits", "high edits"),
         pch = c(8, 19, 19), col = c("red", "#0279EE", "#FF9400"),
         cex = 0.7, bg = "white")
  dev.off()

  svg(file.path(OUT_DIR, paste0("dendrogram_", grna_id, ".svg")), width = 12, height = 8)
  par(mar = c(5, 4, 4, 8))
  plot(dend_colored,
       main = paste("Distance tree - gRNA", grna_id, "(on-target highlighted)"),
       ylab = "Distance", xlab = "Loci",
       cex.main = 1.2)
  abline(h = c(0.3, 0.5, 0.7), lty = 2, col = "gray60")
  points(ont_label_pos, 0, pch = 8, col = "red", cex = 2, lwd = 2)
  text(ont_label_pos, -0.02, "on-target", col = "red", cex = 0.7, srt = 90, adj = 1)
  legend("topright", legend = c("on-target", "low edits", "high edits"),
         pch = c(8, 19, 19), col = c("red", "#0279EE", "#FF9400"),
         cex = 0.7, bg = "white")
  dev.off()
  cat("  Saved dendrogram figures\n")

  # ── Mapper graph ──
  cat("  Computing Mapper graph...\n")
  mp <- mapper_features(
    X_guide, guide_edits,
    num_intervals = 10, percent_overlap = 50, num_bins = 8
  )
  mapper_obj <- mp$mapper_obj_raw
  if (is.null(mapper_obj)) {
    pca1 <- prcomp(scale(X_guide), rank. = 1)$x[, 1]
    mapper_obj <- mapper1D(
      distance_matrix = as.matrix(dist(scale(X_guide))),
      filter_values = pca1,
      num_intervals = 10, percent_overlap = 50, num_bins = 8
    )
  }

  # Compute Mapper node centroids in PC space
  node_centroids <- t(sapply(mapper_obj$points_in_vertex, function(idx) {
    if (length(idx) == 0) return(c(NA, NA))
    c(mean(pc_df$PC1[idx]), mean(pc_df$PC2[idx]))
  }))
  node_sizes <- sapply(mapper_obj$points_in_vertex, length)

  # Build edge data frame for overlay
  adj_mat <- mapper_obj$adjacency
  edges <- which(adj_mat == 1, arr.ind = TRUE)
  edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
  if (nrow(edges) > 0) {
    edge_df <- data.frame(
      x = node_centroids[edges[, 1], 1],
      y = node_centroids[edges[, 1], 2],
      xend = node_centroids[edges[, 2], 1],
      yend = node_centroids[edges[, 2], 2]
    )
  } else {
    edge_df <- data.frame(x = numeric(0), y = numeric(0), xend = numeric(0), yend = numeric(0))
  }

  centroid_df <- data.frame(
    PC1 = node_centroids[, 1],
    PC2 = node_centroids[, 2],
    size = node_sizes
  )
  centroid_df <- centroid_df[!is.na(centroid_df$PC1), ]

  # ── PCA biplot with Mapper overlay ──
  cat("  Generating PCA biplot with Mapper overlay...\n")

  p_biplot <- ggplot() +
    geom_segment(data = edge_df, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "gray50", alpha = 0.4, linewidth = 0.5) +
    geom_point(data = centroid_df, aes(x = PC1, y = PC2, size = size),
               shape = 1, color = "gray30", alpha = 0.6) +
    scale_size_continuous(range = c(2, 8), name = "Mapper\nvertex size") +
    geom_point(data = pc_df, aes(x = PC1, y = PC2, color = edits),
               alpha = 0.5, size = 1.5) +
    scale_color_viridis_c(name = "Edit count", trans = "log1p") +
    geom_point(data = pc_df[pc_df$is_ontarget, ], aes(x = PC1, y = PC2),
               shape = 17, color = "red", size = 4) +
    ggrepel::geom_text_repel(
      data = pc_df[pc_df$is_ontarget, ], aes(x = PC1, y = PC2, label = "on-target"),
      size = 3, color = "red", nudge_y = 0.5, max.overlaps = 5
    ) +
    geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.15, "cm")), color = "darkred", alpha = 0.7) +
    ggrepel::geom_text_repel(
      data = loadings_df, aes(x = PC1, y = PC2, label = feature),
      size = 2.5, color = "darkred", max.overlaps = 5
    ) +
    labs(
      title = paste("PCA biplot with Mapper overlay - gRNA", grna_id),
      x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
    ) +
    theme_minimal(base_family = "sans") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  save_figure(p_biplot, paste0("pca_biplot_mapper_", grna_id), width = 9, height = 7)
}

cat("\n========================================\n")
cat("All visualizations complete.\n")
cat("Saved to:", OUT_DIR, "\n")
cat("========================================\n")
