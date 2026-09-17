# 09b_pair_scorer_viz.R — Visualization for pair scorer results
#
# Produces:
#   1. Predicted vs actual scatter (test set)
#   2. Per-guide performance comparison
#   3. SHAP feature importance by feature group
#   4. Feature group ablation bar chart
#   5. Mapper graph visualization (optional, from TDA module)

library(ggplot2)
library(SHAPforxgboost)

# ---------------------------------------------------------------------------
# 1. Predicted vs actual scatter
# ---------------------------------------------------------------------------

#' Plot predicted vs actual log_edit_count for test set.
#'
#' @param predictions Data frame from train_pair_scorer$predictions.
#' @param model_name String, for title.
#' @param out_path String, path to save figure (NULL = return plot object).
#' @return ggplot object (invisibly).
plot_pred_vs_actual <- function(predictions, model_name = "Full model",
                                 out_path = NULL) {
  df <- predictions
  df$grna <- ifelse(df$grna_1732 == 1, "gRNA 1732",
             ifelse(df$grna_4894 == 1, "gRNA 4894", "No gRNA"))

  p <- ggplot(df, aes(x = actual, y = predicted, color = grna)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("gRNA 1732" = "#0279EE", "gRNA 4894" = "#FF9400",
                                   "No gRNA" = "gray60")) +
    labs(
      x = "Actual log1p(edit count)",
      y = "Predicted log1p(edit count)",
      color = "Guide RNA",
      title = paste("Predicted vs Actual —", model_name)
    ) +
    theme_minimal(base_family = "sans") +
    theme(legend.position = "bottom")

  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 6, height = 5, dpi = 150)
    # Also save SVG
    svg_path <- sub("\\.png$", ".svg", out_path)
    ggsave(svg_path, p, width = 6, height = 5)
  }

  invisible(p)
}

# ---------------------------------------------------------------------------
# 2. Per-guide performance comparison
# ---------------------------------------------------------------------------

#' Plot per-guide Spearman correlation across models.
#'
#' @param summary_df Data frame from run_baseline_comparison$summary.
#' @param out_path String, path to save figure.
#' @return ggplot object (invisibly).
plot_per_guide_performance <- function(summary_df, out_path = NULL) {
  # Reshape for plotting
  plot_df <- data.frame(
    model = rep(summary_df$model, 3),
    metric = rep(c("Spearman (all)", "Spearman (1732)", "Spearman (4894)"),
                 each = nrow(summary_df)),
    value = c(summary_df$spearman_all,
              summary_df$spearman_1732,
              summary_df$spearman_4894)
  )
  plot_df$metric <- factor(plot_df$metric,
                            levels = c("Spearman (all)", "Spearman (1732)", "Spearman (4894)"))
  model_levels <- c("sequence_only", "structural_only", "chromatin_only", "full", "full_chromatin", "full_chromatin_intx")
  model_levels <- model_levels[model_levels %in% unique(plot_df$model)]
  plot_df$model <- factor(plot_df$model, levels = model_levels)

  model_colors <- c("sequence_only" = "#75A025", "structural_only" = "#FD9BED",
                     "chromatin_only" = "#FF9400", "full" = "#0279EE",
                     "full_chromatin" = "#E9ED4C", "full_chromatin_intx" = "#000000")
  model_colors <- model_colors[model_levels]

  p <- ggplot(plot_df, aes(x = model, y = value, fill = model)) +
    geom_col(width = 0.6) +
    facet_wrap(~ metric) +
    scale_fill_manual(values = model_colors) +
    labs(
      x = "Model",
      y = "Spearman correlation",
      title = "Per-guide performance comparison"
    ) +
    theme_minimal(base_family = "sans") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 8, height = 4, dpi = 150)
    svg_path <- sub("\\.png$", ".svg", out_path)
    ggsave(svg_path, p, width = 8, height = 4)
  }

  invisible(p)
}

# ---------------------------------------------------------------------------
# 3. SHAP feature importance by group
# ---------------------------------------------------------------------------

#' Compute and plot SHAP feature importance, grouped by feature category.
#'
#' @param model XGBoost model.
#' @param long_df Long-format data frame.
#' @param feature_cols Character vector, feature columns used.
#' @param groups Feature group definitions.
#' @param out_path String, path to save figure.
#' @return ggplot object (invisibly).
plot_shap_by_group <- function(model, long_df, feature_cols, groups, out_path = NULL) {
  # Sample for SHAP (use test set or a sample for speed)
  X <- as.matrix(long_df[, feature_cols, drop = FALSE])
  X[is.na(X)] <- 0

  # Compute SHAP values
  shap <- predict(model, X, predcontrib = TRUE)

  # Mean absolute SHAP per feature
  shap_mean <- colMeans(abs(shap[, -ncol(shap)]))  # exclude BIAS column
  shap_df <- data.frame(feature = names(shap_mean), shap_value = shap_mean)

  # Assign each feature to a group
  assign_group <- function(feat, groups) {
    for (gname in names(groups)) {
      if (feat %in% groups[[gname]]) return(gname)
    }
    "other"
  }
  shap_df$group <- sapply(shap_df$feature, assign_group, groups = groups)

  # Summarize by group
  group_summary <- aggregate(shap_value ~ group, data = shap_df, FUN = sum)
  group_summary <- group_summary[order(-group_summary$shap_value), ]

  # Plot
  p <- ggplot(group_summary, aes(x = reorder(group, shap_value), y = shap_value)) +
    geom_col(fill = "#0279EE", width = 0.6) +
    coord_flip() +
    labs(
      x = "Feature group",
      y = "Sum of mean |SHAP|",
      title = "Feature importance by group (SHAP)"
    ) +
    theme_minimal(base_family = "sans")

  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 6, height = 4, dpi = 150)
    svg_path <- sub("\\.png$", ".svg", out_path)
    ggsave(svg_path, p, width = 6, height = 4)
  }

  invisible(p)
}

# ---------------------------------------------------------------------------
# 4. Feature group ablation bar chart
# ---------------------------------------------------------------------------

#' Plot ablation comparison across the three baseline models.
#'
#' @param summary_df Data frame from run_baseline_comparison$summary.
#' @param out_path String, path to save figure.
#' @return ggplot object (invisibly).
plot_ablation <- function(summary_df, out_path = NULL) {
  model_levels <- c("sequence_only", "structural_only", "chromatin_only", "full", "full_chromatin", "full_chromatin_intx")
  model_levels <- model_levels[model_levels %in% unique(summary_df$model)]
  plot_df <- data.frame(
    model = factor(summary_df$model, levels = model_levels),
    spearman = summary_df$spearman_all,
    r2 = summary_df$r2_all
  )

  # Reshape
  plot_long <- data.frame(
    model = rep(plot_df$model, 2),
    metric = rep(c("Spearman", "R²"), each = nrow(plot_df)),
    value = c(plot_df$spearman, plot_df$r2)
  )

  model_colors <- c("sequence_only" = "#75A025", "structural_only" = "#FD9BED",
                     "chromatin_only" = "#FF9400", "full" = "#0279EE",
                     "full_chromatin" = "#E9ED4C", "full_chromatin_intx" = "#000000")
  model_colors <- model_colors[model_levels]

  p <- ggplot(plot_long, aes(x = model, y = value, fill = model)) +
    geom_col(width = 0.6) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = model_colors) +
    labs(
      x = "Model configuration",
      y = "Metric value",
      title = "Feature group ablation"
    ) +
    theme_minimal(base_family = "sans") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  if (!is.null(out_path)) {
    ggsave(out_path, p, width = 7, height = 4, dpi = 150)
    svg_path <- sub("\\.png$", ".svg", out_path)
    ggsave(svg_path, p, width = 7, height = 4)
  }

  invisible(p)
}

# ---------------------------------------------------------------------------
# 5. Mapper graph visualization
# ---------------------------------------------------------------------------

#' Plot Mapper graph with vertices colored by mean edit count.
#'
#' @param mapper_graph igraph object from compute_tda_features.
#' @param mapper_obj Mapper object (from TDAmapper).
#' @param edits Numeric vector, edit counts.
#' @param title String.
#' @param out_path String, path to save figure.
#' @return igraph plot (invisibly).
plot_mapper_graph <- function(mapper_graph, mapper_obj, edits, title = "Mapper graph",
                               out_path = NULL) {
  if (!is.null(out_path)) {
    png(out_path, width = 800, height = 600)
  }

  mean_edits <- sapply(mapper_obj$points_in_vertex, function(idx) {
    if (length(idx) == 0) return(0)
    mean(edits[idx], na.rm = TRUE)
  })

  palette <- colorRampPalette(c("#0279EE", "#E9ED4C", "#FF9400"))(100)
  min_ed <- min(mean_edits, na.rm = TRUE)
  max_ed <- max(mean_edits, na.rm = TRUE)
  if (max_ed > min_ed) {
    scale_ed <- round((mean_edits - min_ed) / (max_ed - min_ed) * 99) + 1
  } else {
    scale_ed <- rep(50, length(mean_edits))
  }

  plot(
    mapper_graph,
    layout = layout_with_fr(mapper_graph),
    vertex.color = palette[scale_ed],
    vertex.size = log(sapply(mapper_obj$points_in_vertex, length) + 1) * 2.5,
    vertex.label = NA,
    edge.color = "gray60",
    main = title
  )

  if (!is.null(out_path)) {
    dev.off()
    # Also save SVG
    svg_path <- sub("\\.png$", ".svg", out_path)
    svg(svg_path, width = 8, height = 6)
    plot(
      mapper_graph,
      layout = layout_with_fr(mapper_graph),
      vertex.color = palette[scale_ed],
      vertex.size = log(sapply(mapper_obj$points_in_vertex, length) + 1) * 2.5,
      vertex.label = NA,
      edge.color = "gray60",
      main = title
    )
    dev.off()
  }

  invisible(mapper_graph)
}
