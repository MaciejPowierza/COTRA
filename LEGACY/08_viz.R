
# 08_viz.R — plotting helpers

plot_outlier_scatter <- function(pc_scores, title = "PC scatter") {
  ggplot2::ggplot(pc_scores, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(ggplot2::aes(color = is_anomaly), shape = 16, alpha = 0.6) +
    ggplot2::geom_point(
      data = subset(pc_scores, is_on_target),
      ggplot2::aes(x = PC1, y = PC2),
      shape = 17, color = "black", size = 3
    ) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "salmon", `TRUE` = "steelblue"),
                                name = "Classified as outlier") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "PC1", y = "PC2")
}
