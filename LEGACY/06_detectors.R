
# 06_detectors.R — anomaly / one-class detectors

fit_mahalanobis_detector <- function(X_pca, normal_idx = NULL, alpha = 0.001) {
  X <- as.matrix(X_pca)
  if (is.null(normal_idx)) normal_idx <- seq_len(nrow(X))
  X_norm <- X[normal_idx, , drop = FALSE]
  mu <- colMeans(X_norm)
  Sigma <- cov(X_norm)
  Sigma_inv <- solve(Sigma)
  cutoff <- qchisq(1 - alpha, df = ncol(X))
  list(mu = mu, Sigma_inv = Sigma_inv, cutoff = cutoff, alpha = alpha)
}

score_mahalanobis <- function(detector, X_pca) {
  X <- as.matrix(X_pca)
  d2 <- mahalanobis(X, center = detector$mu, cov = solve(detector$Sigma_inv), inverted = TRUE)
  data.frame(anomaly_score = d2, is_anomaly = d2 > detector$cutoff)
}

fit_ocsvm <- function(X_pca, normal_idx = NULL, nu = 0.05, kernel = "radial", gamma = "auto") {
  X <- as.matrix(X_pca)
  if (is.null(normal_idx)) normal_idx <- seq_len(nrow(X))
  X_norm <- X[normal_idx, , drop = FALSE]
  if (identical(gamma, "auto")) gamma <- 1 / ncol(X_norm)

  e1071::svm(
    x = X_norm, y = NULL,
    type = "one-classification",
    kernel = kernel, nu = nu, gamma = gamma,
    scale = TRUE
  )
}

score_ocsvm <- function(model, X_pca) {
  X <- as.matrix(X_pca)
  dv <- attr(predict(model, X, decision.values = TRUE), "decision.values")
  # Larger dv = more normal; invert so larger = more anomalous
  scores <- -as.numeric(dv)
  data.frame(anomaly_score = scores, is_anomaly = scores > 0)
}

fit_iso_forest <- function(X_pca, ntrees = 500, sample_size = 256) {
  isotree::isolation.forest(as.matrix(X_pca), ntrees = ntrees, sample_size = sample_size)
}

score_iso_forest <- function(model, X_pca, quantile_cut = 0.95) {
  scores <- as.numeric(predict(model, as.matrix(X_pca)))
  thr <- as.numeric(stats::quantile(scores, probs = quantile_cut, na.rm = TRUE))
  data.frame(anomaly_score = scores, is_anomaly = scores >= thr, threshold = thr)
}
