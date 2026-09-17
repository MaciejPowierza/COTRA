# 03b_preprocessing.R — optional early preprocessing policies

preprocess_feature_matrix <- function(df_feat,
                                      preproc_mode = c("none","basic","aggressive"),
                                      corr_cutoff = 0.95,
                                      feature_set = NULL,
                                      verbose = FALSE) {
  preproc_mode <- match.arg(preproc_mode)

  if (is.null(feature_set)) feature_set <- get_feature_set(df_feat)

  X <- as.matrix(df_feat[, feature_set, drop = FALSE])
  orig_feature_names <- colnames(X)

  summary_row <- data.frame(
    preproc_mode = preproc_mode,
    n_features_input = ncol(X),
    n_removed_zero_var = 0L,
    n_removed_nzv = 0L,
    n_removed_linear_combo = 0L,
    n_removed_high_corr = 0L,
    n_features_output = ncol(X),
    stringsAsFactors = FALSE
  )

  if (preproc_mode == "none") {
    return(list(
      df_feat = df_feat,
      X = X,
      kept_features = colnames(X),
      removed_features = character(0),
      summary_row = summary_row
    ))
  }

  # -------------------------
  # basic: zero-var + nzv + exact linear combos + scale later in DR methods
  # aggressive: basic + high correlation pruning
  # -------------------------

  keep_names <- colnames(X)

  # zero variance
  sds <- apply(X, 2, sd, na.rm = TRUE)
  zero_var <- names(sds)[!is.finite(sds) | sds == 0]
  if (length(zero_var) > 0) {
    keep_names <- setdiff(keep_names, zero_var)
  }
  summary_row$n_removed_zero_var <- length(zero_var)

  X1 <- X[, keep_names, drop = FALSE]

  # near-zero variance
  nzv_idx <- caret::nearZeroVar(X1)
  nzv_names <- if (length(nzv_idx) > 0) colnames(X1)[nzv_idx] else character(0)
  if (length(nzv_names) > 0) {
    keep_names <- setdiff(keep_names, nzv_names)
  }
  summary_row$n_removed_nzv <- length(nzv_names)

  X2 <- X[, keep_names, drop = FALSE]

  # exact linear combinations
  lc <- tryCatch(caret::findLinearCombos(X2), error = function(e) NULL)
  lc_names <- character(0)
  if (!is.null(lc) && !is.null(lc$remove) && length(lc$remove) > 0) {
    lc_names <- colnames(X2)[lc$remove]
    keep_names <- setdiff(keep_names, lc_names)
  }
  summary_row$n_removed_linear_combo <- length(lc_names)

  X3 <- X[, keep_names, drop = FALSE]

  # aggressive only: remove highly correlated features
  hc_names <- character(0)
  if (preproc_mode == "aggressive" && ncol(X3) > 1) {
    cor_mat <- suppressWarnings(cor(X3, use = "pairwise.complete.obs"))
    cor_mat[!is.finite(cor_mat)] <- 0
    hc_idx <- caret::findCorrelation(cor_mat, cutoff = corr_cutoff)
    if (length(hc_idx) > 0) {
      hc_names <- colnames(X3)[hc_idx]
      keep_names <- setdiff(keep_names, hc_names)
    }
  }
  summary_row$n_removed_high_corr <- length(hc_names)

  X_out <- X[, keep_names, drop = FALSE]

  removed_features <- setdiff(orig_feature_names, keep_names)

  df_out <- df_feat
  # replace only the feature block, keep metadata columns intact
  df_out <- cbind(
    df_out[, setdiff(colnames(df_out), feature_set), drop = FALSE],
    as.data.frame(X_out, check.names = FALSE)
  )

  summary_row$n_features_output <- ncol(X_out)

  if (verbose) {
    print(summary_row)
  }

  list(
    df_feat = df_out,
    X = X_out,
    kept_features = keep_names,
    removed_features = removed_features,
    summary_row = summary_row
  )
}

compute_feature_diagnostics <- function(df_feat, feature_set = NULL) {
  if (is.null(feature_set)) feature_set <- get_feature_set(df_feat)

  X <- as.matrix(df_feat[, feature_set, drop = FALSE])

  # zero-var removed only for safe PCA diagnostics
  sds <- apply(X, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  X <- X[, keep, drop = FALSE]

  if (ncol(X) == 0) {
    return(data.frame(
      n_features = 0,
      rank_qr = NA_real_,
      effective_dim_participation_ratio = NA_real_,
      pcs_for_90_var = NA_real_,
      pcs_for_95_var = NA_real_,
      kaiser_count = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  X_scaled <- scale(X, center = TRUE, scale = TRUE)
  X_scaled[!is.finite(X_scaled)] <- 0

  rank_est <- qr(X_scaled)$rank

  pr <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
  lambda <- pr$sdev^2
  var_explained <- lambda / sum(lambda)
  cumvar <- cumsum(var_explained)

  k_90 <- which(cumvar >= 0.90)[1]
  k_95 <- which(cumvar >= 0.95)[1]
  k_kaiser <- sum(lambda > mean(lambda))
  eff_dim <- (sum(lambda)^2) / sum(lambda^2)

  data.frame(
    n_features = ncol(X_scaled),
    rank_qr = rank_est,
    effective_dim_participation_ratio = eff_dim,
    pcs_for_90_var = k_90,
    pcs_for_95_var = k_95,
    kaiser_count = k_kaiser,
    stringsAsFactors = FALSE
  )
}
