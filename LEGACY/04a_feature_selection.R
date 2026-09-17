# 04a_feature_selection.R — unsupervised feature selection layer

sanitize_selected_matrix <- function(Xsel) {
  Xsel <- as.matrix(Xsel)
  storage.mode(Xsel) <- "double"
  Xsel[!is.finite(Xsel)] <- NA_real_

  for (j in seq_len(ncol(Xsel))) {
    col <- Xsel[, j]
    m <- mean(col, na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    col[is.na(col)] <- m
    Xsel[, j] <- col
  }

  sds <- apply(Xsel, 2, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 0
  Xsel <- Xsel[, keep, drop = FALSE]

  if (ncol(Xsel) == 0L) return(NULL)

  Xsel <- scale(Xsel, center = TRUE, scale = TRUE)
  Xsel[!is.finite(Xsel)] <- 0
  Xsel
}

select_features <- function(df_feat,
                            method = c("none", "topvar", "lscore", "fosmod", "mcfs", "block_reps"),
                            k_select = NULL,
                            fs_preprocess = "center",
                            fs_graph_type = c("knn", 10),
                            fs_lscore_t = 1,
                            fs_mcfs_K = NULL,
                            fs_mcfs_lambda = 1,
                            fs_mcfs_t = 10,
                            block_rep_per_block = 1) {
  method <- match.arg(method)

  feature_set <- get_feature_set(df_feat)
  X <- as.matrix(df_feat[, feature_set, drop = FALSE])

  # remove constant/invalid columns before selection
  sds <- apply(X, 2, sd, na.rm = TRUE)
  keep0 <- is.finite(sds) & sds > 0
  X2 <- X[, keep0, drop = FALSE]
  feat2 <- colnames(X2)

  if (ncol(X2) == 0L) stop("select_features: no valid features available.")

  if (is.null(k_select)) {
    k_select <- ncol(X2)
  }
  k_target <- min(as.integer(k_select), ncol(X2))

  if (method == "none") {
    Xsel <- sanitize_selected_matrix(X2)
    return(list(
      X = Xsel,
      selected_features = colnames(Xsel),
      model = NULL,
      selection_method = "none",
      k_select_req = k_select,
      k_select_eff = ncol(Xsel)
    ))
  }

  if (method == "topvar") {
    vars <- apply(X2, 2, var, na.rm = TRUE)
    ord <- order(vars, decreasing = TRUE)
    sel_idx <- ord[seq_len(k_target)]
    Xsel <- sanitize_selected_matrix(X2[, sel_idx, drop = FALSE])

   print("SELECTED VARS inside TOPVAR")
   print(head(Xsel))

    return(list(
      X = Xsel,
      selected_features = colnames(Xsel),
      model = NULL,
      selection_method = "topvar",
      feature_scores = vars[sel_idx],
      k_select_req = k_select,
      k_select_eff = ncol(Xsel)
    ))
  }

  if (method == "lscore") {
    fit <- Rdimtools::do.lscore(
      X2,
      ndim = k_target,
      preprocess = fs_preprocess,
      type = fs_graph_type,
      t = fs_lscore_t
    )
    sel_idx <- fit$featidx
    Xsel <- sanitize_selected_matrix(X2[, sel_idx, drop = FALSE])

    return(list(
      X = Xsel,
      selected_features = colnames(Xsel),
      model = fit,
      selection_method = "lscore",
      feature_scores = if (!is.null(fit$lscore)) fit$lscore[sel_idx] else NULL,
      k_select_req = k_select,
      k_select_eff = ncol(Xsel)
    ))
  }

  if (method == "fosmod") {
    fit <- Rdimtools::do.fosmod(
      X2,
      ndim = k_target,
      preprocess = fs_preprocess
    )
    sel_idx <- fit$featidx
    Xsel <- sanitize_selected_matrix(X2[, sel_idx, drop = FALSE])

    return(list(
      X = Xsel,
      selected_features = colnames(Xsel),
      model = fit,
      selection_method = "fosmod",
      k_select_req = k_select,
      k_select_eff = ncol(Xsel)
    ))
  }

  if (method == "mcfs") {
    K_used <- if (is.null(fs_mcfs_K)) max(2L, min(10L, round(sqrt(nrow(X2))))) else as.integer(fs_mcfs_K)

    fit <- Rdimtools::do.mcfs(
      X2,
      ndim = k_target,
      type = fs_graph_type,
      preprocess = fs_preprocess,
      K = K_used,
      lambda = fs_mcfs_lambda,
      t = fs_mcfs_t
    )
    sel_idx <- fit$featidx
    Xsel <- sanitize_selected_matrix(X2[, sel_idx, drop = FALSE])

    return(list(
      X = Xsel,
      selected_features = colnames(Xsel),
      model = fit,
      selection_method = "mcfs",
      mcfs_K = K_used,
      mcfs_lambda = fs_mcfs_lambda,
      mcfs_t = fs_mcfs_t,
      k_select_req = k_select,
      k_select_eff = ncol(Xsel)
    ))
  }

  if (method == "block_reps") {
    Xdf <- as.data.frame(X2)
    routing_df <- route_features(Xdf)
    b_sp <- build_blocks_spearman(Xdf, routing_df, minClusterSize = 6)
    b_j  <- build_blocks_jaccard(Xdf, routing_df)
    b_h  <- build_blocks_hellinger(Xdf, routing_df)

    blocks_all <- c(b_sp$blocks, b_j$blocks, b_h$blocks)
    blocks_all <- blocks_all[!vapply(blocks_all, is.null, logical(1))]
    blocks_all <- blocks_all[vapply(blocks_all, length, integer(1)) > 0]

    reps <- character(0)
    for (bn in names(blocks_all)) {
      feats <- intersect(blocks_all[[bn]], colnames(Xdf))
      if (length(feats) == 0L) next

      vars_b <- apply(Xdf[, feats, drop = FALSE], 2, var, na.rm = TRUE)
      ord_b <- order(vars_b, decreasing = TRUE)
      kk <- min(block_rep_per_block, length(ord_b))
      reps <- c(reps, feats[ord_b[seq_len(kk)]])
    }

    reps <- unique(reps)
    if (length(reps) == 0L) stop("block_reps: no valid representative features selected.")

    vars_rep <- apply(Xdf[, reps, drop = FALSE], 2, var, na.rm = TRUE)
    ord_rep <- order(vars_rep, decreasing = TRUE)
    reps <- reps[ord_rep[seq_len(min(k_target, length(ord_rep)))]]

    Xsel <- sanitize_selected_matrix(Xdf[, reps, drop = FALSE])

    return(list(
      X = Xsel,
      selected_features = colnames(Xsel),
      model = list(routing_df = routing_df, blocks_all = blocks_all),
      selection_method = "block_reps",
      block_rep_per_block = block_rep_per_block,
      k_select_req = k_select,
      k_select_eff = ncol(Xsel)
    ))
  }

  stop("Unknown selection method.")
}
