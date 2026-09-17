# --- utilities (no extra packages) ---
auc_roc <- function(score, y01) {
  o <- order(score, decreasing = TRUE)
  y <- y01[o]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score[o], ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auc_pr <- function(score, y01) {
  o <- order(score, decreasing = TRUE)
  y <- y01[o]
  P <- sum(y == 1)
  if (P == 0) return(NA_real_)
  tp <- cumsum(y == 1)
  fp <- cumsum(y == 0)
  recall <- tp / P
  precision <- tp / pmax(1, tp + fp)
  recall <- c(0, recall)
  precision <- c(1, precision)
  sum((recall[-1] - recall[-length(recall)]) * (precision[-1] + precision[-length(precision)]) / 2)
}

confusion_from_flag <- function(flag, y_pos) {
  tp <- sum(flag & y_pos)
  fp <- sum(flag & !y_pos)
  tn <- sum(!flag & !y_pos)
  fn <- sum(!flag & y_pos)
  list(tp=tp, fp=fp, tn=tn, fn=fn)
}

run_one_pass <- function(df_feat,
                         dataset_i = NA_integer_,
                         selection_method = c("none", "topvar", "lscore", "fosmod", "mcfs", "block_reps"),
                         dim_method = c("none","weighted_pca","plain_pca","hpca", "nmf", "rgcca", "sgcca", "bpca","fa", "lpe", "ica"),
                         detector = c("ocsvm","isoforest","mahalanobis"),
			 hpca_feature_mode = c("super_blocks"),
                         normal_idx = NULL,
                         weight_col = NULL,
                         k_select = NULL,
			 k_dr,
                         iso_quantile = 0.95,
                         maha_alpha = 0.001,
                         ocsvm_nu = 0.05,
		         nmf_nonneg_mode = c("shift","drop_negative_features"),
                         nmf_feature_subset = c("all","nonnegative_only"),
                         nmf_seed = 123,
                         preproc_mode = c("none", "basic", "aggressive"),
                         corr_cutoff = 0.95,
                         edits_col = "edits",
                         ontarget_policy = c("first","all"),
                         positive_class = c("offtarget","ontarget"),
                         fs_preprocess = "center",
                         fs_graph_type = c("knn", 10),
                         fs_lscore_t = 1,
                         fs_mcfs_K = NULL,
                         fs_mcfs_lambda = 1,
                         fs_mcfs_t = 10,
                         block_rep_per_block = 1,
                         rgcca_scheme = "factorial",
                         rgcca_tau = "optimal",
                         rgcca_sparsity = 0.8,
                         rgcca_connection_mode = "full",
                         rgcca_scores_use = c("first","first2","all"),
                         rgcca_superblock = FALSE,
                         lpe_preprocess = "center",
                         lpe_numk = NULL,
                         ica_alg_typ = "parallel",
                         ica_fun = "logcosh",
                         ica_alpha = 1.0,
                         ica_method = "C",
                         ica_row_norm = FALSE,
                         ica_maxit = 200,
                         ica_tol = 1e-4) {

  print(dim_method)
  selection_method <- match.arg(selection_method)
  dim_method <- match.arg(dim_method)
  detector <- match.arg(detector)
  hpca_feature_mode <- match.arg(hpca_feature_mode)
  nmf_nonneg_mode <- match.arg(nmf_nonneg_mode)
  nmf_feature_subset <- match.arg(nmf_feature_subset)
  ontarget_policy <- match.arg(ontarget_policy)
  positive_class <- match.arg(positive_class)
  preproc_mode <- match.arg(preproc_mode)

  if (dim_method %in% c("rgcca", "sgcca")) {
    rgcca_scores_use <- match.arg(rgcca_scores_use, c("first","first2","all"))
  }

  #1) Feature selector

  sel <- select_features(
    df_feat = df_feat,
    method = selection_method,
    k_select = k_select,
    fs_preprocess = fs_preprocess,
    fs_graph_type = fs_graph_type,
    fs_lscore_t = fs_lscore_t,
    fs_mcfs_K = fs_mcfs_K,
    fs_mcfs_lambda = fs_mcfs_lambda,
    fs_mcfs_t = fs_mcfs_t,
    block_rep_per_block = block_rep_per_block
  )

  # 1b) Dim reduction -> 2 PCs
  #dr <- dim_reduce(df_feat, method = dim_method, weight_col = weight_col, k=k, hpca_feature_mode = hpca_feature_mode)
  #dr <- dim_reduce(
  #    df_feat,
  #    method = dim_method,
  #    weight_col = weight_col,
  #    k = k,
  #    hpca_feature_mode = hpca_feature_mode,
  #    nmf_nonneg_mode = nmf_nonneg_mode,
  #    nmf_feature_subset = nmf_feature_subset,
  #    nmf_seed = nmf_seed,
  #    preproc_mode = preproc_mode,
  #    corr_cutoff = corr_cutoff,
  #    rgcca_scheme = rgcca_scheme,
  #    rgcca_tau = rgcca_tau,
  #    rgcca_sparsity = rgcca_sparsity,
  #    rgcca_connection_mode = rgcca_connection_mode,
  #    rgcca_scores_use = rgcca_scores_use,
  #    rgcca_superblock = rgcca_superblock,
  #    lpe_preprocess = lpe_preprocess,
  #    lpe_numk = lpe_numk,
  #    ica_alg_typ = ica_alg_typ,
  #    ica_fun = ica_fun,
  #    ica_alpha = ica_alpha,
  #    ica_method = ica_method,
  #    ica_row_norm = ica_row_norm,
  #    ica_maxit = ica_maxit,
  #    ica_tol = ica_tol)

  dr <- dim_reduce(
    df_feat,
    method = dim_method,
    weight_col = weight_col,
    k_dr = k_dr,
    X_input = sel$X,
    feature_names_input = colnames(sel$X),
    hpca_feature_mode = hpca_feature_mode,
    nmf_nonneg_mode = nmf_nonneg_mode,
    nmf_feature_subset = nmf_feature_subset,
    nmf_seed = nmf_seed,
    preproc_mode = preproc_mode,
    corr_cutoff = corr_cutoff,
    rgcca_scheme = rgcca_scheme,
    rgcca_tau = rgcca_tau,
    rgcca_sparsity = rgcca_sparsity,
    rgcca_connection_mode = rgcca_connection_mode,
    rgcca_scores_use = rgcca_scores_use,
    rgcca_superblock = rgcca_superblock,
    lpe_preprocess = lpe_preprocess,
    lpe_numk = lpe_numk,
    ica_alg_typ = ica_alg_typ,
    ica_fun = ica_fun,
    ica_alpha = ica_alpha,
    ica_method = ica_method,
    ica_row_norm = ica_row_norm,
    ica_maxit = ica_maxit,
    ica_tol = ica_tol
  )

  pcs <- dr$pcs
  print(head(pcs))
  if(dim_method == "hpca"){
  #print("whoa!")
  k <- dr$k_eff}
  k_eff_used <- if (!is.null(dr$k_eff)) dr$k_eff else ncol(pcs)
  print("k_eff_used")
  print(k_eff_used)
  #print(k)
  #print(head(pcs))

  prep_n_in <- NA_real_
  prep_n_out <- NA_real_
  prep_zero_var <- NA_real_
  prep_nzv <- NA_real_
  prep_lincombo <- NA_real_
  prep_highcorr <- NA_real_
  prep_effdim <- NA_real_
  prep_rank <- NA_real_

  if (!is.null(dr$preprocessing) && !is.null(dr$preprocessing$summary_row)) {
    prep_n_in <- dr$preprocessing$summary_row$n_features_input
    prep_n_out <- dr$preprocessing$summary_row$n_features_output
    prep_zero_var <- dr$preprocessing$summary_row$n_removed_zero_var
    prep_nzv <- dr$preprocessing$summary_row$n_removed_nzv
    prep_lincombo <- dr$preprocessing$summary_row$n_removed_linear_combo
    prep_highcorr <- dr$preprocessing$summary_row$n_removed_high_corr
  }

  if (!is.null(dr$feature_diagnostics) && nrow(dr$feature_diagnostics) > 0) {
    prep_effdim <- dr$feature_diagnostics$effective_dim_participation_ratio
    prep_rank <- dr$feature_diagnostics$rank_qr
  }

  # 2) Detector
  if (detector == "mahalanobis") {
    det_obj <- fit_mahalanobis_detector(pcs, normal_idx = normal_idx, alpha = maha_alpha)
    det_scores <- score_mahalanobis(det_obj, pcs)

    # signed distance consistent with your original script:
    # positive => beyond cutoff => anomaly
    det_scores$signed_distance <- det_scores$anomaly_score - det_obj$cutoff
    det_scores$is_anomaly <- det_scores$signed_distance > 0

  } else if (detector == "ocsvm") {
    det_obj <- fit_ocsvm(pcs, normal_idx = normal_idx, nu = ocsvm_nu)
    det_scores <- score_ocsvm(det_obj, pcs)

    # In your current implementation anomaly_score = -decision_value, threshold 0.
    # This matches your original score_ocsvm_distance() signed_distance exactly.
    det_scores$signed_distance <- det_scores$anomaly_score
    det_scores$is_anomaly <- det_scores$signed_distance > 0

  } else {
    det_obj <- fit_iso_forest(pcs)
    det_scores <- score_iso_forest(det_obj, pcs, quantile_cut = iso_quantile)

    # signed distance consistent with your original script:
    # signed_distance = score - threshold
    if (!("threshold" %in% names(det_scores))) {
      stop("score_iso_forest() must return a 'threshold' column to construct signed_distance.")
    }
    det_scores$signed_distance <- det_scores$anomaly_score - det_scores$threshold
    det_scores$is_anomaly <- det_scores$signed_distance > 0
  }

  # Use signed_distance everywhere downstream (AUC, on-target score, etc.)
  score_vec <- det_scores$signed_distance


  # 3) On-target definition
  is_on_target <- if (edits_col %in% names(df_feat)) df_feat[[edits_col]] == 0 else rep(FALSE, nrow(df_feat))
  ont_idx <- which(is_on_target)

  # 4) On-target score(s) (direct analogue of your old results_df columns)
  ont_scores <- NA
  if (length(ont_idx) > 0) {
    if (ontarget_policy == "first") ont_scores <- score_vec[ont_idx[1]]
    else ont_scores <- score_vec[ont_idx]
  }

  # 5) Define positives for AUC/PR-AUC
  y_pos <- if (positive_class == "offtarget") !is_on_target else is_on_target
  y01 <- as.integer(y_pos)

  roc_auc <- auc_roc(score_vec, y01)
  pr_auc  <- auc_pr(score_vec, y01)

  # 6) Threshold metrics if available
  tpr <- fpr <- precision <- recall <- NA_real_
  cm <- NULL
  if ("is_anomaly" %in% names(det_scores)) {
    flag <- as.logical(det_scores$is_anomaly)
    cm <- confusion_from_flag(flag, y_pos)
    tpr <- if ((cm$tp + cm$fn) > 0) cm$tp / (cm$tp + cm$fn) else NA_real_
    fpr <- if ((cm$fp + cm$tn) > 0) cm$fp / (cm$fp + cm$tn) else NA_real_
    precision <- if ((cm$tp + cm$fp) > 0) cm$tp / (cm$tp + cm$fp) else NA_real_
    recall <- tpr
  }

  # 7) Per-row outputs
  #pc_scores <- cbind(
  #  data.frame(PC1 = pcs[,1], PC2 = pcs[,2]),
  #  det_scores,
  #  is_on_target = is_on_target
  #)

  pc_df <- data.frame(
    PC1 = pcs[, 1],
    PC2 = if (ncol(pcs) >= 2) pcs[, 2] else rep(NA_real_, nrow(pcs))
  )

  pc_scores <- cbind(
    pc_df,
    det_scores,
    is_on_target = is_on_target
  )
  
  if_anomaly <-pc_scores[pc_scores$is_on_target, c("is_anomaly")]
  #print(if_anomaly)
  #print(class(if_anomaly))
  
  ###===Beginning of the insertion, 09.02===
  ranks <- rank(-score_vec, ties.method = "average")
  #print("SCORE VEC")
  #print(score_vec)
  #print(nrow(df_feat))
  #print(length(score_vec))
  #print("RANKS")
  #print(ranks)
  
  rank_on_target <- NA_real_
  rr <- NA_real_
  rank_frac <- NA_real_
  hit_at_1 <- NA_integer_
  hit_at_5 <- NA_integer_
  hit_at_10 <- NA_integer_
  
  if (length(ont_idx) > 0) {
    # If there's exactly one on-target, this is just one number.
    # If there are multiple on-targets (rare), use the best (minimum) rank as "how fast you find any on-target".
    ont_ranks <- ranks[ont_idx]

    if (ontarget_policy == "first") {
      rank_on_target <- ont_ranks[1]
    } else {
      rank_on_target <- min(ont_ranks, na.rm = TRUE)
    }

    rank_frac <- rank_on_target / length(score_vec)
    rr <- 1 / rank_on_target
    hit_at_1  <- as.integer(rank_on_target <= 1)
    hit_at_5  <- as.integer(rank_on_target <= 5)
    hit_at_10 <- as.integer(rank_on_target <= 10)
  }
  ###===End of the insertion, 09.02===

  ###===Beginning of the insertion, 08.04===

  rgcca_mean_ncomp <- NA_real_
  rgcca_median_ncomp <- NA_real_
  rgcca_min_ncomp <- NA_real_
  rgcca_max_ncomp <- NA_real_
  rgcca_sum_ncomp <- NA_real_
  rgcca_n_blocks <- NA_real_

  if (dim_method %in% c("rgcca", "sgcca")) {
    if (!is.null(dr$mean_ncomp_rgcca))   rgcca_mean_ncomp <- dr$mean_ncomp_rgcca
    if (!is.null(dr$median_ncomp_rgcca)) rgcca_median_ncomp <- dr$median_ncomp_rgcca
    if (!is.null(dr$min_ncomp_rgcca))    rgcca_min_ncomp <- dr$min_ncomp_rgcca
    if (!is.null(dr$max_ncomp_rgcca))    rgcca_max_ncomp <- dr$max_ncomp_rgcca
    if (!is.null(dr$sum_ncomp_rgcca))    rgcca_sum_ncomp <- dr$sum_ncomp_rgcca
    if (!is.null(dr$n_blocks_rgcca))     rgcca_n_blocks <- dr$n_blocks_rgcca
  }

  ###===End of the insertion, 08.04===

  ###===Beginning of the insertion, 10.04===

  rgcca_scores_mode_used <- NA_character_

  if (dim_method %in% c("rgcca", "sgcca") && !is.null(dr$rgcca_scores_use)) {
    rgcca_scores_mode_used <- dr$rgcca_scores_use
  }

  ###===End of the insertion, 10.04===

  # 8) One-row metrics (easy to bind across datasets / methods)
  metrics_row <- data.frame(
    dataset_i = dataset_i,
    dim_method = dim_method,
    detector = detector,
    positive_class = positive_class,
    k_pcs = k_dr,
    selection_method = selection_method, 
    k_select_req = if (is.null(k_select)) NA_real_ else k_select,
    k_select_eff = if (!is.null(sel$k_select_eff)) sel$k_select_eff else NA_real_,
    k_dr_req = k_dr,
    k_eff = k_eff_used,
    n = nrow(pcs),
    n_on_target = sum(is_on_target),
    on_target_score_signed_distance = if (ontarget_policy == "first") ont_scores else NA_real_,
    on_target_score_mean = if (ontarget_policy == "all" && length(ont_idx) > 0) mean(ont_scores) else NA_real_,
    on_target_score_sd   = if (ontarget_policy == "all" && length(ont_idx) > 1) sd(ont_scores) else NA_real_,
    roc_auc = roc_auc,
    pr_auc  = pr_auc,
    tpr = tpr,
    fpr = fpr,
    precision = precision,
    recall = recall,
    if_anomaly = if_anomaly,
    rank_on_target = rank_on_target,
    rank_frac = rank_frac, 
    rr = rr,
    hit_at_1 = hit_at_1,
    hit_at_5 = hit_at_5,
    hit_at_10 = hit_at_10,
    preproc_mode = preproc_mode,
    prep_n_features_input = prep_n_in,
    prep_n_features_output = prep_n_out,
    prep_removed_zero_var = prep_zero_var,
    prep_removed_nzv = prep_nzv,
    prep_removed_linear_combo = prep_lincombo,
    prep_removed_high_corr = prep_highcorr,
    prep_rank_qr = prep_rank,
    prep_effective_dim = prep_effdim,
    rgcca_scores_use = rgcca_scores_mode_used,
    rgcca_n_blocks = rgcca_n_blocks,
    rgcca_mean_ncomp = rgcca_mean_ncomp,
    rgcca_median_ncomp = rgcca_median_ncomp,
    rgcca_min_ncomp = rgcca_min_ncomp,
    rgcca_max_ncomp = rgcca_max_ncomp,
    rgcca_sum_ncomp = rgcca_sum_ncomp,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )

  list(
    pcs = pcs,
    dimred = dr,
    detector_model = det_obj,
    pc_scores = pc_scores,
    metrics_row = metrics_row,
    metrics = list(
      ontarget_idx = ont_idx,
      ontarget_scores = ont_scores,
      roc_auc = roc_auc,
      pr_auc = pr_auc,
      confusion = cm
    )
  )
}

#===Beginning of the insertion, 13.04===

make_failed_metrics_row <- function(df_feat, dim_method, selection_method, detector, k_dr, k_select, preproc_mode,
                                    positive_class = NA_character_, error_message = NA_character_) {
  data.frame(
    dim_method = dim_method,
    detector = detector,
    positive_class = positive_class,
    k_pcs = k,
    selection_method = selection_method,
    k_select_req = if (is.null(k_select)) NA_real_ else k_select,
    k_select_eff = NA_real_,
    k_dr_req = k_dr,
    k_eff = NA_real_,
    n = nrow(df_feat),
    n_on_target = if ("edits" %in% names(df_feat)) sum(df_feat$edits == 0) else NA_real_,
    on_target_score_signed_distance = NA_real_,
    on_target_score_mean = NA_real_,
    on_target_score_sd = NA_real_,
    roc_auc = NA_real_,
    pr_auc = NA_real_,
    tpr = NA_real_,
    fpr = NA_real_,
    precision = NA_real_,
    recall = NA_real_,
    if_anomaly = NA,
    rank_on_target = NA_real_,
    rank_frac = NA_real_,
    rr = NA_real_,
    hit_at_1 = NA_integer_,
    hit_at_5 = NA_integer_,
    hit_at_10 = NA_integer_,
    preproc_mode = preproc_mode,
    prep_n_features_input = NA_real_,
    prep_n_features_output = NA_real_,
    prep_removed_zero_var = NA_real_,
    prep_removed_nzv = NA_real_,
    prep_removed_linear_combo = NA_real_,
    prep_removed_high_corr = NA_real_,
    prep_rank_qr = NA_real_,
    prep_effective_dim = NA_real_,
    rgcca_scores_use = NA_character_,
    rgcca_n_blocks = NA_real_,
    rgcca_mean_ncomp = NA_real_,
    rgcca_median_ncomp = NA_real_,
    rgcca_min_ncomp = NA_real_,
    rgcca_max_ncomp = NA_real_,
    rgcca_sum_ncomp = NA_real_,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

#===End of the insertion, 13.04===

#===Beginning of the insertion, 17.02===

run_experiment_metrics <- function(datasets,
                                   selection_methods = c("none", "topvar", "lscore", "fosmod", "mcfs", "block_reps"),
                                   k_select_grid = c(10, 20, 50, 100),
                                   dim_methods = c("none","weighted_pca","plain_pca","hpca","nmf","rgcca","sgcca","bpca","fa","lpe","ica"),
                                   detectors = c("ocsvm","isoforest"),
                                   preproc_modes = c("none","basic","aggressive"),
                                   hpca_feature_mode = "super_blocks",
                                   rgcca_scores_uses = c("first","first2","all"),
                                   k_dr_grid = c(2:16),
                                   seed = 123,
                                   shuffle = TRUE,
                                   ...) {
  set.seed(seed)
  ord <- seq_along(datasets)
  if (shuffle) ord <- sample(ord)

  rows <- list()
  kk <- 1

  make_failed_row <- function(i, pm, sm, ks, dm, det, kdr, msg) {
    data.frame(
      dataset_i = i,
      dim_method = dm,
      detector = det,
      positive_class = NA_character_,
      k_pcs = kdr,
      selection_method = sm,
      k_select_req = if (is.null(ks)) NA_real_ else ks,
      k_select_eff = NA_real_,
      k_dr_req = kdr,
      k_eff = NA_real_,
      n = nrow(datasets[[i]]),
      n_on_target = if ("edits" %in% names(datasets[[i]])) sum(datasets[[i]]$edits == 0) else NA_real_,
      on_target_score_signed_distance = NA_real_,
      on_target_score_mean = NA_real_,
      on_target_score_sd = NA_real_,
      roc_auc = NA_real_,
      pr_auc = NA_real_,
      tpr = NA_real_,
      fpr = NA_real_,
      precision = NA_real_,
      recall = NA_real_,
      if_anomaly = NA,
      rank_on_target = NA_real_,
      rank_frac = NA_real_,
      rr = NA_real_,
      hit_at_1 = NA_integer_,
      hit_at_5 = NA_integer_,
      hit_at_10 = NA_integer_,
      preproc_mode = pm,
      prep_n_features_input = NA_real_,
      prep_n_features_output = NA_real_,
      prep_removed_zero_var = NA_real_,
      prep_removed_nzv = NA_real_,
      prep_removed_linear_combo = NA_real_,
      prep_removed_high_corr = NA_real_,
      prep_rank_qr = NA_real_,
      prep_effective_dim = NA_real_,
      rgcca_scores_use = NA_character_,
      rgcca_n_blocks = NA_real_,
      rgcca_mean_ncomp = NA_real_,
      rgcca_median_ncomp = NA_real_,
      rgcca_min_ncomp = NA_real_,
      rgcca_max_ncomp = NA_real_,
      rgcca_sum_ncomp = NA_real_,
      error_message = msg,
      stringsAsFactors = FALSE
    )
  }

  for (i in ord) {
    for (pm in preproc_modes) {
      for (sm in selection_methods) {

        k_select_values <- if (sm == "none") list(NULL) else as.list(k_select_grid)

        for (ks in k_select_values) {
          for (dm in dim_methods) {

            score_modes <- if (dm %in% c("rgcca", "sgcca")) rgcca_scores_uses else "first"

            for (rsm in score_modes) {
              for (det in detectors) {
                for (kdr in k_dr_grid) {

                  out <- tryCatch(
                    run_one_pass(
                      datasets[[i]],
                      dataset_i = i,
                      selection_method = sm,
                      dim_method = dm,
                      detector = det,
                      k_select = ks,
                      k_dr = kdr,
                      preproc_mode = pm,
                      rgcca_scores_use = rsm,
                      ...
                    ),
                    error = function(e) e
                  )

                  if (inherits(out, "error")) {
                    rows[[kk]] <- make_failed_row(
                      i = i,
                      pm = pm,
                      sm = sm,
                      ks = ks,
                      dm = dm,
                      det = det,
                      kdr = kdr,
                      msg = conditionMessage(out)
                    )
                  } else {
                    out$metrics_row$dataset_i <- i
                    out$metrics_row$error_message <- NA_character_
                    rows[[kk]] <- out$metrics_row
                  }

                  kk <- kk + 1
                }
              }
            }
          }
        }
      }
    }
  }

  do.call(rbind, rows)
  #print("AAAAaaaaarrrgggh!")
}

#=============================================================EOF=============================================================================

run_experiment_metrics2 <- function(datasets,
                                   selection_methods = c("none", "topvar", "lscore", "fosmod", "mcfs", "block_reps"),
                                   k_select_grid = c(10, 20, 50, 100),
                                   dim_methods = c("none","weighted_pca","plain_pca","hpca","nmf","rgcca","sgcca","bpca","fa","lpe","ica"),
                                   detectors = c("ocsvm","isoforest"),
                                   preproc_modes = c("none","basic","aggressive"),
                                   hpca_feature_mode = "super_blocks",
                                   rgcca_scores_uses = c("first","first2","all"),
                                   k_dr_grid = c(2:16),
                                   seed = 123,
                                   shuffle = TRUE,
                                   ...) {
  set.seed(seed)
  ord <- seq_along(datasets)
  if (shuffle) ord <- sample(ord)

  rows <- list()
  kk <- 1
  for (i in ord) {
    for (pm in preproc_modes) {
    for (sm in selection_methods) {

      k_select_values <- if (sm == "none") NA else k_select_grid

      for (ks in k_select_values) {
      for (dm in dim_methods) {

        score_modes <- if (dm %in% c("rgcca", "sgcca")) rgcca_scores_uses else "first"

        for (rsm in score_modes) {
          for (det in detectors) {
            for (kdr in k_dr_grid) {
                out <- tryCatch(
                  run_one_pass(
                    datasets[[i]],
                    selection_method = sm,
                    dim_method = dm,
                    detector = det,
                    k_select = if (is.na(ks)) NULL else ks,
                    k_dr = kdr,
                    preproc_mode = pm,
                    rgcca_scores_use = rsm,
                    ...
                  ),
                  error = function(e) e
                )

                if (inherits(out, "error")) {
                  rows[[kk]] <- data.frame(
                    dataset_i = i,
                    preproc_mode = pm,
                    selection_method = sm,
                    k_select_req = if (is.na(ks)) NA_real_ else ks,
                    dim_method = dm,
                    detector = det,
                    k_dr_req = kdr,
                    k_eff = NA_real_,
                    error_message = conditionMessage(out),
                    stringsAsFactors = FALSE
                  )
                } else {
                  out$metrics_row$dataset_i <- i
                  out$metrics_row$error_message <- NA_character_
                  rows[[kk]] <- out$metrics_row
                }

                kk <- kk + 1
                }
              }
            }
          }
        }
      }
    }
  }
  do.call(rbind, rows)
}

#===End of the insertion, 17.02===

#run_experiment_metrics2 <- function(datasets,
#                                   dim_methods = c("weighted_pca","plain_pca","hpca"),
#                                   detectors = c("ocsvm","isoforest"),
#                                   seed = 123,
#                                   shuffle = TRUE,
#                                   ...) {
#  set.seed(seed)
#  ord <- seq_along(datasets)
#  if (shuffle) ord <- sample(ord)

#  rows <- list()
#  k <- 1
#  for (i in ord) {
#    for (dm in dim_methods) {
#      for (det in detectors) {
#        out <- run_one_pass(datasets[[i]], dim_method = dm, detector = det, ...)
#        out$metrics_row$dataset_i <- i
#        rows[[k]] <- out$metrics_row
#        k <- k + 1
#      }
#    }
#  }
#  do.call(rbind, rows)
#}
