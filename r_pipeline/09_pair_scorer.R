# 09_pair_scorer.R — Pair scorer: data assembly, training, evaluation
#
# Assembles long-format training data from (locus, sample) pairs,
# concatenates all feature groups, trains XGBoost regressor on log1p(edit_count),
# and evaluates with Spearman, R², and Recall@k.
# Includes baseline comparison: sequence-only vs structural-only vs full.

library(xgboost)
library(ggplot2)

# ---------------------------------------------------------------------------
# Sample metadata extraction
# ---------------------------------------------------------------------------

#' Parse sample column names to extract experimental context.
#'
#' @param sample_col String, e.g. "1732_BE3_4h_DNA1_1.breakends.noEnds.FILTERED"
#' @return Named list with: grna_id, editor_type, timepoint, has_grna, is_replicate
parse_sample_metadata <- function(sample_col) {
  # Remove suffix
  name <- sub("\\.breakends\\.noEnds\\.FILTERED$", "", sample_col)

  # Determine gRNA
  if (grepl("^1732", name)) {
    grna_id <- "1732"
    has_grna <- 1L
    rest <- sub("^1732_", "", name)
  } else if (grepl("^4894", name)) {
    grna_id <- "4894"
    has_grna <- 1L
    rest <- sub("^4894_", "", name)
  } else {
    grna_id <- "none"
    has_grna <- 0L
    rest <- sub("^(Neg_|NC_)", "", name)
  }

  # Determine editor type
  if (grepl("^BE3x2", rest)) {
    editor_type <- "BE3x2"
  } else if (grepl("^BE3st", rest)) {
    editor_type <- "BE3st"
  } else if (grepl("^BE3", rest)) {
    editor_type <- "BE3"
  } else if (grepl("^Cas9", rest)) {
    editor_type <- "Cas9"
  } else {
    editor_type <- "none"
  }

  # Determine timepoint
  tp <- NA_real_
  if (grepl("1h", name)) tp <- 1
  else if (grepl("2h", name) || grepl("_2$", name) || grepl("_2\\.", name)) tp <- 2
  else if (grepl("4h", name)) tp <- 4
  else if (grepl("24h", name)) tp <- 24

  # Is replicate (has DNA identifier)
  is_replicate <- as.integer(grepl("DNA", name))

  list(
    grna_id = grna_id,
    editor_type = editor_type,
    timepoint = tp,
    has_grna = has_grna,
    is_replicate = is_replicate
  )
}

# ---------------------------------------------------------------------------
# Assemble long-format training data
# ---------------------------------------------------------------------------

#' Assemble long-format training data from wide-format feature matrix + sample columns.
#'
#' @param feat_df Data frame, one row per locus, with sequence features + structural features.
#' @param sample_cols Character vector, sample column names.
#' @param feature_groups List with named vectors of feature column names:
#'   $sequence (Group A), $alignment (Group B), $pair (Group C),
#'   $per_guide_struct (Group D), $global_struct (Group E)
#' @return Data frame in long format: one row per (locus, sample) pair.
assemble_long_data <- function(feat_df, sample_cols, feature_groups) {
  n_loci <- nrow(feat_df)

  # All feature columns (static per locus)
  all_feat_cols <- unlist(feature_groups)

  # Parse metadata for each sample
  meta_list <- lapply(sample_cols, parse_sample_metadata)

  rows <- list()
  kk <- 1

  for (j in seq_along(sample_cols)) {
    sc <- sample_cols[j]
    meta <- meta_list[[j]]
    counts <- feat_df[[sc]]

    # Build base data frame with locus features (repeated for each sample)
    base <- feat_df[, all_feat_cols, drop = FALSE]

    # Add experimental context (Group F)
    base$editor_type <- meta$editor_type
    base$timepoint <- meta$timepoint
    base$has_grna <- meta$has_grna
    base$grna_id <- meta$grna_id
    base$is_replicate <- meta$is_replicate
    base$sample_name <- sc

    # Add target
    base$edit_count <- counts
    base$log_edit_count <- log1p(counts)

    # Add locus ID
    base$locus_id <- feat_df$cluster_id

    rows[[kk]] <- base
    kk <- kk + 1
  }

  long_df <- do.call(rbind, rows)
  rownames(long_df) <- NULL

  # One-hot encode categorical variables
  long_df$editor_Cas9  <- as.integer(long_df$editor_type == "Cas9")
  long_df$editor_BE3   <- as.integer(long_df$editor_type == "BE3")
  long_df$editor_BE3st <- as.integer(long_df$editor_type == "BE3st")
  long_df$editor_BE3x2 <- as.integer(long_df$editor_type == "BE3x2")
  long_df$editor_none  <- as.integer(long_df$editor_type == "none")

  long_df$grna_1732 <- as.integer(long_df$grna_id == "1732")
  long_df$grna_4894 <- as.integer(long_df$grna_id == "4894")
  long_df$grna_none <- as.integer(long_df$grna_id == "none")

  # Remove character columns that shouldn't be in the feature matrix
  long_df$editor_type <- NULL
  long_df$grna_id <- NULL
  long_df$sample_name <- NULL

  long_df
}

# ---------------------------------------------------------------------------
# Define feature groups for baseline comparison
# ---------------------------------------------------------------------------

#' Define which columns belong to each feature group.
#'
#' @param long_df Long-format data frame.
#' @return List with character vectors of column names for each group.
define_feature_groups <- function(long_df) {
  cols <- colnames(long_df)

  # Group A: Sequence-derived (344-dim)
  group_a <- cols[cols %in% c("gc_content", "skew_AT", "skew_GC", "cpg_cnt",
                               grep("^freq_", cols, value = TRUE))]

  # Group B: Alignment (mismatch profile + mismatch_count + PAM + orientation)
  group_b <- cols[cols %in% c(grep("^mm_pos_", cols, value = TRUE),
                               "mismatch_count", "PAM_mismatch", "hit_orientation")]

  # Group C: On-target pair features
  group_c <- cols[cols %in% c("pair_seq_identity", "pair_flanking_kmer_similarity",
                               "pair_tree_cophenetic", "pair_same_chromosome",
                               "pair_genomic_distance", "pair_gc_difference",
                               "pair_cpg_difference")]

  # Group D: Per-guide structural (TDA + tree + latent)
  group_d <- cols[cols %in% c(
    grep("^H[01]_", cols, value = TRUE),       # persistence stats
    grep("^land_H", cols, value = TRUE),        # landscapes
    grep("^mapper_", cols, value = TRUE),       # Mapper features
    grep("^tree_", cols, value = TRUE),         # tree features (per-guide)
    grep("^pg_PC", cols, value = TRUE)          # per-guide latent projection
  )]

  # Group E: Global structural
  group_e <- cols[cols %in% c(
    grep("^global_H[01]_", cols, value = TRUE),
    grep("^global_land", cols, value = TRUE),
    grep("^global_tree_", cols, value = TRUE),
    grep("^global_PC", cols, value = TRUE)
  )]

  # Group F: Experimental context
  group_f <- cols[cols %in% c("editor_Cas9", "editor_BE3", "editor_BE3st",
                               "editor_BE3x2", "editor_none", "timepoint",
                               "has_grna", "grna_1732", "grna_4894", "grna_none",
                               "is_replicate")]

  # Group G: Chromatin accessibility + functional annotation
  group_g <- cols[cols %in% c(
    grep("^atac_", cols, value = TRUE),
    grep("^func_", cols, value = TRUE)
  )]

  # Group H: Chromatin × mismatch interaction features
  group_h <- cols[cols %in% c(
    grep("^intx_", cols, value = TRUE)
  )]

  list(
    sequence = group_a,
    alignment = group_b,
    pair = group_c,
    per_guide_struct = group_d,
    global_struct = group_e,
    experimental = group_f,
    chromatin = group_g,
    interaction = group_h
  )
}

# ---------------------------------------------------------------------------
# Train and evaluate XGBoost pair scorer
# ---------------------------------------------------------------------------

#' Train XGBoost regressor and evaluate.
#'
#' @param long_df Long-format data frame.
#' @param feature_cols Character vector, columns to use as features.
#' @param test_loci Integer vector, locus_ids for test set.
#' @param n_folds Integer, CV folds for hyperparameter tuning.
#' @param seed Integer, random seed.
#' @return List with: predictions, metrics, model
train_pair_scorer <- function(long_df, feature_cols, test_loci,
                               n_folds = 5, seed = 123) {
  set.seed(seed)

  # Split by locus
  train_df <- long_df[!long_df$locus_id %in% test_loci, ]
  test_df  <- long_df[long_df$locus_id %in% test_loci, ]

  # Prepare matrices
  X_train <- as.matrix(train_df[, feature_cols, drop = FALSE])
  y_train <- train_df$log_edit_count
  X_test  <- as.matrix(test_df[, feature_cols, drop = FALSE])
  y_test  <- test_df$log_edit_count

  # Handle NA: replace with 0 (XGBoost handles NA natively, but matrix conversion may introduce them)
  X_train[is.na(X_train)] <- 0
  X_test[is.na(X_test)] <- 0

  # XGBoost DMatrix
  dtrain <- xgb.DMatrix(data = X_train, label = y_train)
  dtest  <- xgb.DMatrix(data = X_test, label = y_test)

  # Hyperparameter grid (modest for "start simple")
  params_list <- list(
    list(max_depth = 4, eta = 0.1,  subsample = 0.8, colsample_bytree = 0.8),
    list(max_depth = 6, eta = 0.1,  subsample = 0.8, colsample_bytree = 0.8),
    list(max_depth = 8, eta = 0.05, subsample = 0.8, colsample_bytree = 0.6)
  )

  best_score <- Inf
  best_params <- NULL
  best_nrounds <- 100

  cat("  Hyperparameter tuning (", length(params_list), " configs, ",
      n_folds, "-fold CV)...\n", sep = "")

  for (i in seq_along(params_list)) {
    p <- params_list[[i]]
    p$objective <- "reg:squarederror"
    p$eval_metric <- "rmse"
    p$seed <- seed

    cv <- xgb.cv(
      params = p,
      data = dtrain,
      nrounds = 200,
      nfold = n_folds,
      early_stopping_rounds = 15,
      verbose = 0,
      prediction = TRUE
    )

    # Find best iteration (handle xgboost version differences)
    best_iter <- cv$best_iteration
    if (is.null(best_iter) || length(best_iter) == 0) {
      # Manual: find iteration with minimum test RMSE
      best_iter <- which.min(cv$evaluation_log$test_rmse_mean) - 1L
    }
    cv_score <- cv$evaluation_log$test_rmse_mean[best_iter + 1L]
    if (length(cv_score) == 0 || !is.finite(cv_score)) {
      cv_score <- min(cv$evaluation_log$test_rmse_mean, na.rm = TRUE)
      best_iter <- which.min(cv$evaluation_log$test_rmse_mean) - 1L
    }

    cat("    Config", i, ": depth=", p$max_depth, " eta=", p$eta,
        " -> RMSE=", round(cv_score, 4), " at round", best_iter, "\n")

    if (cv_score < best_score) {
      best_score <- cv_score
      best_params <- p
      best_nrounds <- best_iter
    }
  }

  cat("  Best config: depth=", best_params$max_depth, " eta=", best_params$eta,
      " nrounds=", best_nrounds, "\n")

  # Train final model on full training set
  watchlist <- list(train = dtrain, test = dtest)
  model <- xgb.train(
    params = best_params,
    data = dtrain,
    nrounds = best_nrounds + 1,
    watchlist = watchlist,
    verbose = 0
  )

  # Predictions
  pred_train <- predict(model, dtrain)
  pred_test  <- predict(model, dtest)

  # Metrics
  metrics <- compute_metrics(pred_test, y_test, test_df)

  list(
    model = model,
    predictions = data.frame(
      locus_id = test_df$locus_id,
      actual = y_test,
      predicted = pred_test,
      grna_1732 = test_df$grna_1732,
      grna_4894 = test_df$grna_4894,
      has_grna = test_df$has_grna
    ),
    metrics = metrics,
    best_params = best_params,
    best_nrounds = best_nrounds
  )
}

# ---------------------------------------------------------------------------
# Evaluation metrics
# ---------------------------------------------------------------------------

#' Compute evaluation metrics.
#'
#' @param pred Numeric vector, predicted log_edit_count.
#' @param actual Numeric vector, actual log_edit_count.
#' @param test_df Data frame, test set with metadata.
#' @return Named numeric vector.
compute_metrics <- function(pred, actual, test_df) {
  # Spearman correlation (primary metric, matches Kaufmann et al.)
  spearman_all <- cor(pred, actual, method = "spearman")

  # R²
  ss_res <- sum((actual - pred)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  r2_all <- 1 - ss_res / ss_tot

  # Per-guide Spearman
  spearman_1732 <- if (any(test_df$grna_1732 == 1)) {
    idx <- test_df$grna_1732 == 1
    cor(pred[idx], actual[idx], method = "spearman")
  } else NA_real_

  spearman_4894 <- if (any(test_df$grna_4894 == 1)) {
    idx <- test_df$grna_4894 == 1
    cor(pred[idx], actual[idx], method = "spearman")
  } else NA_real_

  # Recall@k: fraction of high-count loci (count > median nonzero) in top-k predictions
  # Per guide, average across samples
  recall_at_50 <- compute_recall_at_k(pred, actual, test_df, k = 50)
  recall_at_100 <- compute_recall_at_k(pred, actual, test_df, k = 100)

  c(
    spearman_all = spearman_all,
    r2_all = r2_all,
    spearman_1732 = spearman_1732,
    spearman_4894 = spearman_4894,
    recall_at_50 = recall_at_50,
    recall_at_100 = recall_at_100
  )
}

#' Compute Recall@k: fraction of true high-edit loci found in top-k predictions.
#'
#' @param pred Numeric vector, predictions.
#' @param actual Numeric vector, actual values.
#' @param test_df Data frame with locus_id and grna columns.
#' @param k Integer, number of top predictions to consider.
#' @return Numeric, mean recall across guides.
compute_recall_at_k <- function(pred, actual, test_df, k = 50) {
  recalls <- c()

  for (grna in c("1732", "4894")) {
    grna_col <- paste0("grna_", grna)
    if (!grna_col %in% colnames(test_df)) next
    idx <- test_df[[grna_col]] == 1
    if (sum(idx) == 0) next

    pred_g <- pred[idx]
    actual_g <- actual[idx]
    locus_g <- test_df$locus_id[idx]

    # Define "true high-edit" as above median of nonzero actual values
    nonzero <- actual_g[actual_g > 0]
    if (length(nonzero) < 5) next
    threshold <- median(nonzero)
    true_high <- locus_g[actual_g > threshold]

    # Top-k predicted loci (by locus, using max prediction across samples)
    locus_pred <- tapply(pred_g, locus_g, max)
    top_k_loci <- names(sort(locus_pred, decreasing = TRUE))[1:min(k, length(locus_pred))]

    recall <- sum(true_high %in% top_k_loci) / length(true_high)
    recalls <- c(recalls, recall)
  }

  if (length(recalls) == 0) return(NA_real_)
  mean(recalls)
}

# ---------------------------------------------------------------------------
# Baseline comparison: sequence-only vs structural-only vs full
# ---------------------------------------------------------------------------

#' Run baseline comparison with five feature configurations.
#'
#' Configurations:
#'   1. sequence_only: sequence + alignment + experimental
#'   2. structural_only: pair + per_guide_struct + global_struct + experimental
#'   3. chromatin_only: chromatin + experimental
#'   4. full: all groups except chromatin + experimental
#'   5. full_chromatin: all groups including chromatin + experimental
#'
#' @param long_df Long-format data frame.
#' @param groups Feature group definitions (from define_feature_groups).
#' @param test_loci Integer vector, locus_ids for test set.
#' @param seed Integer, random seed.
#' @return Data frame with metrics for each model.
run_baseline_comparison <- function(long_df, groups, test_loci, seed = 123) {
  configs <- list(
    list(
      name = "sequence_only",
      cols = c(groups$sequence, groups$alignment, groups$experimental)
    ),
    list(
      name = "structural_only",
      cols = c(groups$pair, groups$per_guide_struct, groups$global_struct, groups$experimental)
    ),
    list(
      name = "chromatin_only",
      cols = c(groups$chromatin, groups$experimental)
    ),
    list(
      name = "full",
      cols = c(groups$sequence, groups$alignment, groups$pair,
               groups$per_guide_struct, groups$global_struct, groups$experimental)
    ),
    list(
      name = "full_chromatin",
      cols = c(groups$sequence, groups$alignment, groups$pair,
               groups$per_guide_struct, groups$global_struct,
               groups$chromatin, groups$experimental)
    ),
    list(
      name = "full_chromatin_intx",
      cols = c(groups$sequence, groups$alignment, groups$pair,
               groups$per_guide_struct, groups$global_struct,
               groups$chromatin, groups$interaction, groups$experimental)
    )
  )

  results <- list()

  for (cfg in configs) {
    cat("\n=== Training:", cfg$name, "(", length(cfg$cols), " features) ===\n")
    res <- train_pair_scorer(
      long_df = long_df,
      feature_cols = cfg$cols,
      test_loci = test_loci,
      seed = seed
    )

    res_row <- data.frame(
      model = cfg$name,
      n_features = length(cfg$cols),
      t(res$metrics),
      stringsAsFactors = FALSE
    )
    results[[cfg$name]] <- list(metrics_row = res_row, full_result = res)
    cat("  Spearman:", round(res$metrics["spearman_all"], 4),
        " R²:", round(res$metrics["r2_all"], 4), "\n")
  }

  list(
    summary = do.call(rbind, lapply(results, `[[`, "metrics_row")),
    details = results
  )
}
