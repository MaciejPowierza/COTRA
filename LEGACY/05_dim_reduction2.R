
# 05_dim_reduction.R — PCA variants: plain, weighted, HPCA

#===Beginning of the insertion, 17.02===

# 05_dim_reduction.R — PCA variants: plain, weighted, HPCA

dim_reduce <- function(df_feat,
                       method = c("none","plain_pca","weighted_pca","hpca", "nmf", "rgcca", "sgcca", "bpca", "fa", "lpe", "ica",
                                  "topvar", "lscore", "fosmod", "mcfs", "block_reps"),
                       weight_col = NULL,
                       k_dr = 2,
                       X_input = NULL,
                       feature_names_input = NULL,
		       hpca_feature_mode = c("super", "super_blocks"),
		       nmf_nonneg_mode = c("shift","drop_negative_features"),
                       nmf_feature_subset = c("all","nonnegative_only"),
                       nmf_seed = 123,
                       preproc_mode = c("none","basic","aggressive"),
                       corr_cutoff = 0.95,
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
                       ica_tol = 1e-4,
                       fs_preprocess = "center",
                       fs_graph_type = c("knn", 10),
                       fs_lscore_t = 1,
                       fs_mcfs_K = NULL,
                       fs_mcfs_lambda = 1,
                       fs_mcfs_t = 10,
                       block_rep_per_block = 1,
                       block_rep_stat = "variance") {

  method <- match.arg(method)
  hpca_feature_mode <- match.arg(hpca_feature_mode)
  nmf_nonneg_mode <- match.arg(nmf_nonneg_mode)
  nmf_feature_subset <- match.arg(nmf_feature_subset)
 
  if (method %in% c("rgcca", "sgcca")) {
    rgcca_scores_use <- match.arg(rgcca_scores_use)
  }

###===CUT AND PASTE INTO THE run_one_pass AS THE FIRST FACTOR EXAMINED===  

  preproc_mode <- match.arg(preproc_mode)

  prep <- preprocess_feature_matrix(
    df_feat = df_feat,
    preproc_mode = preproc_mode,
    corr_cutoff = corr_cutoff
    )

  df_feat2 <- prep$df_feat
  feature_set <- get_feature_set(df_feat2)
  diag_row <- compute_feature_diagnostics(df_feat2, feature_set = feature_set)

###===END OF CUT AND PASTE===

  take_k <- function(M, k_dr, prefix = "PC") {
    if (is.null(M)) stop("dim_reduce: PCA scores object is NULL")
    kk <- min(k_dr, ncol(M))
    M2 <- M[, seq_len(kk), drop = FALSE]
    colnames(M2) <- paste0("PC", seq_len(ncol(M2)))
    M2
  }

  #===Beginning of the insertion, 22.04===

  finalize_selected_features <- function(Xsel, k_req, method_name, prep, diag_row, model = NULL, extra = list()) {
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

    if (ncol(Xsel) == 0L) {
      stop(sprintf("%s: no valid selected features after sanitization.", method_name))
    }

    Xsel <- scale(Xsel, center = TRUE, scale = TRUE)
    Xsel[!is.finite(Xsel)] <- 0

    out <- list(
      pcs = Xsel,
      model = model,
      k_req = k_req,
      k_target = min(k_req, ncol(Xsel)),
      k_eff = ncol(Xsel),
      selected_features = colnames(Xsel),
      preprocessing = prep,
      feature_diagnostics = diag_row
    )

    if (length(extra) > 0) {
      out <- c(out, extra)
    }
    out
  }

  #===End of the insertion, 22.04===

  #===NO PROCESSING===
  
  #Beginning of the insertion, 09.03
  if (method == "none") {
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    # remove constant / invalid columns, like in plain_pca
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    # optional but strongly recommended: center/scale for full-dimensional detectors
    X2 <- scale(X2, center = TRUE, scale = TRUE)
    X2[!is.finite(X2)] <- 0
	
	print("X2")
	print(head(X2))

    return(list(
      pcs = X2,
      model = NULL,
      k_req = ncol(X2),
      k_target = ncol(X2),
      k_eff = ncol(X2),
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }
  #End of the insertion, 09.03


  #===DIMENSIONALITY REDUCTION===

  if (method == "plain_pca") {
    #print("aaapa!")
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    sds <- apply(X, 2, sd, na.rm = TRUE)
    const_cols <- names(sds)[is.na(sds) | sds == 0]
    if (length(const_cols) > 0) print(const_cols)

    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    pr <- prcomp(X2, center = TRUE, scale. = TRUE)
    pcs <- take_k(pr$x, k_dr)
	#print("PCs inside dim_reduce - plain_PCA")
	#print(head(pcs))
    return(list(pcs = pcs, 
                model = pr,
                preprocessing = prep,
                feature_diagnostics = diag_row))
  }

  if (method == "bpca") {
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    # drop constant / invalid columns, like plain_pca
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    if (ncol(X2) < 2) stop("bpca: too few non-constant features after filtering.")
    if (nrow(X2) < 2) stop("bpca: too few rows.")

    # Rdimtools handles preprocessing internally if requested.
    fit <- Rdimtools::do.bpca(
      X2,
      ndim = min(k_dr, ncol(X2), nrow(X2)),
      maxiter = 500,
      reltol = 1e-4
    )

    pcs <- as.matrix(fit$Y)

    print("head of pcs")
    print(head(pcs))

    pcs[!is.finite(pcs)] <- 0

    sds_pcs <- apply(pcs, 2, sd, na.rm = TRUE)
    keep_pcs <- is.finite(sds_pcs) & sds_pcs > 0
    pcs <- pcs[, keep_pcs, drop = FALSE]

    if (ncol(pcs) == 0) {
      stop("bpca: embedding has no valid columns after sanitization.")
    }
 
    if (ncol(pcs) > k_dr) {
      pcs <- pcs[, seq_len(k_dr), drop = FALSE]
    }
    colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))

    return(list(
      pcs = pcs,
      model = fit,
      k_req = k_dr,
      k_target = min(k_dr, ncol(X2), nrow(X2)),
      k_eff = ncol(pcs),
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }

  if (method == "fa") {
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    # drop constant / invalid columns, like plain_pca and bpca
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    if (ncol(X2) < 2) stop("fa: too few non-constant features after filtering.")
    if (nrow(X2) < 2) stop("fa: too few rows.")

    fit <- Rdimtools::do.fa(
      X2,
      ndim = min(k_dr, ncol(X2), nrow(X2)),
      maxiter = 50,
      tolerance = 1e-8
    )

    pcs <- as.matrix(fit$Y)

    # sanitize embedding
    pcs[!is.finite(pcs)] <- NA_real_

    sds_pcs <- apply(pcs, 2, sd, na.rm = TRUE)
    keep_pcs <- is.finite(sds_pcs) & sds_pcs > 0
    pcs <- pcs[, keep_pcs, drop = FALSE]

    if (ncol(pcs) == 0) {
      stop("fa: embedding has no valid columns after sanitization.")
    }

    if (ncol(pcs) > k_dr) {
      pcs <- pcs[, seq_len(k_dr), drop = FALSE]
    }

    colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))

    return(list(
      pcs = pcs,
      model = fit,
      k_req = k_dr,
      k_target = min(k_dr, ncol(X2), nrow(X2)),
      k_eff = ncol(pcs),
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }

  if (method == "lpe") {
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    # drop constant / invalid columns
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    if (ncol(X2) < 2) stop("lpe: too few non-constant features after filtering.")
    if (nrow(X2) < 3) stop("lpe: too few rows.")

    ndim_target <- min(k_dr, ncol(X2), nrow(X2))
    #numk_used <- if (is.null(lpe_numk)) max(ceiling(nrow(X2) / 10), 2) else as.integer(lpe_numk)
    numk_used <- if (is.null(lpe_numk)) min(2L, nrow(X2) - 1L) else as.integer(lpe_numk)
    numk_used <- max(2L, min(numk_used, nrow(X2) - 1L))
    print("NUMK_USED")
    print(numk_used)

    fit <- Rdimtools::do.lpe(
      X2,
      ndim = ndim_target,
      preprocess = lpe_preprocess,
      numk = numk_used
    )

    pcs <- as.matrix(fit$Y)
    pcs[!is.finite(pcs)] <- NA_real_

    sds_pcs <- apply(pcs, 2, sd, na.rm = TRUE)
    keep_pcs <- is.finite(sds_pcs) & sds_pcs > 0
    pcs <- pcs[, keep_pcs, drop = FALSE]

    if (ncol(pcs) == 0) {
      stop("lpe: embedding has no valid columns after sanitization.")
    }

    if (ncol(pcs) > k_dr) {
      pcs <- pcs[, seq_len(k_dr), drop = FALSE]
    }

    colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))

    return(list(
      pcs = pcs,
      model = fit,
      k_req = k_dr,
      k_target = ndim_target,
      k_eff = ncol(pcs),
      lpe_preprocess = lpe_preprocess,
      lpe_numk = numk_used,
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }
  
  if (method == "weighted_pca") {
    if (is.null(weight_col)) stop("weighted_pca requires weight_col.")
    if (!(weight_col %in% names(df_feat2))) {
      stop(sprintf("weighted_pca: weight_col '%s' not found. Available: %s",
                   weight_col, paste(head(names(df_feat2), 30), collapse = ", ")))
    }

    w <- suppressWarnings(as.numeric(df_feat2[[weight_col]]))
    if (sum(is.finite(w)) == 0) stop("weighted_pca: weights have no finite values.")
    w[!is.finite(w)] <- NA_real_

    keep <- !is.na(w) & w > 0
    if (sum(keep) < 3) stop("weighted_pca: too few positive weights.")

    X <- df_feat2[keep, feature_set, drop = FALSE]
    w <- w[keep]

    ncp <- min(k_dr, nrow(X) - 1, ncol(X))
    pca <- FactoMineR::PCA(X, row.w = w, scale.unit = TRUE, ncp = ncp, graph = FALSE)

    pcs <- take_k(pca$ind$coord, k_dr)
    return(list(pcs = pcs, 
                model = pca,
                preprocessing = prep,
                feature_diagnostics = diag_row))
}

  #Beginning of the insertion, 26.03===
  
  if (method == "nmf") {
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    # keep finite, non-constant columns
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    # minimal nonnegativity fix:
    # shift only columns that contain negatives
    mins <- apply(X2, 2, min, na.rm = TRUE)
    for (j in seq_len(ncol(X2))) {
      if (is.finite(mins[j]) && mins[j] < 0) X2[, j] <- X2[, j] - mins[j]
    }

    X2[!is.finite(X2)] <- 0

    rank_k <- min(k_dr, nrow(X2), ncol(X2))
    fit <- NMF::nmf(X2, rank = rank_k)

    pcs <- take_k(NMF::basis(fit), rank_k)
	if (nrow(pcs) != nrow(df_feat2)) {
        stop(sprintf("NMF scores have %d rows, but df_feat has %d rows. Check whether basis() or coef() should be used.",
             nrow(pcs), nrow(df_feat2)))
}
    return(list(
      pcs = pcs,
      model = fit,
      k_req = k_dr,
      k_target = rank_k,
      k_eff = ncol(pcs),
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }

  #End of the insertion, 26.03===

  #Beginning of the insertion, 14.04===

  if (method == "ica") {
    X <- as.matrix(df_feat2[, feature_set, drop = FALSE])

    # drop constant / invalid columns
    sds <- apply(X, 2, sd, na.rm = TRUE)
    keep <- !is.na(sds) & sds > 0
    X2 <- X[, keep, drop = FALSE]

    if (ncol(X2) < 2) stop("ica: too few non-constant features after filtering.")
    if (nrow(X2) < 2) stop("ica: too few rows.")

    ncomp_target <- min(k_dr, ncol(X2), nrow(X2))

    fit <- fastICA::fastICA(
      X = X2,
      n.comp = ncomp_target,
      alg.typ = ica_alg_typ,
      fun = ica_fun,
      alpha = ica_alpha,
      method = ica_method,
      row.norm = ica_row_norm,
      maxit = ica_maxit,
      tol = ica_tol,
      verbose = FALSE
    )

    # Independent components are in S
    pcs <- as.matrix(fit$S)
    pcs[!is.finite(pcs)] <- NA_real_

    sds_pcs <- apply(pcs, 2, sd, na.rm = TRUE)
    keep_pcs <- is.finite(sds_pcs) & sds_pcs > 0
    pcs <- pcs[, keep_pcs, drop = FALSE]

    if (ncol(pcs) == 0) {
      stop("ica: embedding has no valid columns after sanitization.")
    }

    if (ncol(pcs) > k_dr) {
      pcs <- pcs[, seq_len(k_dr), drop = FALSE]
    }

    colnames(pcs) <- paste0("IC", seq_len(ncol(pcs)))

    return(list(
      pcs = pcs,
      model = fit,
      k_req = k_dr,
      k_target = ncomp_target,
      k_eff = ncol(pcs),
      ica_alg_typ = ica_alg_typ,
      ica_fun = ica_fun,
      ica_alpha = ica_alpha,
      ica_method = ica_method,
      ica_row_norm = ica_row_norm,
      ica_maxit = ica_maxit,
      ica_tol = ica_tol,
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }

  #End of the insertion, 14.04===

  if (method == "hpca") {
    X <- df_feat2[, feature_set, drop = FALSE]

    routing_df <- route_features(X)
    b_sp <- build_blocks_spearman(X, routing_df, minClusterSize = 6)
    b_j  <- build_blocks_jaccard(X, routing_df)
    b_h  <- build_blocks_hellinger(X, routing_df)

    blocks_all <- c(b_sp$blocks, b_j$blocks, b_h$blocks)
    names(blocks_all) <- paste0("block", seq_along(blocks_all))

    # Drop empty / NULL blocks defensively
    blocks_all <- blocks_all[!vapply(blocks_all, is.null, logical(1))]
    blocks_all <- blocks_all[vapply(blocks_all, length, integer(1)) > 0]

	# block sizes (number of variables/features per block)
	block_sizes <- vapply(blocks_all, length, integer(1))

	# global cap by sample size
	k_target <- min(k_dr, nrow(X) - 1, ncol(X))
	print("k_target")
	print(k_target)
	
	#hp <- run_block_hpca(X, blocks_all, ncomp_per_block = k_target, ncomp_per_block_adjusting = FALSE)

	#if (is.null(hp$X_pca)) stop("hpca produced no components")
	#	k_eff <- min(k_target, ncol(hp$X_pca))
	#	pcs <- hp$X_pca[, seq_len(k_eff), drop = FALSE]
	
	hp <- run_block_hpca(X, blocks_all, ncomp_per_block = k_target, ncomp_per_block_adjusting = FALSE)

    if (hpca_feature_mode == "super") {
        if (is.null(hp$X_pca)) stop("hpca produced no super components")
        k_eff <- min(k_target, ncol(hp$X_pca))
        pcs <- hp$X_pca[, seq_len(k_eff), drop = FALSE]
    } else {
	    print("inside the super_blocks section")
        if (is.null(hp$X_combined)) stop("hpca produced no combined components")
        # Keep ALL combined dims (super first, then block PCs). k still controls per-layer dimension via k_target/k_hpca.
        pcs <- hp$X_combined
        k_eff <- ncol(pcs)
    }

	return(list(pcs = pcs,
               model = hp, 
               k_req = k_dr, 
               k_target = k_target, 
               k_eff = k_eff, 
               hpca_feature_mode = hpca_feature_mode,
               preprocessing = prep,
               feature_diagnostics = diag_row))
  }
  #Beginning of the insertion, 08.04

  if (method %in% c("rgcca", "sgcca")) {
    X <- df_feat2[, feature_set, drop = FALSE]

    rg_blocks <- build_rgcca_blocks_from_routing(X)
    block_list <- rg_blocks$block_list

    if (length(block_list) < 2L) {
      stop("RGCCA/SGCCA requires at least two non-empty routed blocks.")
    }

    Cmat <- make_rgcca_connection(
      block_list,
      mode = rgcca_connection_mode
    )

    # component count per block: same target k, but capped by block dimension
    #ncomp_vec <- vapply(block_list, function(B) {
    #  max(1L, min(k, ncol(B), nrow(B)))
    #}, integer(1))

   ncomp_vec <- vapply(block_list, function(B) {
     as.integer(max(1L, min(as.integer(k_dr), nrow(B), ncol(B), qr(B)$rank)))
   }, integer(1))

    if (method == "rgcca") {
      fit <- RGCCA::rgcca(
        blocks = block_list,
        connection = Cmat,
        method = "rgcca",
        tau = rgcca_tau,
        ncomp = ncomp_vec,
        scheme = rgcca_scheme,
        scale = TRUE,
        superblock = rgcca_superblock,
        verbose = FALSE
      )
    } else {
      # SGCCA: sparsity can be scalar or vector
      sparsity_vec <- rep(rgcca_sparsity, length(block_list))
      fit <- RGCCA::rgcca(
        blocks = block_list,
        connection = Cmat,
        method = "sgcca",
        sparsity = sparsity_vec,
        ncomp = ncomp_vec,
        scheme = rgcca_scheme,
        scale = TRUE,
        superblock = rgcca_superblock,
        verbose = FALSE
      )
    }

    pcs <- extract_rgcca_scores(
      fit,
      use = rgcca_scores_use,
      prefix = toupper(method)
    )

    k_eff <- ncol(pcs)

    print("ncomp_vec")
    print(ncomp_vec)

    return(list(
      pcs = pcs,
      model = fit,
      routing_df = rg_blocks$routing_df,
      routed_blocks = rg_blocks$blocks,
      block_list = block_list,
      connection = Cmat,
      ncomp_vec = ncomp_vec,
      n_blocks_rgcca = length(ncomp_vec),
      mean_ncomp_rgcca = mean(ncomp_vec),
      median_ncomp_rgcca = median(ncomp_vec),
      min_ncomp_rgcca = min(ncomp_vec),
      max_ncomp_rgcca = max(ncomp_vec),
      sum_ncomp_rgcca = sum(ncomp_vec),
      k_req = k_dr,
      k_target = k_dr,
      k_eff = k_eff,
      rgcca_scores_use = rgcca_scores_use,
      preprocessing = prep,
      feature_diagnostics = diag_row
    ))
  }
 
  #End of the insertion, 08.04
}

#===End of the insertion, 17.02===

