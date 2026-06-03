
# 04_routing_blocks_hpca.R — feature routing + block building + HPCA
# Refactors the "insertion, 26.11." section. fileciteturn0file0

feature_gini <- function(x) ineq::Gini(x + 1e-12)

feature_entropy <- function(x) {
  x <- x[!is.na(x)]
  s <- sum(x)
  if (s <= 0) return(NA_real_)
  p <- x / s
  p <- p[p > 0]
  -sum(p * log(p))
}

assess_feature <- function(x) {
  N <- length(x)
  list(
    zero_fraction = sum(x == 0, na.rm = TRUE) / N,
    support       = sum(x != 0, na.rm = TRUE),
    var           = stats::var(x, na.rm = TRUE),
    sd            = stats::sd(x, na.rm = TRUE),
    mean          = mean(x, na.rm = TRUE),
    gini          = feature_gini(x),
    entropy       = feature_entropy(x)
  )
}

choose_feature_route <- function(stats, N,
                                 thr_zero_spearman_max = 0.6,
                                 thr_support_spearman_min = 0.4,
                                 thr_gini_spearman_max = 0.7,
                                 thr_zero_drop_min = 0.95,
                                 thr_var_min = 1e-8) {
  zfrac   <- stats$zero_fraction
  support <- stats$support
  gini    <- stats$gini
  var_x   <- stats$var

  if (is.na(var_x) || var_x < thr_var_min || zfrac >= thr_zero_drop_min) return("drop")

  if (!is.na(gini) &&
      zfrac   < thr_zero_spearman_max &&
      support > thr_support_spearman_min * N &&
      gini    < thr_gini_spearman_max) return("spearman")

  if (zfrac >= thr_zero_spearman_max) return("jaccard")
  "hellinger"
}

route_features <- function(X, ...) {
  X <- as.data.frame(X)
  N <- nrow(X)

  feats <- colnames(X)
  rows <- lapply(feats, function(f) {
    st <- assess_feature(X[[f]])
    route <- choose_feature_route(st, N, ...)
    data.frame(feature = f, route = route, st, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  res$route <- factor(res$route, levels = c("spearman","hellinger","jaccard","drop"))
  res
}

build_blocks_spearman <- function(X, routing_df, deepSplit = 0, minClusterSize = 4, do_plot = FALSE) {
  sel <- routing_df$route == "spearman"
  if (!any(sel)) return(list(blocks = list(), groups = NULL, hc = NULL))
  feats <- routing_df$feature[sel]
  X_sp <- X[, feats, drop = FALSE]

  cor_sp <- cor(X_sp, method = "spearman", use = "pairwise.complete.obs")
  D_sp <- as.dist(sqrt(2 - 2 * cor_sp))
  hc_sp <- hclust(D_sp, method = "average")
  if (do_plot) plot(hc_sp, main = "Spearman feature dendrogram")

  groups_sp <- dynamicTreeCut::cutreeDynamic(
    dendro = hc_sp, distM = as.matrix(D_sp),
    deepSplit = deepSplit, minClusterSize = minClusterSize
  )

  valid <- setdiff(unique(groups_sp), 0L)
  blocks <- split(colnames(X_sp), groups_sp)[as.character(valid)]
  list(blocks = blocks, groups = groups_sp, hc = hc_sp)
}

build_blocks_jaccard <- function(X, routing_df, deepSplit = 2, minClusterSize = 10) {
  sel <- routing_df$route == "jaccard"
  if (!any(sel)) return(list(blocks = list(), groups = NULL, hc = NULL))
  feats <- routing_df$feature[sel]
  X_j <- X[, feats, drop = FALSE]

  X_bin <- (X_j > 0) * 1
  D_j <- vegan::vegdist(t(X_bin), method = "jaccard")
  hc_j <- hclust(D_j, method = "average")

  groups_j <- dynamicTreeCut::cutreeDynamic(
    dendro = hc_j, distM = as.matrix(D_j),
    deepSplit = deepSplit, minClusterSize = minClusterSize
  )

  valid <- setdiff(unique(groups_j), 0L)
  blocks <- split(colnames(X_j), groups_j)[as.character(valid)]
  list(blocks = blocks, groups = groups_j, hc = hc_j)
}

build_blocks_hellinger <- function(X, routing_df, deepSplit = 1, minClusterSize = 10) {
  sel <- routing_df$route == "hellinger"
  if (sum(sel) < 2L) return(list(blocks = list(), groups = NULL, hc = NULL))
  feats <- routing_df$feature[sel]
  X_h <- X[, feats, drop = FALSE]

  X_hel <- vegan::decostand(X_h, "hellinger")
  D_h <- dist(t(X_hel))
  hc_h <- hclust(D_h, method = "average")

  groups_h <- dynamicTreeCut::cutreeDynamic(
    dendro = hc_h, distM = as.matrix(D_h),
    deepSplit = deepSplit, minClusterSize = minClusterSize
  )

  valid <- setdiff(unique(groups_h), 0L)
  blocks <- split(colnames(X_h), groups_h)[as.character(valid)]
  list(blocks = blocks, groups = groups_h, hc = hc_h)
}

#===Beginning of the insertion, 26.02===

# ---- Block similarity helpers (for block_list: list of n x p matrices) ----

block_similarity_rv <- function(block_list, center = TRUE, scale. = TRUE) {
  bn <- names(block_list)# bn: list of blocks' names, equal the number of blocks (e.g. 18)
  B <- length(block_list)# B: length of the block list; numeric (e.g. 18)
  S <- matrix(NA_real_, B, B, dimnames = list(bn, bn))# S: empty matrix filled with NA_real_s, with the dimensions length of the blocks' list x length of the blocks' list

  prep <- function(X) {
    X <- as.matrix(X)
    X[!is.finite(X)] <- NA_real_
    # simple NA handling: replace NA with column means (minimal + stable)
	#mean imputation
    for (j in seq_len(ncol(X))) {
      col <- X[, j]
      m <- mean(col, na.rm = TRUE)
      col[is.na(col)] <- m
      X[, j] <- col
    }
    scale(X, center = center, scale = scale.)# It becomes the component of the Z list (list of scaled matrices)
  }
  
  prep <- function(X) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    X[!is.finite(X)] <- NA_real_

    # mean imputation, but handle all-NA columns
    for (j in seq_len(ncol(X))) {
      col <- X[, j]
      m <- mean(col, na.rm = TRUE)
      if (!is.finite(m)) m <- 0  # <- kluczowe (all-NA => NaN)
      col[is.na(col)] <- m
      X[, j] <- col
    }

    # drop zero-variance / non-finite sd columns BEFORE scale()
    sds <- apply(X, 2, sd)
    keep <- is.finite(sds) & sds > 0
    X <- X[, keep, drop = FALSE]
    if (ncol(X) == 0) return(NULL)

    Z <- scale(X, center = TRUE, scale = TRUE)

    # final guard
    Z[!is.finite(Z)] <- 0
    Z
  }

  Z <- lapply(block_list, prep)# Z: scaled matrices from the block_list; list

  # Precompute within-block Frobenius norms of crossprod(Z)
  denom_part <- vapply(Z, function(A) sum(crossprod(A)^2), numeric(1))#numeric vector of length equal the number of blocks (e.g. 18)

  for (i in seq_len(B)) {
    S[i, i] <- 1
    for (j in seq.int(i + 1, B)) {#error: with i <- 18, seq.int produces the sequence: 19, 18; 19 is beyond limits
	  if(j <= B){
        num <- sum(crossprod(Z[[i]], Z[[j]])^2)  # || Z_i^T Z_j ||_F^2
        den <- sqrt(denom_part[i] * denom_part[j])

        rv <- if (den > 0) num / den else NA_real_
        S[i, j] <- rv
        S[j, i] <- rv}
	  else {
	    #print("ooopa!")
	    #num <- sum(crossprod(Z[[i]], Z[[j-1]])^2)  # || Z_i^T Z_j ||_F^2
        #den <- sqrt(denom_part[i] * denom_part[j-1])}
	  }
    }
  }
  S
}

block_similarity_subspace <- function(block_list, k, center = TRUE, scale. = TRUE) {
  print("INSIDE THE block_similarity_subspace FUNCTION")
  bn <- names(block_list)
  B <- length(block_list)
  S <- matrix(NA_real_, B, B, dimnames = list(bn, bn))

  #prep <- function(X) {#the same helper function as in the block_similarity_rv; OUTPUT: Z: scaled matrices from the block_list; list
  #  X <- as.matrix(X)
  #  X[!is.finite(X)] <- NA_real_
  #  for (j in seq_len(ncol(X))) {
  #    col <- X[, j]
  #    m <- mean(col, na.rm = TRUE)
  #    col[is.na(col)] <- m
  #    X[, j] <- col
  #  }
  #  scale(X, center = center, scale = scale.)
  #}
  
  prep <- function(X) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    X[!is.finite(X)] <- NA_real_

    # mean imputation, but handle all-NA columns
    for (j in seq_len(ncol(X))) {
      col <- X[, j]
      m <- mean(col, na.rm = TRUE)
      if (!is.finite(m)) m <- 0  # <- kluczowe (all-NA => NaN)
      col[is.na(col)] <- m
      X[, j] <- col
    }

    # drop zero-variance / non-finite sd columns BEFORE scale()
    sds <- apply(X, 2, sd)
    keep <- is.finite(sds) & sds > 0
    X <- X[, keep, drop = FALSE]
    if (ncol(X) == 0) return(NULL)

    Z <- scale(X, center = TRUE, scale = TRUE)

    # final guard
    Z[!is.finite(Z)] <- 0
    Z
  }

  # Orthonormal bases in SAMPLE space (n x k) via SVD: U[,1:k]
  Ulist <- lapply(block_list, function(X) {
    Z <- prep(X)
	if (is.null(Z)) return(NULL)
    kk <- min(k, nrow(Z) - 1L, ncol(Z))
    if (kk < 1L) return(NULL)
	if (any(!is.finite(Z))) {
      cat("Non-finite in Z:", sum(!is.finite(Z)), "\n")
      print(which(!is.finite(Z), arr.ind = TRUE)[1:10, , drop = FALSE])
	}
	sds <- apply(Z, 2, sd)
	if (any(!is.finite(sds)) || any(sds == 0)) {
      cat("Bad sds:", sum(!is.finite(sds) | sds == 0), "\n")
	}
    sv <- svd(Z, nu = kk, nv = 0)
	#print(sv$u[, seq_len(kk), drop = FALSE])
    sv$u[, seq_len(kk), drop = FALSE]
  })
  
  #print("Ulist")
  #print(Ulist)

  for (i in seq_len(B)) {
    if (is.null(Ulist[[i]])) next
    S[i, i] <- 1
    for (j in seq.int(i + 1, B)) {
	  if(j <= B){
        if (is.null(Ulist[[j]])) next
        Ui <- Ulist[[i]]
        Uj <- Ulist[[j]]
        kk <- min(ncol(Ui), ncol(Uj))
        # principal angles via singular values of Ui^T Uj
        sv <- svd(crossprod(Ui[, seq_len(kk), drop = FALSE],
                            Uj[, seq_len(kk), drop = FALSE]),
                  nu = 0, nv = 0)$d
        # similarity = mean(cos^2(theta)) = mean(sv^2)
        sim <- mean(pmin(1, pmax(0, sv))^2)
        S[i, j] <- sim
        S[j, i] <- sim
		#print("S[B,B]")
		#print(S[B,B])
		}
		else {
		#print("ooopa!")
		}
    }
  }
  S
}

#===End of the insertion, 26.02===

#===Beginning of the insertion, 26.02a===

merge_blocks_to_meet_k <- function(block_list, k_target, scale. = TRUE,
                                  sim_fun = block_similarity_rv) {
  n <- nrow(block_list[[1]])

  merge_once <- function(block_list) {
    p <- vapply(block_list, ncol, integer(1))
    b_small <- names(which.min(p))

    S <- sim_fun(block_list, center = TRUE, scale. = scale.)
    # exclude self
    srow <- S[b_small, , drop = TRUE]
    srow[b_small] <- -Inf
    b_best <- names(which.max(srow))

    new_name <- paste0(b_small, "__", b_best)
    new_mat <- cbind(block_list[[b_small]], block_list[[b_best]])

    # rebuild list
    keep <- setdiff(names(block_list), c(b_small, b_best))
    out <- block_list[keep]
    out[[new_name]] <- new_mat

    list(block_list = out, merged = c(b_small, b_best), new_name = new_name)
  }

  history <- list()

  repeat {
    p <- vapply(block_list, ncol, integer(1))
    ncomp_blockwise <- pmin(n - 1L, p)
    if (min(ncomp_blockwise) >= k_target) break
    if (length(block_list) <= 1L) break

    step <- merge_once(block_list)
    block_list <- step$block_list
    history[[length(history) + 1L]] <- step
  }

  list(block_list = block_list, history = history)
}

#===End of the insertion, 26.02a===

#===Beginning of the insertion, 18.02===

run_block_hpca <- function(X, blocks, ncomp_per_block = 2, scale. = TRUE, ncomp_per_block_adjusting = TRUE) {
  X <- as.data.frame(X)#great matrix, e.g. 3562x344
  #blocks: list of colnames grouped in blocks
  if (length(blocks) == 0) return(list(hpca_model = NULL, X_pca = NULL, block_list = NULL, ncomp = NULL))

  # Build block_list with only existing features; drop empty blocks
  block_list <- list()
  for (bname in names(blocks)) {
    feats <- intersect(blocks[[bname]], colnames(X))
    if (length(feats) == 0L){
	print(bname)
	next}
    block_list[[bname]] <- as.matrix(X[, feats, drop = FALSE])
  }
  if (length(block_list) == 0L) return(list(hpca_model = NULL, X_pca = NULL, block_list = NULL, ncomp = NULL))

  block_sizes <- vapply(block_list, NCOL, integer(1))#named list of blocks' sizes

  max_by_n <- max(1L, nrow(X) - 1L)#nrow(X) - 1, e.g. 3562-1 = 3561
  
  print("ncomp_per_block before adjusting")
  print(ncomp_per_block)
  k_hpca <- as.integer(ncomp_per_block)

  # --- Build ncomp aligned to block_list ---
  if (ncomp_per_block_adjusting == TRUE){
    print("logical value of the argument ncomp_per_block_adjusting")
    print(ncomp_per_block_adjusting)
    if (length(ncomp_per_block) == 1L) {
      # scalar: same target for all blocks, capped per block
      ncomp <- rep(as.integer(ncomp_per_block), length(block_list))
    } else {
      # vector: align by names if possible, otherwise require exact length
      if (!is.null(names(ncomp_per_block))) {
        ncomp <- as.integer(ncomp_per_block[names(block_list)])
        if (any(is.na(ncomp))) stop("ncomp_per_block names do not cover all block_list names")
      } else {
        if (length(ncomp_per_block) != length(block_list)) {
          stop(sprintf("ncomp_per_block must be length 1 or %d (num blocks); got %d",
                       length(block_list), length(ncomp_per_block)))
        }
        ncomp <- as.integer(ncomp_per_block)
      }
    }
  
    # per-block feasible
    ncomp_blockwise <- pmin(ncomp, block_sizes, max_by_n)

    # shared dimensionality
    k_hpca <- max(1L, min(ncomp_blockwise))}
  
  sim_rv <- block_similarity_rv(block_list, center = TRUE, scale. = scale.)
  #print("sim_rv")
  #print(sim_rv)
  sim_subspace <- block_similarity_subspace(block_list, k = k_hpca, center = TRUE, scale. = scale.)
  
  #print("sim_subspace")
  #print(sim_subspace)
  
  print("k_hpca")
  print(k_hpca)
  
  if (k_hpca > 1L) {
    merged <- merge_blocks_to_meet_k(block_list, k_target = k_hpca, scale. = scale.)
    block_list_2 <- merged$block_list
    merge_history <- merged$history
}

  sanitize_block <- function(M) {
    M <- as.matrix(M)
    storage.mode(M) <- "double"
    M[!is.finite(M)] <- NA_real_

    # impute NA per-column, but handle all-NA columns
    for (j in seq_len(ncol(M))) {
      col <- M[, j]
      m <- mean(col, na.rm = TRUE)
      if (!is.finite(m)) m <- 0  # all-NA => NaN => ustaw 0
      col[is.na(col)] <- m
      M[, j] <- col
    }

    # drop zero-variance / bad columns
    sds <- apply(M, 2, sd)
    keep <- is.finite(sds) & sds > 0
    M <- M[, keep, drop = FALSE]

    if (ncol(M) == 0L) return(NULL)
    M
  }
  
  block_list_2 <- lapply(block_list_2, sanitize_block)
  block_list_2 <- block_list_2[!vapply(block_list_2, is.null, logical(1))]

  if (length(block_list_2) == 0L) {
    stop("After sanitization, all HPCA blocks are empty.")
  }
  
  for (nm in names(block_list_2)) {
    M <- scale(block_list_2[[nm]], center = TRUE, scale = TRUE)
    ok <- tryCatch({ svd(M, nu = 0, nv = 1); TRUE }, error = function(e) e)
    if (!isTRUE(ok)) {
      cat("SVD FAILS in block:", nm, "\n")
      print(ok)
      break
    }
	else {
	  cat("SVD PASSES in block:", nm, "\n")
	}
  }
  
  scaled_blocks <- lapply(block_list_2, function(M) scale(M, TRUE, scale.))
  ranks <- vapply(scaled_blocks, function(M) qr(M)$rank, integer(1))
  sizes <- vapply(block_list_2, ncol, integer(1))
  df <- data.frame(block = names(ranks), ncol = sizes, rank = ranks)
  print(df[order(df$rank), ]) 

  hpca_model <- multiblock::hpca(X = block_list_2, ncomp = k_hpca, scale = scale., init = "random", verbose = FALSE)

  #X_pca <- if (!is.null(hpca_model$scores)) hpca_model$scores else hpca_model$Scores$Common
  X_pca <- tryCatch(scores(hpca_model, block = 0), error = function(e) NULL)
  
  print("X_pca")
  print(ncol(X_pca))
  
  X_blocks <- list()
  for (j in seq_len(length(block_list_2))) {
    Sj <- tryCatch(scores(hpca_model, block = j), error = function(e) NULL)
    if (!is.null(Sj) && ncol(Sj) > 0) {
      X_blocks[[names(block_list_2)[j]]] <- as.matrix(Sj)
    }
  }

  # scale block scores columnwise (recommended for anomaly detection)
  scale_block_scores <- TRUE
  #print("no of X_blocks")
  #print(length(X_blocks))
  print("Number of inner blocks")
  TEMP_diagnostics <- lapply(X_blocks, ncol)
  print(sum(unlist(TEMP_diagnostics)))
  if (length(X_blocks) > 0 && isTRUE(scale_block_scores)) {
    X_blocks <- lapply(X_blocks, function(M) {
      M <- as.matrix(M)
      M <- scale(M, center = TRUE, scale = TRUE)
      M[!is.finite(M)] <- 0
      M
    })
  }

  X_combined <- X_pca
  if (length(X_blocks) > 0) {
    X_combined <- cbind(X_pca, do.call(cbind, X_blocks))
  }

  # Give stable column names
  colnames(X_pca) <- paste0("PC", seq_len(ncol(X_pca)))
  colnames(X_combined) <- c(
    paste0("PC", seq_len(ncol(X_pca))),
    unlist(lapply(names(X_blocks), function(bn) paste0(bn, "_PC", seq_len(ncol(X_blocks[[bn]])))))
  )
  
  if (is.null(X_pca) || ncol(X_pca) == 0) {
    return(list(hpca_model = hpca_model, X_pca = NULL, block_list = block_list, ncomp = k_hpca))
  }

  #colnames(X_pca) <- paste0("PC", seq_len(ncol(X_pca)))
  list(hpca_model = hpca_model, X_pca = X_pca, X_blocks = X_blocks, X_combined = X_combined, block_list = block_list, ncomp = k_hpca, sim_rv = sim_rv, sim_subspace = sim_subspace)
}

#===End of the insertion, 18.02===
