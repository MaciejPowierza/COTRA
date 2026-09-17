# 04b_rgcca_blocks.R — build RGCCA/SGCCA blocks from the existing feature-routing system

sanitize_rgcca_block <- function(M) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  M[!is.finite(M)] <- NA_real_

  for (j in seq_len(ncol(M))) {
    col <- M[, j]
    m <- mean(col, na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    col[is.na(col)] <- m
    M[, j] <- col
  }

  sds <- apply(M, 2, sd)
  keep <- is.finite(sds) & sds > 0
  M <- M[, keep, drop = FALSE]

  if (ncol(M) == 0L) return(NULL)
  M
}

build_rgcca_blocks_from_routing <- function(X,
                                            spearman_minClusterSize = 6,
                                            jaccard_minClusterSize = 10,
                                            hellinger_minClusterSize = 10,
                                            include_route_prefix = TRUE) {
  X <- as.data.frame(X)

  routing_df <- route_features(X)

  b_sp <- build_blocks_spearman(X, routing_df, minClusterSize = spearman_minClusterSize)
  b_j  <- build_blocks_jaccard(X, routing_df, minClusterSize = jaccard_minClusterSize)
  b_h  <- build_blocks_hellinger(X, routing_df, minClusterSize = hellinger_minClusterSize)

  blocks_sp <- b_sp$blocks
  blocks_j  <- b_j$blocks
  blocks_h  <- b_h$blocks

  if (include_route_prefix) {
    if (length(blocks_sp) > 0) names(blocks_sp) <- paste0("sp_", seq_along(blocks_sp))
    if (length(blocks_j)  > 0) names(blocks_j)  <- paste0("jac_", seq_along(blocks_j))
    if (length(blocks_h)  > 0) names(blocks_h)  <- paste0("hel_", seq_along(blocks_h))
  }

  blocks_all <- c(blocks_sp, blocks_j, blocks_h)
  blocks_all <- blocks_all[!vapply(blocks_all, is.null, logical(1))]
  blocks_all <- blocks_all[vapply(blocks_all, length, integer(1)) > 0]

  block_list <- list()
  for (bn in names(blocks_all)) {
    feats <- intersect(blocks_all[[bn]], colnames(X))
    if (length(feats) == 0L) next
    block_list[[bn]] <- sanitize_rgcca_block(X[, feats, drop = FALSE])
  }

  block_list <- block_list[!vapply(block_list, is.null, logical(1))]

  list(
    routing_df = routing_df,
    blocks = blocks_all,
    block_list = block_list
  )
}

make_rgcca_connection <- function(block_list, mode = c("full","within_route_sparse")) {
  mode <- match.arg(mode)
  J <- length(block_list)

  if (J == 0L) stop("No blocks available for RGCCA/SGCCA.")

  C <- matrix(1, J, J)
  diag(C) <- 0
  rownames(C) <- colnames(C) <- names(block_list)

  if (mode == "within_route_sparse") {
    prefixes <- sub("_.*$", "", names(block_list))
    for (i in seq_len(J)) {
      for (j in seq_len(J)) {
        if (i == j) next
        # connect all, but give stronger interpretability to same-route blocks if needed later
        # current implementation keeps binary connectivity
        C[i, j] <- 1
      }
    }
  }

  C
}

#===Beginning of the insertion, 10.04===

extract_rgcca_scores <- function(fit,
                                 use = c("first","first2","all"),
                                 prefix = "RG") {
  use <- match.arg(use)

  Ylist <- fit$Y
  if (is.null(Ylist) || length(Ylist) == 0L) {
    stop("RGCCA returned no block components.")
  }

  if (is.null(names(Ylist))) {
    names(Ylist) <- paste0("block", seq_along(Ylist))
  }

  out <- lapply(seq_along(Ylist), function(j) {
    Yj <- as.matrix(Ylist[[j]])
    if (is.null(dim(Yj))) {
      Yj <- matrix(Yj, ncol = 1)
    }

    if (use == "first") {
      kk <- min(1L, ncol(Yj))
    } else if (use == "first2") {
      kk <- min(2L, ncol(Yj))
    } else {
      kk <- ncol(Yj)
    }

    Yj <- Yj[, seq_len(kk), drop = FALSE]
    colnames(Yj) <- paste0(prefix, "_", names(Ylist)[j], "_C", seq_len(ncol(Yj)))
    Yj
  })

  do.call(cbind, out)
}

#===End of the insertion, 10.04===

#extract_rgcca_scores <- function(fit,
#                                 use = c("all","first","average"),
#                                 prefix = "RG") {
#  use <- match.arg(use)

#  Ylist <- fit$Y
#  if (is.null(Ylist) || length(Ylist) == 0L) stop("RGCCA returned no block components.")

  # standardize names
#  if (is.null(names(Ylist))) names(Ylist) <- paste0("block", seq_along(Ylist))

#  if (use == "first") {
#    out <- lapply(seq_along(Ylist), function(j) {
#      Yj <- as.matrix(Ylist[[j]])
#      Yj <- Yj[, 1, drop = FALSE]
#      colnames(Yj) <- paste0(prefix, "_", names(Ylist)[j], "_C1")
#      Yj
#    })
#    return(do.call(cbind, out))
#  }

#  if (use == "average") {
#    out <- lapply(seq_along(Ylist), function(j) {
#      Yj <- as.matrix(Ylist[[j]])
#      avg <- rowMeans(Yj, na.rm = TRUE)
#      avg <- matrix(avg, ncol = 1)
#      colnames(avg) <- paste0(prefix, "_", names(Ylist)[j], "_AVG")
#      avg
#    })
#    return(do.call(cbind, out))
#  }

#  out <- lapply(seq_along(Ylist), function(j) {
#    Yj <- as.matrix(Ylist[[j]])
#    colnames(Yj) <- paste0(prefix, "_", names(Ylist)[j], "_C", seq_len(ncol(Yj)))
#    Yj
#  })

#  do.call(cbind, out)
#}
