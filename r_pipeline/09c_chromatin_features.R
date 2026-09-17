# 09c_chromatin_features.R — Chromatin accessibility + functional annotation features
#
# Extracts ATAC-seq chromatin accessibility signals and functional genomic
# annotation from columns appended to the input TSV. All features are
# per-locus (static across samples) and prefixed atac_ or func_ for clean
# group assignment in the pair scorer.
#
# Input columns (from the xlsx-derived TSV):
#   Functional annotation (43-53):
#     distance_to_nearest_gene, overlaps_gene, overlaps_exon,
#     overlaps_promoter, functional_class
#   ATAC-seq (54-67):
#     HEK293T_WT_ATAC_rep1_peak_overlap, rep1_nearest_peak_signalValue,
#     rep1_nearest_peak_distance, rep2_peak_overlap,
#     rep2_nearest_peak_signalValue, rep2_nearest_peak_distance,
#     any_rep_peak_overlap, both_reps_peak_overlap,
#     min_nearest_peak_distance, rep1_mean0_signal_200bp,
#     rep2_mean0_signal_200bp, mean0_signal_200bp_mean_reps

# ---------------------------------------------------------------------------
# Helper: convert mixed-type column to numeric
# ---------------------------------------------------------------------------
to_numeric <- function(x) {
  if (is.logical(x)) {
    return(as.numeric(x))
  }
  if (is.character(x)) {
    # Handle "True"/"False" strings
    x <- ifelse(x %in% c("True", "TRUE", "true"), "1",
         ifelse(x %in% c("False", "FALSE", "false"), "0", x))
    return(as.numeric(x))
  }
  as.numeric(x)
}

# ---------------------------------------------------------------------------
# Helper: convert boolean-like column to 0/1 integer
# ---------------------------------------------------------------------------
to_binary <- function(x) {
  if (is.logical(x)) {
    return(as.integer(x))
  }
  if (is.character(x)) {
    return(as.integer(x %in% c("True", "TRUE", "true", "1")))
  }
  as.integer(as.logical(x))
}

# ---------------------------------------------------------------------------
# Compute chromatin + functional annotation features
# ---------------------------------------------------------------------------
#'
#' @param df Data frame (feat_df) containing the ATAC and functional annotation
#'   columns from the xlsx-derived TSV.
#' @return Data frame with ~21 chromatin features, row-aligned with df.
#'   Columns prefixed atac_ or func_.
compute_chromatin_features <- function(df) {
  n <- nrow(df)
  chrom <- data.frame(row.names = seq_len(n))

  # ---- Functional annotation features ----

  # distance_to_nearest_gene → log1p
  if ("distance_to_nearest_gene" %in% names(df)) {
    d <- to_numeric(df$distance_to_nearest_gene)
    d[is.na(d)] <- 0
    chrom$func_distance_to_nearest_gene <- log1p(d)
  }

  # overlaps_gene, overlaps_exon, overlaps_promoter → 0/1
  for (col in c("overlaps_gene", "overlaps_exon", "overlaps_promoter")) {
    feat_name <- paste0("func_", col)
    if (col %in% names(df)) {
      chrom[[feat_name]] <- to_binary(df[[col]])
    } else {
      chrom[[feat_name]] <- 0L
    }
  }

  # functional_class → one-hot (4 categories)
  fc_classes <- c("intergenic", "genic_non_exonic", "exonic", "promoter")
  if ("functional_class" %in% names(df)) {
    fc <- df$functional_class
    for (cls in fc_classes) {
      chrom[[paste0("func_class_", cls)]] <- as.integer(fc == cls)
    }
  } else {
    for (cls in fc_classes) {
      chrom[[paste0("func_class_", cls)]] <- 0L
    }
  }

  # ---- ATAC-seq features ----

  # Peak overlap (binary)
  atac_binary_cols <- list(
    list(src = "HEK293T_WT_ATAC_rep1_peak_overlap",  dst = "atac_rep1_peak_overlap"),
    list(src = "HEK293T_WT_ATAC_rep2_peak_overlap",  dst = "atac_rep2_peak_overlap"),
    list(src = "HEK293T_WT_ATAC_any_rep_peak_overlap",  dst = "atac_any_rep_peak_overlap"),
    list(src = "HEK293T_WT_ATAC_both_reps_peak_overlap", dst = "atac_both_reps_peak_overlap")
  )
  for (pair in atac_binary_cols) {
    if (pair$src %in% names(df)) {
      chrom[[pair$dst]] <- to_binary(df[[pair$src]])
    } else {
      chrom[[pair$dst]] <- 0L
    }
  }

  # Nearest peak signal value → log1p (string → numeric)
  atac_signal_cols <- list(
    list(src = "HEK293T_WT_ATAC_rep1_nearest_peak_signalValue", dst = "atac_rep1_peak_signal"),
    list(src = "HEK293T_WT_ATAC_rep2_nearest_peak_signalValue", dst = "atac_rep2_peak_signal")
  )
  for (pair in atac_signal_cols) {
    if (pair$src %in% names(df)) {
      v <- to_numeric(df[[pair$src]])
      v[is.na(v)] <- 0
      chrom[[pair$dst]] <- log1p(v)
    } else {
      chrom[[pair$dst]] <- 0
    }
  }

  # Distances → log1p
  atac_dist_cols <- list(
    list(src = "HEK293T_WT_ATAC_rep1_nearest_peak_distance", dst = "atac_rep1_peak_distance"),
    list(src = "HEK293T_WT_ATAC_rep2_nearest_peak_distance", dst = "atac_rep2_peak_distance"),
    list(src = "HEK293T_WT_ATAC_min_nearest_peak_distance",  dst = "atac_min_peak_distance")
  )
  for (pair in atac_dist_cols) {
    if (pair$src %in% names(df)) {
      v <- to_numeric(df[[pair$src]])
      v[is.na(v)] <- 0
      chrom[[pair$dst]] <- log1p(v)
    } else {
      chrom[[pair$dst]] <- 0
    }
  }

  # 200bp mean signal → log1p (string → numeric)
  atac_200bp_cols <- list(
    list(src = "HEK293T_WT_ATAC_rep1_mean0_signal_200bp",          dst = "atac_rep1_signal_200bp"),
    list(src = "HEK293T_WT_ATAC_rep2_mean0_signal_200bp",          dst = "atac_rep2_signal_200bp"),
    list(src = "HEK293T_WT_ATAC_mean0_signal_200bp_mean_reps",     dst = "atac_mean_signal_200bp")
  )
  for (pair in atac_200bp_cols) {
    if (pair$src %in% names(df)) {
      v <- to_numeric(df[[pair$src]])
      v[is.na(v)] <- 0
      chrom[[pair$dst]] <- log1p(v)
    } else {
      chrom[[pair$dst]] <- 0
    }
  }

  # Ensure all columns are numeric
  for (col in names(chrom)) {
    chrom[[col]] <- as.numeric(chrom[[col]])
  }

  cat("  Chromatin features:", ncol(chrom), "columns\n")
  cat("  Functional:", sum(grepl("^func_", names(chrom))), "ATAC:",
      sum(grepl("^atac_", names(chrom))), "\n")

  chrom
}
