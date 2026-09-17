# 09b_pair_features.R — Pair-specific feature computation
#
# Computes features for both pair types:
#   Pair Type 1: (gRNA, candidate site) → mismatch profile, PAM mismatch
#   Pair Type 2: (on-target, off-target) → sequence identity, flanking similarity,
#                genomic distance, GC/CpG differences
#
# These features augment the 344-dim sequence features and structural context
# to form the complete pair scorer input.

library(Biostrings)

# ---------------------------------------------------------------------------
# gRNA sequences (from the dataset description)
# ---------------------------------------------------------------------------

GRNA_SEQUENCES <- list(
  "1732" = list(
    spacer = "TGCTCGAGTGGGTCCCCGTG",
    pam    = "AGG",
    full   = "TGCTCGAGTGGGTCCCCGTGAGG"
  ),
  "4894" = list(
    spacer = "GAGGACGAGATGTAAGAGGC",
    pam    = "TGG",
    full   = "GAGGACGAGATGTAAGAGGCTGG"
  )
)

# ---------------------------------------------------------------------------
# Pair Type 1: (gRNA, candidate site) — mismatch profile
# ---------------------------------------------------------------------------

#' Compute position-wise mismatch profile between gRNA spacer and the best aligned substring.
#'
#' The best_substring in the data is the 20nt spacer alignment (no PAM).
#' We compare it against the 20nt gRNA spacer.
#'
#' @param best_substring Character, the 20nt substring from the data (best alignment to gRNA spacer).
#' @param grna_id String, "1732" or "4894".
#' @return Named integer vector (20-dim, 1=mismatch, 0=match).
compute_mismatch_profile <- function(best_substring, grna_id) {
  grna_spacer <- GRNA_SEQUENCES[[grna_id]]$spacer

  if (is.na(best_substring) || nchar(best_substring) < 20) {
    return(setNames(rep(NA_integer_, 20), paste0("mm_pos_", 1:20)))
  }

  # Compare 20nt substring to 20nt gRNA spacer
  sub20 <- substr(best_substring, 1, 20)
  grna20 <- substr(grna_spacer, 1, 20)

  # Position-wise comparison
  mm <- as.integer(strsplit(sub20, "")[[1]] != strsplit(grna20, "")[[1]])
  setNames(mm, paste0("mm_pos_", 1:20))
}

#' Compute PAM mismatch from total mismatch count and spacer mismatch sum.
#'
#' Since best_substring contains only the 20nt spacer (no PAM), we cannot
#' read PAM mismatches directly from the mismatch profile. Instead we derive
#' them: PAM_mismatch = 1 if total mismatch_count > spacer mismatch sum.
#'
#' @param mismatch_count Integer, total mismatches from the data (includes PAM).
#' @param spacer_mismatch_sum Integer, sum of mm_pos_1..mm_pos_20.
#' @return Integer (0 or 1).
compute_pam_mismatch <- function(mismatch_count, spacer_mismatch_sum) {
  as.integer(mismatch_count > spacer_mismatch_sum)
}

# ---------------------------------------------------------------------------
# Pair Type 2: (on-target, off-target) — similarity features
# ---------------------------------------------------------------------------

#' Compute on-target vs off-target pair features.
#'
#' @param offtarget_row Data frame row for the off-target locus.
#' @param ontarget_row Data frame row for the on-target locus.
#' @param offtarget_features Numeric vector, 344-dim sequence features for off-target.
#' @param ontarget_features Numeric vector, 344-dim sequence features for on-target.
#' @param grna_id String, "1732" or "4894".
#' @param tree_cophenetic Numeric, cophenetic distance from tree (if available).
#' @return Named numeric vector (7 values).
compute_ontarget_pair_features <- function(offtarget_row, ontarget_row,
                                            offtarget_features, ontarget_features,
                                            grna_id, tree_cophenetic = NA_real_) {
  grna_spacer <- GRNA_SEQUENCES[[grna_id]]$spacer

  # Sequence identity: fraction of matching positions in 20nt spacer alignment
  off_sub <- offtarget_row[[paste0("best_substring_", grna_id)]]
  if (is.na(off_sub) || nchar(off_sub) < 20) {
    seq_identity <- NA_real_
  } else {
    off20 <- substr(off_sub, 1, 20)
    grna20 <- substr(grna_spacer, 1, 20)
    seq_identity <- mean(strsplit(off20, "")[[1]] == strsplit(grna20, "")[[1]])
  }

  # Flanking k-mer similarity: cosine similarity of full sequence feature vectors
  # (originally 4-mer only, but preprocessing may remove features, so use all)
  cos_sim <- sum(offtarget_features * ontarget_features) /
             (sqrt(sum(offtarget_features^2)) * sqrt(sum(ontarget_features^2)))
  if (!is.finite(cos_sim)) cos_sim <- 0

  # Same chromosome
  same_chrom <- as.integer(offtarget_row$chrom == ontarget_row$chrom)

  # Genomic distance (if same chromosome)
  if (same_chrom == 1) {
    genomic_dist <- log10(abs(offtarget_row$start - ontarget_row$start) + 1)
  } else {
    genomic_dist <- 0
  }

  # GC content difference (feature 1 of 344-dim vector)
  gc_diff <- offtarget_features[1] - ontarget_features[1]

  # CpG count difference (feature 4 of 344-dim vector)
  cpg_diff <- offtarget_features[4] - ontarget_features[4]

  c(
    pair_seq_identity = seq_identity,
    pair_flanking_kmer_similarity = cos_sim,
    pair_tree_cophenetic = tree_cophenetic,
    pair_same_chromosome = same_chrom,
    pair_genomic_distance = genomic_dist,
    pair_gc_difference = gc_diff,
    pair_cpg_difference = cpg_diff
  )
}

# ---------------------------------------------------------------------------
# Master function: compute all pair features for a set of loci
# ---------------------------------------------------------------------------

#' Compute all pair features for every locus relative to a given gRNA.
#'
#' @param feat_df Data frame with all loci (output of feature_extraction + preprocessing).
#' @param feature_set Character vector, feature column names.
#' @param grna_id String, "1732" or "4894".
#' @param mismatch_col String, column name for mismatch count to this gRNA.
#' @param tree_cophenetic Numeric vector (n), cophenetic distances from per-guide tree.
#' @return Data frame (n rows) with:
#'   - 20 mismatch profile columns (mm_pos_1 .. mm_pos_20)
#'   - mismatch_count, PAM_mismatch, hit_orientation
#'   - 7 on-target pair feature columns
compute_all_pair_features <- function(feat_df, feature_set, grna_id,
                                       mismatch_col, tree_cophenetic = NULL) {
  n <- nrow(feat_df)
  best_sub_col <- paste0("best_substring_", grna_id)
  hit_orient_col <- paste0("hit_orientation_", grna_id)

  # Identify on-target
  ont_idx <- which(feat_df[[mismatch_col]] == 0)
  if (length(ont_idx) == 0) {
    warning("No on-target for gRNA ", grna_id, ". Using min-mismatch locus.")
    ont_idx <- which.min(feat_df[[mismatch_col]])[1]
  }
  ont_idx <- ont_idx[1]
  ont_row <- feat_df[ont_idx, ]
  ont_features <- as.numeric(ont_row[, feature_set, drop = FALSE])

  # On-target flanking features for pair comparison
  ont_feat_vec <- as.numeric(feat_df[ont_idx, feature_set, drop = FALSE])

  # --- Mismatch profiles (Pair Type 1) ---
  mm_profiles <- t(sapply(seq_len(n), function(i) {
    compute_mismatch_profile(feat_df[[best_sub_col]][i], grna_id)
  }))

  mm_count <- feat_df[[mismatch_col]]

  # PAM mismatch: derived from total mismatch_count minus spacer mismatch sum.
  # best_substring is 20nt (spacer only), so PAM mismatches are not in the profile.
  spacer_mm_sum <- rowSums(mm_profiles, na.rm = TRUE)
  pam_mm <- compute_pam_mismatch(mm_count, spacer_mm_sum)

  hit_orient <- ifelse(feat_df[[hit_orient_col]] == "+", 1, -1)

  # --- On-target pair features (Pair Type 2) ---
  cat("  Computing on-target pair features for", n, "loci...\n")
  pair_features <- t(sapply(seq_len(n), function(i) {
    off_feat_vec <- as.numeric(feat_df[i, feature_set, drop = FALSE])
    tcoph <- if (!is.null(tree_cophenetic)) tree_cophenetic[i] else NA_real_
    compute_ontarget_pair_features(
      offtarget_row = feat_df[i, ],
      ontarget_row = ont_row,
      offtarget_features = off_feat_vec,
      ontarget_features = ont_feat_vec,
      grna_id = grna_id,
      tree_cophenetic = tcoph
    )
  }))

  # Assemble
  result <- data.frame(
    mm_profiles,
    mismatch_count = mm_count,
    PAM_mismatch = pam_mm,
    hit_orientation = hit_orient,
    pair_features,
    stringsAsFactors = FALSE
  )

  result
}
