# 09d_interaction_features.R — Chromatin × mismatch interaction features
#
# Captures how chromatin accessibility modulates the relationship between
# sequence mismatch profile and editing efficiency. Computed per-guide
# (using the correct gRNA's mismatch profile) and stored in shared intx_*
# columns, same pattern as mm_pos_*.
#
# Three categories:
#   1. Mismatch position aggregations (seed / distal / PAM counts)
#   2. Product interactions (mismatch summaries × key chromatin features)
#   3. Ratio features (accessibility normalized by mismatch burden)
#   4. Conditional features (baseline accessibility, closed+high-mismatch)

# ---------------------------------------------------------------------------
# Compute interaction features for a set of loci (single guide)
# ---------------------------------------------------------------------------
#'
#' @param mm_profiles Matrix (n × 20) of binary mismatch indicators.
#' @param mismatch_count Integer vector (n), total mismatches.
#' @param pam_mismatch Integer vector (n), binary (any PAM mismatch).
#' @param chromatin_df Data frame (n rows) with atac_* and func_* columns.
#' @return Data frame (n rows × 20 cols) prefixed intx_.
compute_interaction_features <- function(mm_profiles, mismatch_count, pam_mismatch,
                                          chromatin_df) {
  n <- nrow(mm_profiles)

  # ---- Step 1: Mismatch position aggregations ----
  # Seed: PAM-proximal 8 nt of spacer (positions 13-20)
  seed_mm <- rowSums(mm_profiles[, 13:20, drop = FALSE], na.rm = TRUE)
  # Distal: PAM-distal 12 nt of spacer (positions 1-12)
  distal_mm <- rowSums(mm_profiles[, 1:12, drop = FALSE], na.rm = TRUE)
  # PAM: derived from total mismatch_count minus spacer mismatch sum
  # (best_substring is 20nt spacer only; PAM is not in the profile)
  spacer_mm_sum <- rowSums(mm_profiles, na.rm = TRUE)
  pam_mm_count <- pmax(mismatch_count - spacer_mm_sum, 0)

  # ---- Extract key chromatin features ----
  signal <- chromatin_df$atac_mean_signal_200bp
  peak_overlap <- chromatin_df$atac_any_rep_peak_overlap
  peak_distance <- chromatin_df$atac_min_peak_distance
  overlaps_exon <- chromatin_df$func_overlaps_exon
  overlaps_promoter <- chromatin_df$func_overlaps_promoter
  rep1_peak_signal <- chromatin_df$atac_rep1_peak_signal
  rep2_peak_signal <- chromatin_df$atac_rep2_peak_signal

  # Handle NULL / missing columns
  if (is.null(signal)) signal <- rep(0, n)
  if (is.null(peak_overlap)) peak_overlap <- rep(0, n)
  if (is.null(peak_distance)) peak_distance <- rep(0, n)
  if (is.null(overlaps_exon)) overlaps_exon <- rep(0, n)
  if (is.null(overlaps_promoter)) overlaps_promoter <- rep(0, n)
  if (is.null(rep1_peak_signal)) rep1_peak_signal <- rep(0, n)
  if (is.null(rep2_peak_signal)) rep2_peak_signal <- rep(0, n)

  # ---- Step 2: Product interactions (12) ----
  intx_mm_count_x_signal       <- mismatch_count * signal
  intx_mm_count_x_peak_overlap <- mismatch_count * peak_overlap
  intx_mm_count_x_peak_distance <- mismatch_count * peak_distance
  intx_seed_mm_x_signal        <- seed_mm * signal
  intx_seed_mm_x_peak_overlap  <- seed_mm * peak_overlap
  intx_seed_mm_x_peak_distance <- seed_mm * peak_distance
  intx_distal_mm_x_signal      <- distal_mm * signal
  intx_distal_mm_x_peak_overlap <- distal_mm * peak_overlap
  intx_pam_mm_x_signal         <- pam_mismatch * signal
  intx_pam_mm_x_peak_overlap   <- pam_mismatch * peak_overlap
  intx_mm_count_x_exonic       <- mismatch_count * overlaps_exon
  intx_mm_count_x_promoter     <- mismatch_count * overlaps_promoter

  # ---- Step 3: Ratio features (3) ----
  # Accessibility normalized by mismatch burden
  intx_signal_per_mm <- signal / (1 + mismatch_count)
  # Peak signal per mismatch (average of both reps)
  intx_peak_signal_per_mm <- (rep1_peak_signal + rep2_peak_signal) /
                              (2 * (1 + mismatch_count))
  # Accessibility weighted by match fraction
  intx_accessibility_x_match <- signal * (1 - mismatch_count / 20)

  # ---- Step 4: Conditional features (2) ----
  # Signal at on-target-like sites (mismatch_count == 0)
  intx_peak_mm0_signal <- signal * as.integer(mismatch_count == 0)
  # Closed chromatin + high mismatch (expected near-zero editing)
  intx_closed_high_mm <- (1 - peak_overlap) * as.integer(mismatch_count > 3)

  # ---- Assemble ----
  result <- data.frame(
    intx_seed_mismatch_count       = seed_mm,
    intx_distal_mismatch_count     = distal_mm,
    intx_pam_mismatch_count        = pam_mm_count,
    intx_mm_count_x_signal         = intx_mm_count_x_signal,
    intx_mm_count_x_peak_overlap   = intx_mm_count_x_peak_overlap,
    intx_mm_count_x_peak_distance  = intx_mm_count_x_peak_distance,
    intx_seed_mm_x_signal          = intx_seed_mm_x_signal,
    intx_seed_mm_x_peak_overlap    = intx_seed_mm_x_peak_overlap,
    intx_seed_mm_x_peak_distance   = intx_seed_mm_x_peak_distance,
    intx_distal_mm_x_signal        = intx_distal_mm_x_signal,
    intx_distal_mm_x_peak_overlap  = intx_distal_mm_x_peak_overlap,
    intx_pam_mm_x_signal           = intx_pam_mm_x_signal,
    intx_pam_mm_x_peak_overlap     = intx_pam_mm_x_peak_overlap,
    intx_mm_count_x_exonic         = intx_mm_count_x_exonic,
    intx_mm_count_x_promoter       = intx_mm_count_x_promoter,
    intx_signal_per_mm             = intx_signal_per_mm,
    intx_peak_signal_per_mm        = intx_peak_signal_per_mm,
    intx_accessibility_x_match     = intx_accessibility_x_match,
    intx_peak_mm0_signal           = intx_peak_mm0_signal,
    intx_closed_high_mm            = intx_closed_high_mm,
    stringsAsFactors = FALSE
  )

  # Ensure numeric
  for (col in names(result)) {
    result[[col]] <- as.numeric(result[[col]])
  }

  result
}
