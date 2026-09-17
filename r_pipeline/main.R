# main.R — COTRA remodeled pipeline entry point
# CRISPR On-Target Ranking Architecture → Structural Off-Target Pair Scorer
#
# Pipeline: data → features → structural representations (TDA + tree + latent)
#           → pair features → XGBoost pair scorer → evaluation
#
# Edit the CONFIG section only.

# ============ CONFIGURATION ============
# EDIT THESE PATHS to match your local setup

# R library path (where you installed R packages, e.g. via .libPaths())
# Set to your conda R library or remove this line if packages are in default path
.libPaths("/path/to/your/R/library")

# Path to the TSV data file (output of 00_data_prep.py)
DATA_PATH <- "FINAL_RESULTS_1732_4894.functional_annotation.FINAL.with_HEK293T_WT_ATAC.tsv"

# Output directory for results (long_data.csv, test_loci.csv, etc.)
OUTPUT_DIR <- "../results"
# ========================================

# -------------------------
# CONFIG
# -------------------------
CONFIG <- list(
  data_path = DATA_PATH,
  sep = "\t",

  # gRNA definitions
  grnas = list(
    "1732" = list(
      mismatch_col = "mismatch_count_1732",
      sample_prefix = "1732"
    ),
    "4894" = list(
      mismatch_col = "mismatch_count_4894",
      sample_prefix = "4894"
    )
  ),

  # Sample column range (for identifying all sample columns)
  sample_col_start = "1732_1.breakends.noEnds.FILTERED",
  sample_col_end   = "Neg_ctrl_1h_DNA1.breakends.noEnds.FILTERED",

  # Feature extraction
  seq_col = "centered_sequence",

  # Structural representation parameters
  tda_n_landmarks = 200,
  tda_n_landscape_samples = 20,
  tda_mapper_params = list(num_intervals = 10, percent_overlap = 50, num_bins = 8),

  tree_cut_heights = c(0.3, 0.5, 0.7),

  # Latent projection
  latent_method = "hpca",  # "pca" or "hpca"
  n_pcs_per_guide = 16,
  n_pcs_global = 16,
  use_umap_global = FALSE,

  # Pair scorer
  test_frac = 0.2,
  cv_folds = 5,
  seed = 123,

  # Output
  output_dir = OUTPUT_DIR
)

# -------------------------
# Load modules
# -------------------------
source("00_setup.R")
source("01_io.R")
source("02_sampling.R")
source("03_features.R")
source("03b_preprocessing.R")
source("04_routing_blocks_hpca2.R")
source("04b_rgcca_blocks.R")
source("05_dim_reduction2.R")
source("05b_tda.R")
source("05c_distance_tree.R")
source("05d_global_structure.R")
source("09b_pair_features.R")
source("09c_chromatin_features.R")
source("09d_interaction_features.R")
source("09_pair_scorer.R")
source("09b_pair_scorer_viz.R")

# -------------------------
# Run
# -------------------------
set.seed(CONFIG$seed)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n========================================\n")
cat("COTRA Remodeled Pipeline\n")
cat("Structural Off-Target Pair Scorer\n")
cat("========================================\n\n")

# --- Step 1: Read data ---
cat("Step 1: Reading data...\n")
bed_df <- read_bed_like(CONFIG$data_path, sep = CONFIG$sep)
cat("  Loaded", nrow(bed_df), "loci\n")

# Identify sample columns
all_cols <- colnames(bed_df)
sample_cols <- all_cols[which(all_cols == CONFIG$sample_col_start):which(all_cols == CONFIG$sample_col_end)]
cat("  Found", length(sample_cols), "sample columns\n\n")

# --- Step 2: Feature extraction ---
cat("Step 2: Extracting sequence features (344-dim)...\n")
feat_df <- feature_extraction(bed_df, seq_col = CONFIG$seq_col, draw_plot = FALSE)
feature_set <- get_feature_set(feat_df)
cat("  Feature set:", length(feature_set), "features\n\n")

# --- Step 3: Preprocessing ---
cat("Step 3: Preprocessing feature matrix...\n")
prep <- preprocess_feature_matrix(feat_df, preproc_mode = "basic", corr_cutoff = 0.95)
feat_df <- prep$df_feat
feature_set <- get_feature_set(feat_df)
cat("  Features after preprocessing:", length(feature_set), "\n\n")

# --- Step 4: Per-guide structural representations ---
cat("Step 4: Computing per-guide structural representations...\n")

per_guide_struct <- list()

for (grna_id in names(CONFIG$grnas)) {
  cat("\n  --- gRNA", grna_id, "---\n")
  grna_cfg <- CONFIG$grnas[[grna_id]]
  mm_col <- grna_cfg$mismatch_col

  # Select loci relevant to this gRNA (nonzero in at least one matching sample)
  grna_samples <- sample_cols[grepl(paste0("^", grna_cfg$sample_prefix), sample_cols)]
  has_edits <- rowSums(feat_df[, grna_samples, drop = FALSE]) > 0
  # Also include on-target (mismatch=0)
  is_ontarget <- feat_df[[mm_col]] == 0
  guide_idx <- which(has_edits | is_ontarget)

  cat("  Loci for this guide:", length(guide_idx), "\n")

  if (length(guide_idx) < 10) {
    cat("  Too few loci, skipping\n")
    next
  }

  guide_feat <- feat_df[guide_idx, , drop = FALSE]
  X_guide <- as.matrix(guide_feat[, feature_set, drop = FALSE])

  # Total edits across this guide's samples (for Mapper coloring)
  guide_edits <- rowSums(guide_feat[, grna_samples, drop = FALSE])

  # a) Latent projection (hPCA or PCA)
  if (CONFIG$latent_method == "hpca") {
    cat("  a) Latent projection (hPCA)...\n")
    dr_result <- dim_reduce(df_feat = guide_feat, method = "hpca",
                            k_dr = CONFIG$n_pcs_per_guide)
    pg_latent <- dr_result$pcs
  } else {
    cat("  a) Latent projection (PCA)...\n")
    X_scaled <- scale(X_guide)
    X_scaled[!is.finite(X_scaled)] <- 0
    pca_guide <- prcomp(X_scaled, rank. = CONFIG$n_pcs_per_guide)
    pg_latent <- pca_guide$x[, seq_len(min(CONFIG$n_pcs_per_guide, ncol(pca_guide$x))), drop = FALSE]
  }
  colnames(pg_latent) <- paste0("pg_PC", seq_len(ncol(pg_latent)))

  # b) TDA features
  cat("  b) TDA features...\n")
  tda_res <- compute_tda_features(
    X = X_guide,
    edits = guide_edits,
    n_landmarks = CONFIG$tda_n_landmarks,
    n_landscape_samples = CONFIG$tda_n_landscape_samples,
    mapper_params = CONFIG$tda_mapper_params,
    seed = CONFIG$seed
  )
  tda_global_df <- expand_tda_global(tda_res, nrow(guide_feat))
  tda_mapper_df <- tda_res$mapper_features

  # c) Distance tree
  cat("  c) Distance tree...\n")
  tree_res <- build_distance_tree(
    df = guide_feat,
    mismatch_col = mm_col,
    feature_cols = feature_set,
    cut_heights = CONFIG$tree_cut_heights
  )
  tree_features <- tree_res$features

  # Store per-guide results
  per_guide_struct[[grna_id]] <- list(
    indices = guide_idx,
    latent = pg_latent,
    tda_global = tda_global_df,
    tda_mapper = tda_mapper_df,
    tree = tree_features,
    tree_cophenetic = tree_res$cophenetic_dist,
    tda_full = tda_res
  )
}

cat("\n")

# --- Step 5: Global structural representations ---
cat("Step 5: Computing global structural representations...\n")

# Total edits across all samples for Mapper coloring
all_edits <- rowSums(feat_df[, sample_cols, drop = FALSE])

global_struct <- compute_global_structure(
  feat_df = feat_df,
  feature_set = feature_set,
  edits = all_edits,
  n_landmarks = CONFIG$tda_n_landmarks,
  n_landscape_samples = CONFIG$tda_n_landscape_samples,
  n_pcs = CONFIG$n_pcs_global,
  use_umap = CONFIG$use_umap_global,
  seed = CONFIG$seed,
  latent_method = CONFIG$latent_method
)
cat("\n")

# --- Step 6: Merge structural features into feat_df ---
cat("Step 6: Merging structural features...\n")

# Initialize structural columns with NA
n_total <- nrow(feat_df)

# Global features (same for all loci)
for (col in colnames(global_struct$tda_global)) {
  feat_df[[paste0("global_", col)]] <- global_struct$tda_global[[col]]
}
for (col in colnames(global_struct$tree_global)) {
  feat_df[[paste0("global_", col)]] <- global_struct$tree_global[[col]]
}
for (col in colnames(global_struct$latent_global)) {
  feat_df[[col]] <- global_struct$latent_global[[col]]
}

# Per-guide features (only for loci belonging to that guide)
for (grna_id in names(per_guide_struct)) {
  pg <- per_guide_struct[[grna_id]]
  idx <- pg$indices

  # Latent
  for (col in colnames(pg$latent)) {
    feat_df[[col]] <- NA_real_
    feat_df[[col]][idx] <- pg$latent[, col]
  }

  # TDA global (broadcast)
  for (col in colnames(pg$tda_global)) {
    feat_df[[col]] <- NA_real_
    feat_df[[col]][idx] <- pg$tda_global[[col]]
  }

  # TDA mapper (per-locus)
  for (col in colnames(pg$tda_mapper)) {
    feat_df[[col]] <- NA_real_
    feat_df[[col]][idx] <- pg$tda_mapper[[col]]
  }

  # Tree (per-locus)
  for (col in colnames(pg$tree)) {
    feat_df[[col]] <- NA_real_
    feat_df[[col]][idx] <- pg$tree[[col]]
  }
}

cat("  Structural features merged\n\n")

# --- Step 7a: Compute chromatin + functional annotation features ---
cat("Step 7a: Computing chromatin accessibility + functional annotation features...\n")

chrom_feat <- compute_chromatin_features(feat_df)
for (col in colnames(chrom_feat)) {
  feat_df[[col]] <- chrom_feat[[col]]
}
cat("  Chromatin features merged into feat_df\n\n")

# --- Step 7b: Compute pair + interaction features ---
cat("Step 7b: Computing pair + interaction features...\n")

# Identify chromatin columns (for interaction computation)
chrom_cols <- c(grep("^atac_", names(feat_df), value = TRUE),
                grep("^func_", names(feat_df), value = TRUE))

pair_features_all <- list()

for (grna_id in names(CONFIG$grnas)) {
  cat("\n  --- Pair + interaction features for gRNA", grna_id, "---\n")
  grna_cfg <- CONFIG$grnas[[grna_id]]
  mm_col <- grna_cfg$mismatch_col

  # Get indices for this guide
  if (grna_id %in% names(per_guide_struct)) {
    guide_idx <- per_guide_struct[[grna_id]]$indices
    tree_coph <- per_guide_struct[[grna_id]]$tree_cophenetic
  } else {
    next
  }

  guide_feat <- feat_df[guide_idx, , drop = FALSE]

  pf <- compute_all_pair_features(
    feat_df = guide_feat,
    feature_set = feature_set,
    grna_id = grna_id,
    mismatch_col = mm_col,
    tree_cophenetic = tree_coph
  )

  # Compute interaction features using this guide's mismatch profile + chromatin
  mm_profile_cols <- paste0("mm_pos_", 1:20)
  mm_profiles <- as.matrix(pf[, mm_profile_cols, drop = FALSE])
  chromatin_for_guide <- feat_df[guide_idx, chrom_cols, drop = FALSE]

  intx <- compute_interaction_features(
    mm_profiles = mm_profiles,
    mismatch_count = pf$mismatch_count,
    pam_mismatch = pf$PAM_mismatch,
    chromatin_df = chromatin_for_guide
  )

  # Combine pair features + interaction features
  pf <- cbind(pf, intx)

  # Store with guide index mapping
  pair_features_all[[grna_id]] <- list(indices = guide_idx, features = pf)
}

# Merge pair features into feat_df
for (grna_id in names(pair_features_all)) {
  pf <- pair_features_all[[grna_id]]
  idx <- pf$indices

  for (col in colnames(pf$features)) {
    if (!(col %in% c("mismatch_count", "hit_orientation"))) {
      # Only add new columns (don't overwrite existing mismatch_count etc.)
      if (!col %in% colnames(feat_df)) {
        feat_df[[col]] <- NA_real_
      }
      feat_df[[col]][idx] <- pf$features[[col]]
    }
  }
}

cat("\n  Pair + interaction features merged\n\n")

# --- Step 8: Assemble long-format training data ---
cat("Step 8: Assembling long-format training data...\n")

# Define feature groups for the long data
# We need to identify which columns belong to which group
# For now, collect all feature columns (excluding metadata and sample columns)

# Get all feature columns (everything that's not metadata or sample counts)
meta_cols <- c("chrom", "start", "end", "strand", "TOTAL_SUMS", "cluster_id",
               "clustered_sequence", "centered_sequence",
               "strand_collapsed_cluster_id", "strand_collapsed_cluster_total",
               "hit_orientation_1732", "best_start_1732", "best_end_1732",
               "best_substring_1732", "mismatch_count_1732",
               "hit_orientation_4894", "best_start_4894", "best_end_4894",
               "best_substring_4894", "mismatch_count_4894",
               "exp_flanking_sequence", "edits", "weight")
sample_cols_in_df <- sample_cols
exclude_cols <- c(meta_cols, sample_cols_in_df)

# Feature columns for the pair scorer
pair_scorer_feat_cols <- colnames(feat_df)[!colnames(feat_df) %in% exclude_cols]

# Build feature groups mapping
feature_groups_for_assembly <- list(
  sequence = intersect(feature_set, pair_scorer_feat_cols),
  alignment = intersect(c(grep("^mm_pos_", pair_scorer_feat_cols, value=TRUE),
                           "mismatch_count", "PAM_mismatch", "hit_orientation"),
                         pair_scorer_feat_cols),
  pair = intersect(c("pair_seq_identity", "pair_flanking_kmer_similarity",
                      "pair_tree_cophenetic", "pair_same_chromosome",
                      "pair_genomic_distance", "pair_gc_difference",
                      "pair_cpg_difference"), pair_scorer_feat_cols),
  per_guide_struct = intersect(c(grep("^H[01]_", pair_scorer_feat_cols, value=TRUE),
                                  grep("^land_H", pair_scorer_feat_cols, value=TRUE),
                                  grep("^mapper_", pair_scorer_feat_cols, value=TRUE),
                                  grep("^tree_", pair_scorer_feat_cols, value=TRUE),
                                  grep("^pg_PC", pair_scorer_feat_cols, value=TRUE)),
                                pair_scorer_feat_cols),
  global_struct = intersect(c(grep("^global_", pair_scorer_feat_cols, value=TRUE)),
                             pair_scorer_feat_cols),
  chromatin = intersect(c(grep("^atac_", pair_scorer_feat_cols, value=TRUE),
                           grep("^func_", pair_scorer_feat_cols, value=TRUE)),
                         pair_scorer_feat_cols),
  interaction = intersect(c(grep("^intx_", pair_scorer_feat_cols, value=TRUE)),
                           pair_scorer_feat_cols)
)

long_df <- assemble_long_data(feat_df, sample_cols, feature_groups_for_assembly)
cat("  Long-format data:", nrow(long_df), "rows\n")

# Define feature groups on the long data
groups <- define_feature_groups(long_df)
cat("  Feature groups:\n")
for (g in names(groups)) cat("    ", g, ":", length(groups[[g]]), "cols\n")
cat("\n")

# --- Step 9: Train/test split by locus ---
cat("Step 9: Splitting train/test by locus (80/20)...\n")
set.seed(CONFIG$seed)
unique_loci <- unique(long_df$locus_id)
test_loci <- sample(unique_loci, size = floor(length(unique_loci) * CONFIG$test_frac))
cat("  Train loci:", length(unique_loci) - length(test_loci), "\n")
cat("  Test loci:", length(test_loci), "\n\n")

# --- Step 9b: Export long-format data for Python neural model ---
cat("Step 9b: Exporting long-format data for neural model...\n")
nn_export_dir <- file.path(CONFIG$output_dir, "nn_data")
dir.create(nn_export_dir, showWarnings = FALSE, recursive = TRUE)

# Replace NA with 0 for export
long_df_export <- long_df
for (col in colnames(long_df_export)) {
  if (is.numeric(long_df_export[[col]])) {
    long_df_export[[col]][is.na(long_df_export[[col]])] <- 0
  }
}
write.csv(long_df_export, file.path(nn_export_dir, "long_data.csv"), row.names = FALSE)

# Export feature group mapping
group_mapping <- data.frame(
  group = rep(names(groups), sapply(groups, length)),
  feature = unlist(groups)
)
write.csv(group_mapping, file.path(nn_export_dir, "feature_groups.csv"), row.names = FALSE)

# Export test loci for reproducible split
write.csv(data.frame(locus_id = test_loci), file.path(nn_export_dir, "test_loci.csv"), row.names = FALSE)
cat("  Exported", nrow(long_df_export), "rows to", nn_export_dir, "\n\n")

# --- Step 10: Baseline comparison ---
cat("Step 10: Running baseline comparison...\n")
cat("  (sequence-only vs structural-only vs chromatin-only vs full vs full+chromatin vs full+chromatin+intx)\n\n")

baseline_results <- run_baseline_comparison(
  long_df = long_df,
  groups = groups,
  test_loci = test_loci,
  seed = CONFIG$seed
)

cat("\n========================================\n")
cat("RESULTS SUMMARY\n")
cat("========================================\n")
print(baseline_results$summary)

# Save results
write.csv(baseline_results$summary,
          file = file.path(CONFIG$output_dir, "baseline_comparison_metrics.csv"),
          row.names = FALSE)

# --- Step 11: Visualizations ---
cat("\nStep 11: Generating visualizations...\n")

# Predicted vs actual (full+chromatin+interaction model)
full_pred <- baseline_results$details$full_chromatin_intx$full_result$predictions
plot_pred_vs_actual(full_pred, "Full + chromatin + interaction model",
  out_path = file.path(CONFIG$output_dir, "pred_vs_actual_full.png"))

# Per-guide performance
plot_per_guide_performance(baseline_results$summary,
  out_path = file.path(CONFIG$output_dir, "per_guide_performance.png"))

# SHAP by group (full+chromatin+interaction model)
full_model <- baseline_results$details$full_chromatin_intx$full_result$model
full_cols <- c(groups$sequence, groups$alignment, groups$pair,
               groups$per_guide_struct, groups$global_struct,
               groups$chromatin, groups$interaction, groups$experimental)
plot_shap_by_group(full_model, long_df[!long_df$locus_id %in% test_loci, ],
                   full_cols, groups,
  out_path = file.path(CONFIG$output_dir, "shap_by_group.png"))

# Ablation chart
plot_ablation(baseline_results$summary,
  out_path = file.path(CONFIG$output_dir, "ablation.png"))

# Mapper graph for gRNA 1732
if ("1732" %in% names(per_guide_struct)) {
  pg <- per_guide_struct[["1732"]]
  grna_samples <- sample_cols[grepl("^1732", sample_cols)]
  guide_edits <- rowSums(feat_df[pg$indices, grna_samples, drop = FALSE])
  plot_mapper_graph(
    mapper_graph = pg$tda_full$mapper_graph,
    mapper_obj = pg$tda_full$mapper_obj,
    edits = guide_edits,
    title = "Mapper graph — gRNA 1732",
    out_path = file.path(CONFIG$output_dir, "mapper_graph_1732.png")
  )
}

cat("\n========================================\n")
cat("Pipeline complete.\n")
cat("Results saved to:", CONFIG$output_dir, "\n")
cat("========================================\n")
