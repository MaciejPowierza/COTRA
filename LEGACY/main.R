
# main.R — entry point (refactored from akd_script_to_inspect.txt)
# This script wires modules together. Edit the CONFIG section only.

#setwd("/home/kinga/UP_projekt_2_semestr/akd_refactor_modules_FEATURE_SELECTION/akd_refactor_modules")

# -------------------------
# CONFIG
# -------------------------
CONFIG <- list(
  data_path = "CROSSSTRAND_revcomp_forward_local_alignment2.csv", # <- change
  sep = ";",
  sample_col_start = "1732_1",
  sample_col_end   = "1732_4h_DNA6_1",
  choose_sample    = "1732_4h_DNA6_1",  # optional: run for a single sample
  seed = 123,
  # Dim reduction: "plain_pca", "weighted_pca", "nmf", "hpca", "rgcca", "sgcca", "bpca", "fa", "lpe", "ica"
  # Feature selection: "topvar", "lscore", "fosmod", "mcfs", "block_reps"
  dim_method = c("bpca","fa","ica"),
  selection_methods = c("topvar","lscore","fosmod","block_reps"),
  k_select_grid = c(20,50,100),
  k_dr_grid = c(2:16),
  fs_preprocess = "center",
  fs_graph_type = c("knn", 10),
  fs_lscore_t = 1,
  fs_mcfs_K = NULL,
  fs_mcfs_lambda = 1,
  fs_mcfs_t = 10,
  block_rep_per_block = 1,
  # Detector: "ocsvm", "isoforest", "mahalanobis"
  detector = "ocsvm",
  # Which rows are "normal" for one-class fitting: a function that returns indices
  normal_idx_fun = function(df) which(df$edits != 0),
  # Outlier thresholding (for iso/maha) when needed
  iso_quantile = 0.95,
  maha_alpha = 0.001,
  ocsvm_nu = 0.05,
  weight_col = "weight", #for wPCA
  nmf_rank = 3,
  nmf_nonneg_mode = "shift",     # or "drop_negative_features"
  nmf_feature_subset = "all",    # or "nonnegative_only"
  nmf_seed = 123,
  #experiment_dim_methods = c("rgcca"),
  #preproc_modes: "none", "basic", "aggressive"
  preproc_modes = c("none","basic","aggressive"),
  corr_cutoff = 0.95,
  #k_grid = c(2),
  rgcca_scheme = "factorial",
  rgcca_tau = "optimal",
  rgcca_sparsity = 0.8,
  rgcca_connection_mode = "full",
  rgcca_scores_use = "first",
  rgcca_scores_uses = c("first","first2","all"),
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
  block_rep_stat = "variance",
  metrics_output_path = paste0(
    "metrics_df_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".csv"
  ),
  selection_method = "none",
  k_select = NULL,
  k_dr = 16
)

# -------------------------
# Load modules
# -------------------------
source("R/00_setup.R")
source("R/01_io.R")
source("R/02_sampling.R")
source("R/03_features.R")
source("R/03b_preprocessing.R")
source("R/04a_feature_selection.R")
source("R/04_routing_blocks_hpca2.R")
source("R/04b_rgcca_blocks.R")
source("R/05_dim_reduction2.R")
source("R/06_detectors.R")
source("R/07_experiment2.R")
source("R/08_viz.R")

# -------------------------
# Run
# -------------------------
set.seed(CONFIG$seed)

bed_df <- read_bed_like(CONFIG$data_path, sep = CONFIG$sep)

# Option A: split into all samples found between columns
list_of_dfs <- split_into_samples(
  bed_df,
  col_start = CONFIG$sample_col_start,
  col_end   = CONFIG$sample_col_end
)

###===Beginning of insertion, 19.01===

results <- lapply(list_of_dfs, function(d) {
  feature_extraction(d, seq_col = "exp_flanking_sequence", draw_plot = FALSE)
})

###===End of insertion, 19.01===

# Option B: choose one sample only
one_df <- choose_one_sample(bed_df, CONFIG$choose_sample)

# Feature extraction on chosen dataset (use one_df or any list_of_dfs[[i]])
feat_df <- feature_extraction(list_of_dfs[[6]], seq_col = "exp_flanking_sequence", draw_plot = TRUE)

#print(colnames(feat_df))

# Run one analysis pass (dim reduction + detector + scoring)
run <- run_one_pass(
  feat_df,
  dim_method = CONFIG$dim_method[1],
  selection_method = CONFIG$selection_method,
  k_select = CONFIG$k_select,
  k_dr = CONFIG$k_dr,
  detector   = CONFIG$detector,
  normal_idx = CONFIG$normal_idx_fun(feat_df),
  weight_col = CONFIG$weight_col,
  iso_quantile = CONFIG$iso_quantile,
  maha_alpha   = CONFIG$maha_alpha,
  ocsvm_nu     = CONFIG$ocsvm_nu,
  #k = 10,
  nmf_nonneg_mode = CONFIG$nmf_nonneg_mode,
  nmf_feature_subset = CONFIG$nmf_feature_subset,
  nmf_seed = CONFIG$nmf_seed,
  preproc_mode = "basic",
  corr_cutoff = 0.95,
  rgcca_scheme = CONFIG$rgcca_scheme,
  rgcca_tau = CONFIG$rgcca_tau,
  rgcca_sparsity = CONFIG$rgcca_sparsity,
  rgcca_connection_mode = CONFIG$rgcca_connection_mode,
  rgcca_scores_use = CONFIG$rgcca_scores_use,
  rgcca_superblock = CONFIG$rgcca_superblock,
  lpe_preprocess = CONFIG$lpe_preprocess,
  lpe_numk = CONFIG$lpe_numk,
  ica_alg_typ = CONFIG$ica_alg_typ,
  ica_fun = CONFIG$ica_fun,
  ica_alpha = CONFIG$ica_alpha,
  ica_method = CONFIG$ica_method,
  ica_row_norm = CONFIG$ica_row_norm,
  ica_maxit = CONFIG$ica_maxit,
  ica_tol = CONFIG$ica_tol
)

print("RUN SUMMARY")
print(run$metrics_row)

# Visualization (PC1/PC2 + predictions + on-target)
plot_outlier_scatter(
  run$pc_scores,
  title = paste(CONFIG$dim_method, "+", CONFIG$detector)
)

#===Beginning of the insertion, 13.01===

#results_df <- run_experiment_summary_ontarget(
#  list_of_dfs,
#  edits_col = "edits",
#  svm_nu = 0.05,
#  seed = 123,
#  shuffle = TRUE,
#  ontarget_policy = "first"
#)
#print(results_df)

#===End of the insertion, 13.01===

metrics_df <- run_experiment_metrics(
  datasets = results,
  dim_methods = CONFIG$dim_method,
  detectors = c("ocsvm","isoforest"),
  edits_col = "edits",
  ontarget_policy = "first",
  positive_class = "offtarget",
  ocsvm_nu = 0.05,
  weight_col = CONFIG$weight_col,
  hpca_feature_mode = "super_blocks",
  preproc_modes = CONFIG$preproc_modes,
  corr_cutoff = CONFIG$corr_cutoff,
  #k_grid = CONFIG$k_grid,
  selection_methods = CONFIG$selection_methods,
  k_select_grid = CONFIG$k_select_grid,
  k_dr_grid = CONFIG$k_dr_grid,
  fs_preprocess = CONFIG$fs_preprocess,
  fs_graph_type = CONFIG$fs_graph_type,
  fs_lscore_t = CONFIG$fs_lscore_t,
  fs_mcfs_K = CONFIG$fs_mcfs_K,
  fs_mcfs_lambda = CONFIG$fs_mcfs_lambda,
  fs_mcfs_t = CONFIG$fs_mcfs_t,
  block_rep_per_block = CONFIG$block_rep_per_block,
  rgcca_scheme = CONFIG$rgcca_scheme,
  rgcca_tau = CONFIG$rgcca_tau,
  rgcca_sparsity = CONFIG$rgcca_sparsity,
  rgcca_connection_mode = CONFIG$rgcca_connection_mode,
  rgcca_scores_uses = CONFIG$rgcca_scores_uses,
  rgcca_superblock = CONFIG$rgcca_superblock,
  lpe_preprocess = CONFIG$lpe_preprocess,
  lpe_numk = CONFIG$lpe_numk,
  ica_alg_typ = CONFIG$ica_alg_typ,
  ica_fun = CONFIG$ica_fun,
  ica_alpha = CONFIG$ica_alpha,
  ica_method = CONFIG$ica_method,
  ica_row_norm = CONFIG$ica_row_norm,
  ica_maxit = CONFIG$ica_maxit,
  ica_tol = CONFIG$ica_tol
)

print(metrics_df)
write.csv(metrics_df, file = CONFIG$metrics_output_path, row.names = FALSE)
