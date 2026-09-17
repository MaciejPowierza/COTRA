
# 03_features.R — feature extraction (344 features) from sequences
# Refactors feature_extraction from the original script. fileciteturn0file0

feature_extraction <- function(df, seq_col, draw_plot = TRUE, highlight_idx = NULL) {
  seqs <- Biostrings::DNAStringSet(df[[seq_col]])

  # GC content
  gc_mat <- Biostrings::letterFrequency(seqs, letters = c("G","C"), as.prob = TRUE)
  gc_content <- rowSums(gc_mat)

  # AT/GC skews
  at <- Biostrings::letterFrequency(seqs, letters = c("A","T"), as.prob = FALSE)
  gc <- Biostrings::letterFrequency(seqs, letters = c("G","C"), as.prob = FALSE)

  skew_AT <- (at[,1] - at[,2]) / pmax(1, (at[,1] + at[,2]))
  skew_GC <- (gc[,1] - gc[,2]) / pmax(1, (gc[,1] + gc[,2]))

  # CpG count
  cpg_cnt <- Biostrings::vcountPattern("CG", seqs)

  # k-mer frequencies (1..4), normalized per-row (bugfix vs original)
  kmer_freqs_list <- list()
  for (k in 1:4) {
    kmer_freq <- Biostrings::oligonucleotideFrequency(seqs, width = k, step = 1)
    rs <- rowSums(kmer_freq)
    kmer_freq <- kmer_freq / pmax(1, rs)
    kmer_freq <- as.data.frame(kmer_freq)
    colnames(kmer_freq) <- paste0("freq_", colnames(kmer_freq))
    kmer_freqs_list[[k]] <- kmer_freq
  }

  out <- cbind(
    df,
    gc_content = gc_content,
    skew_AT = skew_AT,
    skew_GC = skew_GC,
    cpg_cnt = cpg_cnt,
    kmer_freqs_list[[1]], kmer_freqs_list[[2]], kmer_freqs_list[[3]], kmer_freqs_list[[4]]
  )

  if (draw_plot) {
    plot(
      gc_mat[, "G"], gc_mat[, "C"],
      xlab = "G content", ylab = "C content",
      pch = 16, col = "lightgray",
      main = "GC scatterplot"
    )
    if (!is.null(highlight_idx)) {
      points(gc_mat[highlight_idx, "G"], gc_mat[highlight_idx, "C"], pch = 19, col = "red", cex = 1.2)
    }
  }

  out
}

#get_feature_set <- function(df_feat) {
#  start <- which(names(df_feat) == "gc_content")
#  end   <- which(names(df_feat) == "freq_TTTT")
#  if (length(start) == 0 || length(end) == 0) stop("Could not find gc_content or freq_TTTT in df.")
#  colnames(df_feat)[start:end]
#}

get_feature_set <- function(df_feat) {
  feature_candidates <- c(
    "gc_content", "skew_AT", "skew_GC", "cpg_cnt",
    grep("^freq_", names(df_feat), value = TRUE)
  )
  feature_candidates[feature_candidates %in% names(df_feat)]
}

extract_feature_matrix <- function(df_feat, feature_set = NULL) {
  if (is.null(feature_set)) feature_set <- get_feature_set(df_feat)
  X <- df_feat[, feature_set, drop = FALSE]
  X
}
