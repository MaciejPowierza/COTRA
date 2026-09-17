
# 02_sampling.R — selecting samples / splitting into per-sample data.frames

#===Beginning of the insertion, 17.02===

choose_one_sample <- function(df, sample_to_choose, weight_name = "weight") {
  needed <- c("contig","coord_start","coord_end","strand", sample_to_choose,
              "exp_flanking_sequence","edits")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  df_out <- df[, needed]
  df_out <- df_out[df_out[[sample_to_choose]] != 0, , drop = FALSE]

  names(df_out)[names(df_out) == sample_to_choose] <- weight_name

  df_out
}

#===End of the insertion, 17.02===

choose_one_sample2 <- function(df, sample_to_choose) {
  needed <- c("contig","coord_start","coord_end","strand", sample_to_choose, "exp_flanking_sequence","edits")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  df_out <- df[, needed]
  df_out <- df_out[df_out[[sample_to_choose]] != 0, , drop = FALSE]
  df_out
}

split_into_samples <- function(df, col_start, col_end) {
  start <- which(names(df) == col_start)
  end   <- which(names(df) == col_end)
  if (length(start) == 0 || length(end) == 0) stop("Could not find col_start/col_end in data.")
  if (start > end) stop("col_start must appear before col_end.")
  cols_set <- colnames(df)[start:end]

  out_list <- setNames(vector("list", length(cols_set)), cols_set)
  for (nm in cols_set) out_list[[nm]] <- choose_one_sample(df, nm)
  out_list
}

get_on_targets <- function(list_of_dfs) {
  lapply(list_of_dfs, function(d) d[d$edits == 0, , drop = FALSE])
}
