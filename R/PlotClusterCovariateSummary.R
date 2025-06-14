PlotClusterCovariateSummary <- function(se_path, save_dir, assay_name = NULL, filename = "cluster_summary_plot.pdf") {
  # Load packages
  load_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("❌ Package '%s' is required but not installed.", pkg))
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  pkgs <- c("SummarizedExperiment", "Matrix", "ggplot2", "ggbeeswarm", "dplyr", "tidyr")
  invisible(lapply(pkgs, load_pkg))

  if (!file.exists(se_path)) {
    stop("❌ SummarizedExperiment file not found at: ", se_path)
  }

  se <- readRDS(se_path)
  if (is.null(assay_name)) {
    assay_name <- SummarizedExperiment::assayNames(se)[1]
    message("ℹ️ Using default assay: ", assay_name)
  }

  mat <- SummarizedExperiment::assay(se, assay_name)
  col_data <- as.data.frame(SummarizedExperiment::colData(se))
  col_data$cluster <- as.character(col_data$cluster)
  col_data$sample_id <- as.character(col_data$sample_id)

  all_clusters <- sort(unique(col_data$cluster))
  all_samples <- sort(unique(col_data$sample_id))

  result_list <- list()

  for (clust in all_clusters) {
    cols_in_cluster <- which(col_data$cluster == clust)
    sub_mat <- mat[, cols_in_cluster, drop = FALSE]
    sub_coldata <- col_data[cols_in_cluster, , drop = FALSE]
    sample_groups <- split(seq_along(cols_in_cluster), sub_coldata$sample_id)

    for (sid in all_samples) {
      if (!sid %in% names(sample_groups)) next

      vals <- as.numeric(sub_mat[, sample_groups[[sid]], drop = FALSE])
      vals <- vals[vals != 0]

      if (length(vals) == 0) next

      result_list[[length(result_list) + 1]] <- data.frame(
        cluster = clust,
        sample_id = sid,
        mean = mean(vals, na.rm = TRUE),
        median = median(vals, na.rm = TRUE)
      )
    }
  }

  summary_df <- dplyr::bind_rows(result_list) %>%
    tidyr::pivot_longer(cols = c("mean", "median"), names_to = "statistic", values_to = "expression") %>%
    dplyr::mutate(statistic = factor(statistic, levels = c("mean", "median")))
    print(summary_df)

  # Plot
  plot <- ggplot(summary_df, aes(x = cluster, y = expression, fill = statistic)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggbeeswarm::geom_quasirandom(
      dodge.width = 0.8,
      size = 0.2,
      alpha = 0.8,
      color = "black"
    ) +
    # facet_wrap(~statistic, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(
      x = "Cluster", y = "Expression",
      title = "Mean & Median Expression (non-zero) per Cluster-Sample"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Save
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(save_dir, filename), plot, width = 20, height = 12, limitsize = FALSE)
  message(sprintf("✅ Plot saved to: %s", file.path(save_dir, filename)))
}
