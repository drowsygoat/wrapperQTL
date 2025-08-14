#' Plot Cluster–Sample Summary (Mean/Median of Non‑Zero Values)
#'
#' Reads a serialized \code{SummarizedExperiment} (\code{.rds}), computes per
#' \strong{cluster × sample} summaries of a chosen assay (mean and median of
#' non‑zero entries), and saves a combined boxplot with jittered points.
#'
#' @param se_path Character path to an \code{.rds} file containing a
#'   \code{SummarizedExperiment}.
#' @param save_dir Directory where the plot file will be written
#'   (created if it does not exist).
#' @param assay_name Optional character name of the assay to use. If \code{NULL}
#'   (default), the first assay in \code{assayNames(se)} is used.
#' @param filename Output filename for the saved plot (PDF/PNG depending on
#'   extension). Default: \code{"cluster_summary_plot.pdf"}.
#'
#' @return Invisibly returns \code{NULL}. The function is used for its side
#'   effects (writing a plot).
#'
#' @details
#' The input \code{SummarizedExperiment} must contain at least two columns in
#' \code{colData(se)}:
#' \itemize{
#'   \item \code{cluster}: cluster identifier (coerced to character)
#'   \item \code{sample_id}: sample identifier (coerced to character)
#' }
#' For each cluster, samples present in that cluster are considered. For each
#' \code{cluster × sample} subset, all \emph{non‑zero} entries of the assay are
#' taken, and their mean and median are computed. These statistics are reshaped
#' to long format and visualized as boxplots with quasi‑random points on top.
#'
#' If \code{assay_name} is \code{NULL}, the first assay is selected and a note
#' is messaged. The plot is saved with \code{ggplot2::ggsave()} to
#' \code{file.path(save_dir, filename)}.
#'
#' @section Package Requirements:
#' This function uses \pkg{SummarizedExperiment}, \pkg{Matrix}, \pkg{ggplot2},
#' \pkg{ggbeeswarm}, \pkg{dplyr}, and \pkg{tidyr}. The function checks for these
#' packages at runtime.
#'
#' @examples
#' \dontrun{
#' PlotClusterCovariateSummary(
#'   se_path   = "path/to/object.rds",
#'   save_dir  = "results/plots",
#'   assay_name = NULL,
#'   filename   = "cluster_summary_plot.pdf"
#' )
#' }
#'
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_minimal labs theme element_text ggsave
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom dplyr bind_rows mutate
#' @importFrom tidyr pivot_longer
#' @export
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
  plot <- ggplot2::ggplot(summary_df, ggplot2::aes(x = cluster, y = expression, fill = statistic)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggbeeswarm::geom_quasirandom(
      dodge.width = 0.8,
      size = 0.2,
      alpha = 0.8,
      color = "black"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      x = "Cluster", y = "Expression",
      title = "Mean & Median Expression (non-zero) per Cluster-Sample"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Save
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(save_dir, filename), plot, width = 20, height = 12, limitsize = FALSE)
  message(sprintf("✅ Plot saved to: %s", file.path(save_dir, filename)))
}
