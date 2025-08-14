#' Plot PeakMatrix Boxplots per Sample
#'
#' Generates a per-sample boxplot of peak signals from a `SummarizedExperiment` object. 
#' Optionally subsets to the top N peaks per sample before plotting. 
#' Can save the plot to a PDF and optionally return the subsetted object.
#'
#' @param se A `SummarizedExperiment` containing the peak signal assay.
#' @param assay_name Name of the assay to use. Default is `"PeakMatrix"`.
#' @param top_n_peaks Optional integer. If set, selects the top N peaks 
#'   (by signal) per sample before plotting. The function will return 
#'   the subsetted `SummarizedExperiment` invisibly in this case.
#' @param title Optional plot title. If `NULL`, a default title is generated.
#' @param outfile Optional file path. If provided, saves the plot as a high-resolution PDF.
#' @param lim_y Numeric; upper y-axis limit for plotting. Default is `1`.
#'
#' @return 
#' - If `top_n_peaks` is provided: the subsetted `SummarizedExperiment` (invisibly).
#' - If `top_n_peaks` is `NULL`: `NULL` (invisibly).
#'
#' @details
#' The function:
#' 1. Validates that the input is a `SummarizedExperiment` with the specified assay.
#' 2. Optionally subsets to the top N peaks per sample.
#' 3. Converts the peak matrix to a long format (supports dense and sparse matrices).
#' 4. Creates a combined boxplot + quasi-random scatter plot of peak signals.
#'
#' Sparse matrices (`dgCMatrix` or `dgTMatrix`) are converted internally for plotting.
#' Missing values are ignored by `ggplot2`.
#'
#' @examples
#' \dontrun{
#' library(SummarizedExperiment)
#'
#' # Example SE object with assay "PeakMatrix"
#' plotPeakMatrixBoxplot(se, top_n_peaks = 1000, outfile = "peaks.pdf")
#' }
#'
#' @export
plotPeakMatrixBoxplot <- function(
  se,
  assay_name   = "PeakMatrix",
  top_n_peaks  = NULL,
  title        = NULL,
  outfile      = NULL,
  lim_y        = 1
) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("Input must be a SummarizedExperiment object.")
  }
  if (!(assay_name %in% SummarizedExperiment::assayNames(se))) {
    stop(paste("Assay", assay_name, "not found."))
  }

  peak_matrix <- SummarizedExperiment::assay(se, assay_name)

  # Optionally subset to top N peaks per sample
  subsetted <- FALSE
  if (!is.null(top_n_peaks)) {
    message("Selecting top ", top_n_peaks, " peaks per sample...")

    if (inherits(peak_matrix, "dgCMatrix") || inherits(peak_matrix, "dgTMatrix")) {
      peak_matrix_dense <- as.matrix(peak_matrix)
    } else {
      peak_matrix_dense <- peak_matrix
    }

    top_peaks <- unique(as.integer(unlist(
      apply(peak_matrix_dense, 2, function(x) {
        ord <- order(x, decreasing = TRUE)
        ord[seq_len(min(top_n_peaks, length(ord)))]
      })
    )))

    se <- se[as.vector(top_peaks), ]
    subsetted <- TRUE
    peak_matrix <- SummarizedExperiment::assay(se, assay_name)
  }

  # Long format conversion
  if (inherits(peak_matrix, "dgCMatrix") || inherits(peak_matrix, "dgTMatrix")) {
    peak_matrix <- methods::as(peak_matrix, "dgTMatrix")
    long_df <- data.table::data.table(
      Sample     = colnames(se)[peak_matrix@j + 1],
      PeakSignal = peak_matrix@x
    )
  } else {
    df <- as.data.frame(as.matrix(peak_matrix))
    data.table::setDT(df)
    df[, Row := .I]
    long_df <- data.table::melt(df, id.vars = "Row",
                                variable.name = "Sample",
                                value.name = "PeakSignal")
  }

  # Title handling
  if (is.null(title)) {
    title <- paste0(
      "PeakMatrix boxplot",
      if (!is.null(top_n_peaks)) paste0(" (top ", top_n_peaks, " peaks/sample)") else ""
    )
  }

  # Plot
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = .data$Sample, y = .data$PeakSignal)) +
    ggbeeswarm::geom_quasirandom(alpha = 0.2, size = 0.1, color = "grey60") +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = NA, color = "black", alpha = 0.2, width = 0.5) +
    ggplot2::coord_cartesian(ylim = c(0, lim_y)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggplot2::labs(title = title, x = "Sample", y = "Peak signal")

  if (!is.null(outfile)) {
    rasterpdf::raster_pdf(outfile, width = 12, height = 6, res = 300)
    print(p)
    grDevices::dev.off()
  } else {
    print(p)
  }

  if (subsetted) {
    return(invisible(se))
  } else {
    return(invisible(NULL))
  }
}
