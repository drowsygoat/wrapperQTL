#' Plot QTL Index vs. p-value/FDR (with optional -log10 transform)
#'
#' Creates a scatter plot of QTLs sorted by p-value (the "index"), showing
#' one or more metrics across the index: raw p-values, FDR, and optionally a
#' recalculated FDR column. The plot can be rendered on the natural scale or
#' as \eqn{-\log_{10}} values. A rasterized PDF is written using
#' \code{rasterpdf::raster_pdf()} to keep file sizes small while preserving
#' vector text.
#'
#' @param df A data frame containing at least the columns referenced by
#'   \code{p_col}, \code{fdr_col}, and optionally \code{fdr_recal_col}.
#' @param p_col Character scalar. Column name for p-values. Default \code{"p_value"}.
#' @param fdr_col Character scalar. Column name for FDR values. Default \code{"FDR"}.
#' @param fdr_recal_col Optional character scalar. Column name for a recalculated
#'   FDR metric to include in the plot. Default \code{NULL} (omit).
#' @param out_file Optional output filename for the PDF. If \code{NULL}, a name
#'   is derived from the object name of \code{df} as \code{"<df>_qtl_index.pdf"}.
#' @param width,height Numeric. PDF width and height (inches). Defaults \code{8}, \code{4}.
#' @param use_neglog10 Logical. If \code{TRUE}, metrics are plotted as
#'   \eqn{-\log_{10}}; otherwise, raw values are plotted. Default \code{FALSE}.
#'
#' @return The \code{ggplot2} plot object (invisibly saved to a rasterized PDF).
#' @export
plot_p_vals <- function(
  df,
  p_col = "p_value",
  fdr_col = "FDR",
  fdr_recal_col = NULL,
  out_file = NULL,
  width = 8,
  height = 4,
  use_neglog10 = FALSE
) {
  stopifnot(is.data.frame(df))

  p_sym         <- rlang::ensym(p_col)
  fdr_sym       <- rlang::ensym(fdr_col)
  fdr_recal_sym <- if (!is.null(fdr_recal_col)) rlang::ensym(fdr_recal_col) else NULL

  # basic presence checks
  if (!rlang::as_string(p_sym) %in% names(df))
    stop("Column '", rlang::as_string(p_sym), "' not found in df.")
  if (!rlang::as_string(fdr_sym) %in% names(df))
    stop("Column '", rlang::as_string(fdr_sym), "' not found in df.")
  if (!is.null(fdr_recal_sym) && !rlang::as_string(fdr_recal_sym) %in% names(df))
    stop("fdr_recal_col '", rlang::as_string(fdr_recal_sym), "' not found in df.")

  df_name <- deparse(substitute(df))

  if (is.null(out_file)) out_file <- paste0(df_name, "_qtl_index.pdf")

  # build plotting data, optionally including FDR_recal
  dd <- df %>%
    dplyr::arrange(!!p_sym) %>%
    dplyr::mutate(index = dplyr::row_number()) %>%
    {
      if (!is.null(fdr_recal_sym)) {
        dplyr::transmute(.,
          index,
          p_value   = !!p_sym,
          FDR       = !!fdr_sym,
          FDR_recal = !!fdr_recal_sym
        )
      } else {
        dplyr::transmute(.,
          index,
          p_value   = !!p_sym,
          FDR       = !!fdr_sym
        )
      }
    }

  # type checks
  if (!is.numeric(dd$p_value) || !is.numeric(dd$FDR)) {
    stop("Selected columns '", rlang::as_string(p_sym), "' and '",
         rlang::as_string(fdr_sym), "' must be numeric.")
  }
  if (!is.null(fdr_recal_sym) && (!"FDR_recal" %in% names(dd) || !is.numeric(dd$FDR_recal))) {
    stop("Selected fdr_recal_col '", rlang::as_string(fdr_recal_sym), "' must be numeric.")
  }

  # optional -log10 transform, guarding for FDR_recal presence
  if (use_neglog10) {
    eps <- .Machine$double.xmin
    dd <- dd %>%
      dplyr::mutate(
        p_value = -log10(pmax(p_value, eps)),
        FDR     = -log10(pmax(FDR,     eps))
      )
    if ("FDR_recal" %in% names(dd)) {
      dd <- dd %>% dplyr::mutate(FDR_recal = -log10(pmax(FDR_recal, eps)))
    }
    y_lab <- expression(-log[10](value))
  } else {
    y_lab <- "Value"
  }

  # pivot only the columns that exist
  cols_to_pivot <- intersect(c("p_value", "FDR", "FDR_recal"), names(dd))
  dd_long <- dd %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(cols_to_pivot),
      names_to = "metric",
      values_to = "value"
    )

  p <- ggplot2::ggplot(dd_long, ggplot2::aes(x = index, y = value, color = metric)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8) +
    ggplot2::labs(
      x = "QTLs (sorted by p-value)",
      y = y_lab,
      color = "Metric",
      title = paste0("QTL p-values and FDR: ", df_name)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      plot.title = ggplot2::element_text(size = 12, face = "bold")
    )

rasterpdf::raster_pdf(
    filename = out_file,
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  print(p)
  
rasterpdf::dev.off()

  # ggsave(
  #   filename = out_file,
  #   plot = p,
  #   device = "pdf",
  #   width = width,
  #   height = height,
  #   units = "in"
  # )
  message("Saved rasterized PDF via rasterpdf to: ", out_file)
p
}
