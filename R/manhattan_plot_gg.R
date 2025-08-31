#' Manhattan Plot with Optional Region Highlighting
#'
#' Creates a Manhattan-style plot from a data frame using \pkg{ggplot2}, with optional
#' highlighted genomic regions. Chromosome labels are naturally ordered (via
#' \code{gtools::mixedsort}), basepair positions are auto-scaled to Mb when values look
#' like basepairs, and a horizontal significance line can be added.
#'
#' If region columns are not present (or \code{show_regions = FALSE}), the plot is drawn
#' without any highlighted rectangles.
#'
#' @param df A data.frame containing chromosome, position, and p-value columns.
#' @param chr_col Column name (string) for chromosome identifiers. Default \code{"Chromosome"}.
#' @param pos_col Column name (string) for genomic position (bp or Mb). Default \code{"PeakPosition"}.
#' @param p_col Column name (string) for p-values. Default \code{"PeakPvalue"}.
#' @param region_start_col Column name (string) for region start (bp or Mb). Default \code{"RegionStart"}.
#' @param region_end_col Column name (string) for region end (bp or Mb). Default \code{"RegionEnd"}.
#' @param sig.level Numeric p-value threshold for a dashed horizontal line (e.g., \code{5e-8}).
#'   Use \code{NA} (default) to omit the line.
#' @param show_regions Logical; if \code{TRUE} (default) and both region columns are present
#'   in \code{df}, shaded rectangles are drawn over those regions.
#' @param region_alpha Alpha for region shading (0-1). Default \code{0.1}.
#' @param col Two-color vector used to alternate chromosome point colors. Default
#'   \code{c("gray40","gray60")}.
#' @param point_size Point size for scatter layer. Default \code{1}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' The function:
#' \enumerate{
#' \item Strips a leading \code{"chr"} from chromosome labels and orders chromosomes
#'       using natural sort.
#' \item Converts positions to megabases if any positions exceed 1e6.
#' \item Computes cumulative genomic positions to lay chromosomes end-to-end.
#' \item Alternates point colors by chromosome.
#' \item Optionally shades regions if \code{show_regions = TRUE} and both region columns
#'       exist in \code{df}. If region columns are missing, the plot is drawn without regions.
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   Chromosome = rep(paste0("chr", 1:3), each = 100),
#'   PeakPosition = rep(seq(1, 5e6, length.out = 100), 3),
#'   PeakPvalue = runif(300)
#' )
#' p <- manhattan_plot_gg(df, sig.level = 5e-8)
#' print(p)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#' @export
manhattan_plot_gg <- function(
  df,
  chr_col = "chr",
  pos_col = "location",
  p_col   = "FDR_recal",
  region_start_col = "RegionStart",
  region_end_col   = "RegionEnd",
  sig.level = NA,
  show_regions = TRUE,
  region_alpha = 0.1,
  col = c("gray40", "gray60"),
  point_size = 1
) {
  stopifnot(all(c(chr_col, pos_col, p_col) %in% names(df)))

  # Work on a copy
  data <- df

  # Normalize chr labels and order naturally
  data[[chr_col]] <- gsub("^chr", "", data[[chr_col]], ignore.case = TRUE)
  chr_levels <- unique(data[[chr_col]])
  chr_levels <- gtools::mixedsort(chr_levels)
  data[[chr_col]] <- factor(data[[chr_col]], levels = chr_levels, ordered = TRUE)

  # Convert bp -> Mb if values look like basepairs
  if (any(data[[pos_col]] > 1e6, na.rm = TRUE)) {
    data[[pos_col]] <- data[[pos_col]] / 1e6
    if (show_regions && all(c(region_start_col, region_end_col) %in% names(data))) {
      data[[region_start_col]] <- data[[region_start_col]] / 1e6
      data[[region_end_col]]   <- data[[region_end_col]]   / 1e6
    }
  }

  # Cumulative offsets per chromosome
  chr_offsets <- data %>%
    dplyr::group_by(!!rlang::sym(chr_col)) %>%
    dplyr::summarise(chr_len = max(!!rlang::sym(pos_col), na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(offset = dplyr::lag(cumsum(.data$chr_len), default = 0))

  data <- data %>%
    dplyr::left_join(chr_offsets, by = setNames(chr_col, chr_col)) %>%
    dplyr::mutate(
      pos_cum = !!rlang::sym(pos_col) + .data$offset,
      logp    = -log10(!!rlang::sym(p_col))
    )

  # Axis centers for chromosome labels
  axis_df <- data %>%
    dplyr::group_by(!!rlang::sym(chr_col)) %>%
    dplyr::summarise(center = mean(range(.data$pos_cum)), .groups = "drop")

  # Base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$pos_cum, y = .data$logp, color = !!rlang::sym(chr_col))) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = rep(col, length.out = nlevels(data[[chr_col]]))) +
    ggplot2::scale_x_continuous(breaks = axis_df$center, labels = axis_df[[chr_col]]) +
    ggplot2::labs(x = "Chromosome", y = expression(-log[10](p-value))) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )

  # Optional region shading (only if columns exist)
  if (show_regions && all(c(region_start_col, region_end_col) %in% names(data))) {
    regions <- df %>%
      dplyr::distinct(!!rlang::sym(chr_col), !!rlang::sym(region_start_col), !!rlang::sym(region_end_col)) %>%
      dplyr::mutate(
        !!rlang::sym(chr_col) := gsub("^chr", "", !!rlang::sym(chr_col), ignore.case = TRUE),
        !!rlang::sym(chr_col) := factor(!!rlang::sym(chr_col), levels = chr_levels, ordered = TRUE)
      ) %>%
      dplyr::left_join(chr_offsets, by = setNames(chr_col, chr_col)) %>%
      dplyr::mutate(
        start_cum = !!rlang::sym(region_start_col) / ifelse(any(df[[pos_col]] > 1e6), 1e6, 1) + .data$offset,
        end_cum   = !!rlang::sym(region_end_col)   / ifelse(any(df[[pos_col]] > 1e6), 1e6, 1) + .data$offset
      )

    p <- p +
      ggplot2::geom_rect(
        data = regions,
        ggplot2::aes(xmin = .data$start_cum, xmax = .data$end_cum, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "blue",
        alpha = region_alpha
      )
  }

  # Optional significance line
  if (!is.na(sig.level)) {
    p <- p + ggplot2::geom_hline(yintercept = -log10(sig.level), linetype = "dashed", color = "red")
  }

  p
}
