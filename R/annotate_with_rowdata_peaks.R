#' Annotate a tibble with peak coordinates from a SummarizedExperiment
#'
#' Left-joins \code{df} to \code{rowData(se)} by matching \code{df[[key_col]]}
#' to \code{rownames(se)}. Appends four columns to \code{df}:
#' \code{peak_id} (the matched key), \code{peak_chr} (from \code{seqnames}),
#' \code{peak_start} (from \code{start}), and \code{peak_end} (from \code{end}).
#' The number and order of rows in \code{df} are preserved.
#'
#' @param df A data.frame or tibble to annotate (target table).
#' @param se A \code{SummarizedExperiment} whose \code{rowData} contains
#'   \code{seqnames}, \code{start}, and \code{end}; \code{rownames(se)} must be
#'   the peak identifiers to match.
#' @param key_col Character scalar. Column in \code{df} that contains peak IDs
#'   matching \code{rownames(se)}. Default \code{"peak_id"}.
#' @param seq_field,start_field,end_field Character scalars naming the columns
#'   in \code{rowData(se)} to use. Defaults: \code{"seqnames"}, \code{"start"}, \code{"end"}.
#' @param warn_unmatched Logical; warn if some keys in \code{df} do not match any
#'   \code{rownames(se)}. Default \code{TRUE}.
#'
#' @return \code{df} with added columns: \code{peak_id}, \code{peak_chr},
#'   \code{peak_start}, \code{peak_end}.
#'
#' @examples
#' \dontrun{
#' # df has "peak_ID" that matches rownames(se_atac_full)
#' df2 <- annotate_with_rowdata_peaks(
#'   df = df,
#'   se = se_atac_full,
#'   key_col = "peak_ID"  # becomes df2$peak_id
#' )
#' }
#'
#' @importFrom SummarizedExperiment rowData
#' @export
annotate_with_rowdata_peaks <- function(
  df,
  se,
  key_col = "peak_id",
  seq_field = "seqnames",
  start_field = "start",
  end_field = "end",
  warn_unmatched = TRUE
) {
  # ---- checks ----
  if (!is.data.frame(df)) stop("df must be a data.frame/tibble.")
  if (!methods::is(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment.")
  if (!is.character(key_col) || length(key_col) != 1L) stop("key_col must be a single column name.")
  if (!key_col %in% names(df)) stop(sprintf("Column '%s' not found in df.", key_col))

  rn <- rownames(se)
  if (is.null(rn)) stop("rownames(se) are NULL; set rownames(se) to peak IDs before joining.")

  rd <- SummarizedExperiment::rowData(se)
  needed <- c(seq_field, start_field, end_field)
  miss <- setdiff(needed, colnames(rd))
  if (length(miss)) {
    stop(sprintf("Missing columns in rowData(se): %s", paste(miss, collapse = ", ")))
  }

  # ---- left join via match (preserves df's order/rows) ----
  key_vec <- df[[key_col]]
  idx <- match(key_vec, rn)  # NA where no match

  if (warn_unmatched && anyNA(idx)) {
    n_un <- sum(is.na(idx))
    warning(sprintf("annotate_with_rowdata_peaks: %d/%d keys in df[%s] did not match rownames(se).",
                    n_un, length(idx), key_col))
  }

  # Helper to coerce S4Vectors::Rle or other vector-like to base vector
  to_vec <- function(x) if (inherits(x, "Rle")) as.vector(x) else as.vector(x)

  # Pull & map fields
  chr_vec   <- to_vec(rd[[seq_field]])[idx]
  start_vec <- to_vec(rd[[start_field]])[idx]
  end_vec   <- to_vec(rd[[end_field]])[idx]

  # ---- append output columns (overwrite if they already exist) ----
  df[["peak_id"]]    <- as.character(key_vec)
  df[["peak_chr"]]   <- as.character(chr_vec)
  df[["peak_start"]] <- as.integer(start_vec)
  df[["peak_end"]]   <- as.integer(end_vec)

  as_tibble(df)
}
