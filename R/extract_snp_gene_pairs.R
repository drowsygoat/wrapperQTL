
#' Extract SNP-Gene Pairs with Chromosome and Position
#'
#' Extracts unique SNP-gene pairs from a data frame and parses SNP identifiers into chromosome and position columns.
#'
#' @param df A data frame containing SNP and gene identifiers.
#' @param snp_col Name of the column in `df` containing SNP identifiers in the format `"CHR:POS"` or `"CHR:POS_ALT"`. Default is `"snp_id"`.
#' @param gene_col Name of the column in `df` containing gene or peak identifiers. Default is `"peak_id"`.
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item \code{SNP_ID}: The SNP identifier.
#'   \item \code{GENE_ID}: The gene or peak identifier.
#'   \item \code{CHR}: Chromosome parsed from \code{SNP_ID}.
#'   \item \code{POS}: Position parsed from \code{SNP_ID} as an integer.
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Renames the input columns to standardized names (\code{SNP_ID}, \code{GENE_ID}).
#'   \item Ensures uniqueness of SNP-gene pairs.
#'   \item Extracts chromosome and position information from the SNP ID.
#' }
#'
#' @import dplyr
#'
#' @examples
#' df <- data.frame(
#'   snp_id = c("1:12345_A", "2:67890_T"),
#'   peak_id = c("geneA", "geneB")
#' )
#' extract_snp_gene_pairs(df)
#'
#' @export
extract_snp_gene_pairs <- function(df, snp_col = "snp_id", gene_col = "peak_id") {
  df %>%
    dplyr::select(SNP_ID = all_of(snp_col), GENE_ID = all_of(gene_col)) %>%
    distinct() %>%
    mutate(
      CHR = sub(":.*", "", SNP_ID),
      POS = as.integer(sub(".*:(\\d+).*", "\\1", SNP_ID))
    )
}