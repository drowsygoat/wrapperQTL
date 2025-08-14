#' Read Chromosome and Position from a VCF File
#'
#' Extracts chromosome and position columns from a VCF file, excluding header lines.
#'
#' @param vcf_path Path to the VCF file (can be gzipped or plain text).
#'
#' @return A tibble with two columns:
#' \itemize{
#'   \item \code{chr_snp}: Chromosome identifiers, optionally converted using \code{convert_chromosome_ids()}.
#'   \item \code{location_snp}: SNP positions as character strings.
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Uses \code{grep -v '^##'} to skip metadata header lines.
#'   \item Reads only the first two columns from the VCF: chromosome and position.
#'   \item Renames them to \code{chr_snp} and \code{location_snp}.
#'   \item Applies \code{convert_chromosome_ids()} function to standardize chromosome names.
#' }
#'
#' @note This function assumes that a helper function \code{convert_chromosome_ids()} is available in the environment.
#'
#' @import data.table
#' @import tibble
#'
#' @examples
#' \dontrun{
#' vcf_file <- "variants.vcf.gz"
#' chr_pos_tbl <- read_chr_pos_from_vcf(vcf_file)
#' }
#'
#' @export
read_chr_pos_from_vcf <- function(vcf_path) {
  library(data.table)
  library(tibble)

  # Read only lines that are not headers and select the first two columns
  dt <- fread(cmd = paste("grep -v '^##' ", shQuote(vcf_path)), select = 1:2)

  # Rename columns
  setnames(dt, c("chr_snp", "location_snp"))

  dt$chr_snp <- convert_chromosome_ids(dt$chr_snp)
  dt$location_snp <- as.character(dt$location_snp)

  return(as_tibble(dt))
}