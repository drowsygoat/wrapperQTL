
#' Extract Genotypes from VCF by Chromosome and Position
#'
#' Extracts genotype data for specific SNP positions from a VCF file and formats genotypes as unphased allele pairs (e.g., `A/G`, `C/C`).
#'
#' @param vcf_path Path to the VCF file (can be gzipped or plain text).
#' @param snp_df A data frame containing SNP positions to extract. Must include chromosome and position columns.
#' @param chr_col Name of the column in `snp_df` that contains chromosome identifiers. Default is `"CHR"`.
#' @param pos_col Name of the column in `snp_df` that contains position values. Default is `"POS"`.
#' @param Capitalize Logical; whether to capitalize cleaned sample names. Default is `TRUE`.
#' @param pattern Optional integer vector specifying which parts of sample names (split by "_") to retain. Used for sample name cleaning.
#'
#' @return A data frame of genotypes with SNP IDs as rows and cleaned sample names as columns. Genotypes are in unphased format (e.g., `A/T`, `G/G`).
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts VCF lines corresponding to provided chromosome-position pairs.
#'   \item Parses genotype fields (`GT`) and maps 0/1 alleles to REF/ALT bases.
#'   \item Cleans sample names using the provided `pattern` or automatic heuristics.
#' }
#'
#' SNPs are identified by combining chromosome, position, and ALT allele (`CHR:POS_ALT`). Lines without valid genotypes are returned as `NA`.
#'
#' @import data.table
#' @import dplyr
#' @import tibble
#'
#' @examples
#' \dontrun{
#' snp_df <- data.frame(CHR = c("1", "2"), POS = c("12345", "67890"))
#' vcf_file <- "variants.vcf.gz"
#' genotypes <- extract_genotypes_from_vcf_by_pos(vcf_file, snp_df)
#' }
#'
#' @export

extract_genotypes_from_vcf_by_pos <- function(
  vcf_path,
  snp_df,
  chr_col = "CHR",
  pos_col = "POS",
  Capitalize = TRUE,
  pattern = NULL
) {
  library(data.table)
  library(dplyr)
  library(tibble)

  # Helper: Clean sample names
  clean_ids <- function(ids, Capitalize = TRUE, pattern = NULL) {
    ids <- as.character(ids)
    parts <- strsplit(ids, "_")
    cleaned <- vapply(parts, function(x) {
      if (!is.null(pattern) && all(pattern <= length(x))) {
        x_keep <- x[pattern]
      } else {
        n <- length(x)
        if (n %% 2 == 0 && all(x[1:(n / 2)] == x[(n / 2 + 1):n])) {
          x_keep <- x[1:(n / 2)]
        } else {
          x_keep <- x
        }
      }
      paste0(x_keep, collapse = "")
    }, character(1))
    if (Capitalize) toupper(cleaned) else tolower(cleaned)
  }

  if (!file.exists(vcf_path)) stop("VCF file does not exist: ", vcf_path)

  # Step 1: Create grep pattern from CHR + POS
  snp_df_clean <- snp_df %>%
    mutate(
      CHR_CLEAN = trimws(as.character(.data[[chr_col]])),
      POS_CLEAN = trimws(as.character(.data[[pos_col]])),
      grep_pattern = paste0("^", CHR_CLEAN, "\t", POS_CLEAN, "\t")
    )

  grep_pattern <- paste(unique(snp_df_clean$grep_pattern), collapse = "|")

  # Step 2: Extract VCF header
  header_lines <- system(paste0("grep '^#' ", shQuote(vcf_path)), intern = TRUE)
  header_line <- grep("^#CHROM", header_lines, value = TRUE)
  header <- strsplit(sub("^#", "", header_line), "\t")[[1]]
  if (length(header) < 10) stop("Could not parse VCF header correctly.")

  # Step 3: Extract matching data lines
  variant_lines <- system(paste0("grep -P '", grep_pattern, "' ", shQuote(vcf_path)), intern = TRUE)
  if (length(variant_lines) == 0) stop("No matching variants found in the VCF file for provided positions.")

  # Step 4: Read matching VCF lines
  vcf_dt <- fread(text = paste(c(paste(header, collapse = "\t"), variant_lines), collapse = "\n"))
  colnames(vcf_dt)[1] <- "CHROM"  # fix for "#CHROM"

  # Step 5: Build SNP_IDs
  vcf_dt$SNP_ID <- paste0(vcf_dt$CHROM, ":", vcf_dt$POS, "_", vcf_dt$ALT)
  ref <- vcf_dt$REF
  alt <- vcf_dt$ALT
  rownames_tag <- vcf_dt$SNP_ID

  # Step 6: Extract sample columns and clean names
  fixed_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
  sample_cols <- setdiff(names(vcf_dt), fixed_cols)
  cleaned_sample_cols <- clean_ids(sample_cols, Capitalize = Capitalize, pattern = pattern)

  # Step 7: Convert GT to unphased allele format
  gt_index <- 1
  mat <- matrix(NA_character_, nrow = nrow(vcf_dt), ncol = length(sample_cols))
  rownames(mat) <- rownames_tag
  colnames(mat) <- cleaned_sample_cols

  for (i in seq_along(sample_cols)) {
    raw_vals <- vcf_dt[[sample_cols[i]]]
    split_gt <- tstrsplit(raw_vals, ":", fixed = TRUE)
    gt_vals <- split_gt[[gt_index]]

    mat[, i] <- mapply(function(gt, r, a) {
      if (is.na(gt) || gt %in% c("./.", ".|.", ".")) return(NA_character_)
      alleles <- unlist(strsplit(gt, "[/|]"))
      if (length(alleles) != 2) return(NA_character_)
      bases <- sapply(alleles, function(x) {
        if (x == "0") r else if (x == "1") a else NA_character_
      })
      if (any(is.na(bases))) return(NA_character_)
      paste(sort(as.character(bases)), collapse = "/")
    }, gt_vals, ref, alt)
  }

  # Final genotype data frame
  genotype_df <- as.data.frame(mat)

  # Ensure no duplicate SNP_ID column
  if ("SNP_ID" %in% colnames(genotype_df)) {
    genotype_df <- genotype_df[, colnames(genotype_df) != "SNP_ID"]
  }

  genotype_df <- genotype_df %>%
    rownames_to_column("SNP_ID")

  return(genotype_df)
}
