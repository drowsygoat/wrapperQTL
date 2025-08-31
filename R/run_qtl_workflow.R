#' Run the QTL post-processing workflow (+ optional nearest-gene annotation)
#'
#' @param results QTL results data.frame/tibble with columns:
#'   \code{p_value}, \code{chr_snp}, \code{location_snp}, \code{peak_id}, \code{QTL_type}
#'   (and optionally \code{FDR}).
#' @param vcf_file Path to VCF used to filter SNPs (via \code{read_chr_pos_from_vcf}).
#' @param apply_peak_range Logical; run \code{peak_range()} per cis/trans and join back. Default TRUE.
#' @param peak_params List of params passed to \code{peak_range()}.
#'   Default: \code{list(METHOD="DistanceFromPeak", SIGTHRESHOLD=0.01, SEARCHDISTANCE=3e6, LOCALPEAKTHRESHOLD=0.01, MINGAPSIZE=0)}.
#' @param recalc_fdr Logical; if TRUE, add \code{FDR_recal := p.adjust(p_value,"BH")}. Default TRUE.
#' @param annotate_peaks Logical; if TRUE, call \code{annotate_with_rowdata_peaks()} on cis/trans. Default FALSE.
#' @param se SummarizedExperiment for \code{annotate_with_rowdata_peaks()} (required if \code{annotate_peaks=TRUE}).
#' @param key_col Peak ID column in \code{results} for annotation. Default "peak_id".
#' @param nearest_genes Logical; if TRUE, annotate nearest gene by TSS from a GTF. Default FALSE.
#' @param gtf_file Path to GTF (plain or .gz). Required if \code{nearest_genes=TRUE} and \code{tss_gr} is NULL.
#' @param nearest_on Compute nearest gene relative to \code{"peak"} (columns \code{chr, location}) or \code{"snp"} (\code{chr_snp, location_snp}). Default "peak".
#' @param nearest_output What to add: \code{"gene"}, \code{"distance"}, or \code{"both"}. Default "both".
#' @param tss_gr Optional precomputed TSS GRanges (from \code{nearest_gene_from_gtf_build_tss()}) to speed up repeated calls.
#' @param output_dir Directory for outputs when \code{write_outputs=TRUE}. Default ".".
#' @param prefix Filename prefix for outputs. Default "qtl".
#' @param write_outputs Logical; write CSVs and PDFs. Default FALSE.
#'
#' @return List:
#' \itemize{
#'   \item \code{filtered} – VCF-filtered results.
#'   \item \code{peak_counts} – counts per QTL_type.
#'   \item \code{cis_ranges}, \code{trans_ranges} – outputs from \code{peak_range()} (if used).
#'   \item \code{cis}, \code{trans} – final tables (with optional FDR recalculation, rowData & nearest gene annotations).
#'   \item \code{plot_files} – written plot paths (if any).
#'   \item \code{tss} – the TSS GRanges used (if nearest_genes=TRUE).
#' }
#' @export
run_qtl_workflow <- function(
  results,
  vcf_file,
  apply_peak_range = TRUE,
  peak_params = list(
    METHOD = "DistanceFromPeak",
    SIGTHRESHOLD = 0.01,
    SEARCHDISTANCE = 3e6,
    LOCALPEAKTHRESHOLD = 0.01,
    MINGAPSIZE = 0
  ),
  recalc_fdr = TRUE,
  annotate_peaks = FALSE,
  se = NULL,
  key_col = "peak_id",
  nearest_genes = FALSE,
  gtf_file = NULL,
  nearest_on = c("peak","snp"),
  nearest_output = c("both","gene","distance"),
  tss_gr = NULL,
  output_dir = ".",
  prefix = "qtl",
  write_outputs = FALSE
) {
  nearest_on     <- match.arg(nearest_on)
  nearest_output <- match.arg(nearest_output)

  stopifnot(is.data.frame(results), is.character(vcf_file), file.exists(vcf_file))
  req_cols <- c("p_value", "chr_snp", "location_snp", "peak_id", "QTL_type")
  missing <- setdiff(req_cols, names(results))
  if (length(missing)) stop("results missing required columns: ", paste(missing, collapse = ", "))

  if (annotate_peaks && is.null(se)) stop("annotate_peaks=TRUE requires 'se'.")

  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(tibble)
  })

  # --- preprocess ---
  res0 <- results %>%
    mutate(
      minus_log10_p   = -log10(pmax(p_value, .Machine$double.xmin)),
      minus_log10_FDR = if ("FDR" %in% names(.)) -log10(pmax(FDR, .Machine$double.xmin)) else NA_real_,
      chr      = as.character(chr_snp),
      location = as.integer(location_snp)
    ) %>%
    as_tibble()

  # --- VCF filter ---
  vcf_pos  <- read_chr_pos_from_vcf(vcf_file)               # must return chr_snp, location_snp
  filtered <- inner_join(res0, vcf_pos, by = c("chr_snp","location_snp"))

  # --- peak counts + optional plot ---
  peak_counts <- filtered %>%
    group_by(peak_id, QTL_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(QTL_type) %>%
    summarise(num_peaks = n(), .groups = "drop")

  plot_files <- character(0)
  p_counts <- ggplot(peak_counts, aes(QTL_type, num_peaks)) +
    geom_col() +
    labs(title = "Number of peaks per QTL type", x = "QTL type", y = "Number of peaks") +
    theme_minimal()

  if (write_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    f_counts <- file.path(output_dir, paste0(prefix, "_peak_counts_per_qtl_type.pdf"))
    ggsave(f_counts, plot = p_counts, width = 8, height = 6)
    plot_files["peak_counts"] <- f_counts
  }

  # split
  filt_cis   <- filtered %>% filter(QTL_type == "cis")
  filt_trans <- filtered %>% filter(QTL_type == "trans")

  cis_ranges <- trans_ranges <- NULL
  cis_final  <- filt_cis
  trans_final <- filt_trans

  # --- peak_range (optional) ---
  if (apply_peak_range) {
    cis_ranges   <- do.call(peak_range,   c(list(qtl_data = filt_cis),   peak_params))
    trans_ranges <- do.call(peak_range,   c(list(qtl_data = filt_trans), peak_params))

    cis_ranges   <- cis_ranges   %>% rename(chr = Chromosome, location = PeakPosition)
    trans_ranges <- trans_ranges %>% rename(chr = Chromosome, location = PeakPosition)

    cis_final   <- inner_join(cis_ranges,   filt_cis,   by = c("chr","location")) %>% as_tibble()
    trans_final <- inner_join(trans_ranges, filt_trans, by = c("chr","location")) %>% as_tibble()
  }

  # --- FDR recalculation (optional) ---
  if (recalc_fdr) {
    recalc <- function(df, p_col = "p_value", out_col = "FDR_recal") {
      df[[out_col]] <- p.adjust(df[[p_col]], method = "BH"); df
    }
    cis_final   <- recalc(cis_final)
    trans_final <- recalc(trans_final)
  }

  # --- rowData annotation (optional) ---
  if (annotate_peaks) {
    cis_final   <- annotate_with_rowdata_peaks(cis_final,   se = se, key_col = key_col)
    trans_final <- annotate_with_rowdata_peaks(trans_final, se = se, key_col = key_col)
  }

  # --- nearest gene by TSS from GTF (optional) ---
  tss_used <- NULL
  if (nearest_genes) {
    if (is.null(tss_gr)) {
      if (is.null(gtf_file)) stop("nearest_genes=TRUE requires gtf_file (or tss_gr).")
      tss_used <- nearest_gene_from_gtf_build_tss(gtf_file)
    } else {
      tss_used <- tss_gr
    }

    add_nearest <- function(df) {
      if (!nrow(df)) {
        df$nearest_gene <- character(0); df$nearest_TSS_distance_bp <- integer(0); return(df)
      }
      if (nearest_on == "peak") {
        chr_vec <- df$chr
        pos_vec <- df$location
      } else { # snp
        chr_vec <- df$chr_snp
        pos_vec <- df$location_snp
      }

      if (nearest_output == "gene") {
        ng <- nearest_gene_from_gtf(gtf_file = NULL, chr = chr_vec, pos = pos_vec,
                                    output = "gene", tss_gr = tss_used)
        df$nearest_gene <- ng
      } else if (nearest_output == "distance") {
        nd <- nearest_gene_from_gtf(gtf_file = NULL, chr = chr_vec, pos = pos_vec,
                                    output = "distance", tss_gr = tss_used)
        df$nearest_TSS_distance_bp <- as.integer(nd)
      } else { # both
        nb <- nearest_gene_from_gtf(gtf_file = NULL, chr = chr_vec, pos = pos_vec,
                                    output = "both", tss_gr = tss_used)
        # nb has: query_chr, query_pos, gene_name, gene_strand, tss_pos, distance_bp, signed_distance_bp
        df$nearest_gene               <- nb$gene_name
        df$nearest_gene_strand        <- nb$gene_strand
        df$nearest_TSS_pos            <- as.integer(nb$tss_pos)
        df$nearest_TSS_distance_bp    <- as.integer(nb$distance_bp)
        df$nearest_TSS_signed_bp      <- as.integer(nb$signed_distance_bp)
      }
      df
    }

    cis_final   <- add_nearest(cis_final)
    trans_final <- add_nearest(trans_final)
  }

  # --- index plots (optional) ---
  if (write_outputs) {
    f_cis   <- file.path(output_dir, paste0(prefix, "_qtl_index_cis.pdf"))
    f_trans <- file.path(output_dir, paste0(prefix, "_qtl_index_trans.pdf"))
    plot_qtl_index(cis_final,   use_neglog10 = TRUE, fdr_recal_col = if (recalc_fdr) "FDR_recal" else NULL, out_file = f_cis)
    plot_qtl_index(trans_final, use_neglog10 = TRUE, fdr_recal_col = if (recalc_fdr) "FDR_recal" else NULL, out_file = f_trans)
    plot_files["cis_index"]   <- f_cis
    plot_files["trans_index"] <- f_trans
  }

  # --- CSV outputs (optional) ---
  if (write_outputs) {
    if (!is.null(cis_ranges))   write.csv(cis_ranges,   file.path(output_dir, paste0(prefix, "_cis_ranges.csv")),   row.names = FALSE)
    if (!is.null(trans_ranges)) write.csv(trans_ranges, file.path(output_dir, paste0(prefix, "_trans_ranges.csv")), row.names = FALSE)
    write.csv(cis_final,   file.path(output_dir, paste0(prefix, "_cis_results.csv")),   row.names = FALSE)
    write.csv(trans_final, file.path(output_dir, paste0(prefix, "_trans_results.csv")), row.names = FALSE)
  }

  list(
    filtered     = filtered,
    peak_counts  = peak_counts,
    cis_ranges   = cis_ranges,
    trans_ranges = trans_ranges,
    cis          = cis_final,
    trans        = trans_final,
    plot_files   = plot_files,
    tss          = tss_used
  )
}
