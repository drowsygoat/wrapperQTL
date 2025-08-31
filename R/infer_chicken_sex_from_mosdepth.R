#' Infer biological sex from mosdepth Z and W chromosome coverage
#'
#' @title Infer Sex (ZW system) from mosdepth summaries and optionally write outputs
#' @description
#' Reads all `summary.txt` files under a mosdepth output directory, extracts mean
#' coverage for chromosomes Z and W, computes the Z:W coverage ratio per sample,
#' and calls sex under the avian ZW system (female = ZW, male = ZZ).
#'
#' @param mosdepth_dir Character. Path to a directory containing mosdepth outputs.
#' @param out_dir Character. Directory where output files (`covariates_sex.txt` and
#'   `ZW_coverage_ratio.pdf`) will be written. Default: current working directory.
#' @param threshold Numeric. Z:W ratio cutoff (default `1.5`).
#' @param chrom_labels Character vector of length 2 or named list with `"Z"` and `"W"`.
#' @param verbose Logical. Print messages (default `TRUE`).
#'
#' @return A list with:
#' \itemize{
#'   \item `results`: data.frame with per-sample Z, W, ratio, and sex call.
#'   \item `covariates`: one-row data.frame suitable for Matrix eQTL.
#'   \item `files`: character vector of summary files used.
#'   \item `plot_file`: path to saved PDF (if written).
#'   \item `covariate_file`: path to covariate file (if written).
#' }
#'
#' @export
infer_sex_from_mosdepth <- function(
  mosdepth_dir,
  out_dir = ".",
  threshold = 1.5,
  chrom_labels = c("Z", "W"),
  verbose = TRUE
) {
  if (!dir.exists(mosdepth_dir)) {
    stop("Directory does not exist: ", mosdepth_dir)
  }
  if (!dir.exists(out_dir)) {
    if (verbose) message("Creating output directory: ", out_dir)
    dir.create(out_dir, recursive = TRUE)
  }

  summary_files <- list.files(mosdepth_dir,
                              pattern = "summary\\.txt$",
                              recursive = TRUE, full.names = TRUE)
  if (length(summary_files) == 0) {
    stop("No mosdepth summary files found in: ", mosdepth_dir)
  }

  # Handle chrom_labels
  get_label_set <- function(x) {
    if (is.list(x)) {
      if (!all(c("Z","W") %in% names(x)))
        stop("If chrom_labels is a list, must be named list(Z=..., W=...).")
      list(Z = as.character(x$Z), W = as.character(x$W))
    } else {
      x <- as.character(x)
      if (length(x) != 2) stop("chrom_labels must have length 2.")
      names(x) <- c("Z","W")
      list(Z = x[1], W = x[2])
    }
  }
  clabs <- get_label_set(chrom_labels)

  results <- data.frame(Sample = character(),
                        Z = numeric(), W = numeric(),
                        ZW_ratio = numeric(),
                        stringsAsFactors = FALSE)

  for (file in summary_files) {
    sample_name <- basename(dirname(file))
    df <- utils::read.table(file, header = TRUE)
    z_cov <- df$mean[df$chrom %in% clabs$Z]
    w_cov <- df$mean[df$chrom %in% clabs$W]

    if (length(z_cov) == 1 && length(w_cov) == 1 && w_cov > 0) {
      results <- rbind(results, data.frame(
        Sample = sample_name,
        Z = z_cov,
        W = w_cov,
        ZW_ratio = z_cov / w_cov,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Call sex
  results$Sex <- ifelse(results$ZW_ratio > threshold, "female", "male")
  results$Sex_numeric <- ifelse(results$Sex == "female", 1, 0)

  # Output file paths
  plot_file <- file.path(out_dir, "ZW_coverage_ratio.pdf")
  covariate_file <- file.path(out_dir, "covariates_sex.txt")

  # Save PDF
  grDevices::pdf(plot_file, width = 8, height = 5)
  print(
    ggplot2::ggplot(results,
           ggplot2::aes(x = Sample, y = ZW_ratio, fill = Sex)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
      ggplot2::labs(title = "Z:W Chromosome Coverage Ratio (female = ZW)",
           x = "Sample", y = "Z / W Mean Coverage Ratio") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::scale_fill_manual(values = c("female" = "orange", "male" = "steelblue"))
  )
  grDevices::dev.off()
  if (verbose) message("✅ Plot saved to: ", normalizePath(plot_file))

  # Save covariates
  covariate_df <- data.frame(row.names = "sex", t(results$Sex_numeric))
  colnames(covariate_df) <- results$Sample
  utils::write.table(covariate_df, file = covariate_file,
              sep = "\t", quote = FALSE, col.names = NA)
  if (verbose) message("✅ Covariates saved to: ", normalizePath(covariate_file))

  list(results = results,
       covariates = covariate_df,
       files = summary_files,
       plot_file = plot_file,
       covariate_file = covariate_file)
}