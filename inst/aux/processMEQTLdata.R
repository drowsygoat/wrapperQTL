#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(tidyverse)
  library(gtools)   # mixedsort
  library(ggplot2)
  library(ggpubr)
  library(argparse)
})

options(datatable.showProgress = TRUE)

num_cores <- 1
message("Detected cores: ", detectCores())
message("Using cores: ", num_cores)

# -------------------- CLI --------------------
parse_args <- function() {
  parser <- ArgumentParser(description = "Process eQTL files (no hotspots)")
  parser$add_argument("directory", help = "Directory containing the files")
  parser$add_argument("--prefix",
                      default = "defaultprefix",
                      help = "Prefix to select files and name outputs [default: %(default)s]")
  parser$add_argument("--p_value_threshold",
                      default = "",
                      help = "Optional raw p-value prefilter (numeric). Omit or empty to skip.")
  args <- parser$parse_args()

  # normalize threshold
  if (is.null(args$p_value_threshold) || args$p_value_threshold == "") {
    args$p_value_threshold <- NULL
  } else {
    args$p_value_threshold <- as.numeric(args$p_value_threshold)
  }
  args
}

exists_and_not_empty <- function(directory_path) {
  if (!dir.exists(directory_path)) return(FALSE)
  length(list.files(directory_path)) > 0
}

# -------------------- helpers --------------------
read_ntests <- function(directory, pattern, prefix) {
  ntests_file <- list.files(directory, pattern = pattern, full.names = TRUE)
  ntests_file <- grep(prefix, ntests_file, value = TRUE)

  if (length(ntests_file) == 0) {
    stop("No ntests_combined.txt file found in the directory (after prefix filter).")
  }
  n_tests <- fread(ntests_file[1], header = FALSE)
  n_tests <- lapply(n_tests, as.double)
  list(n_tests_cis = n_tests[[3]], n_tests_trans = n_tests[[2]])
}

load_and_filter <- function(file, p_val_threshold = NULL, threshold = "FDR", n_tests_for_p_adjust) {
  message("load_and_filter: ", basename(file))

  num_lines <- R.utils::countLines(file)[[1]]
  if (num_lines <= 2) {
    message("Skipping file due to insufficient lines: ", file)
    return(NULL)
  }

  dt <- fread(
    file,
    header = FALSE,    # file expected tab-separated with header line we will skip
    skip = 1,
    showProgress = TRUE,
    select = c(1, 2, 4, 5)
  )

  setnames(dt, c("V1","V2","V4","V5"), c("snp_id","peak_id","beta","p_value"))

  if (nrow(dt) < 2) {
    message("Skipping file due to insufficient data: ", file)
    return(NULL)
  }

  if (!is.null(p_val_threshold)) {
    dt <- dt[p_value < p_val_threshold]
  }

  if (threshold == "FWER") {
    dt[, FWER := p.adjust(p_value, method = "bonferroni", n = n_tests_for_p_adjust)]
    dt <- dt[FWER < 0.1]
  } else {
    dt[, FDR := p.adjust(p_value, method = "BH", n = n_tests_for_p_adjust)]
    dt <- dt[FDR < 0.1]
  }

  if (nrow(dt) == 0) return(NULL)
  dt
}

process_eqtl_files <- function(directory,
                               results_path,
                               cis_pattern,
                               trans_pattern,
                               n_tests_pattern,
                               threshold = "FDR",
                               prefix = "defaultprefix",
                               p_val_threshold = NULL) {
  message("process_eqtl_files")

  n_tests <- read_ntests(directory, pattern = n_tests_pattern, prefix = prefix)
  n_tests_trans <- n_tests$n_tests_trans
  n_tests_cis   <- n_tests$n_tests_cis

  if (any(c(n_tests_trans, n_tests_cis) <= 1)) {
    stop("Invalid n_tests value (<= 1).")
  }

  trans_files <- list.files(directory, pattern = trans_pattern, full.names = TRUE)
  cis_files   <- list.files(directory, pattern = cis_pattern,   full.names = TRUE)

  trans_files <- grep(prefix, trans_files, value = TRUE)
  cis_files   <- grep(prefix, cis_files,   value = TRUE)

  dir.create(file.path(results_path, "rds"), recursive = TRUE, showWarnings = FALSE)

  # TRANS
  for (f in trans_files) {
    d <- load_and_filter(f, p_val_threshold = p_val_threshold,
                         threshold = threshold, n_tests_for_p_adjust = n_tests_trans)
    if (!is.null(d)) {
      saveRDS(d, file = file.path(results_path, "rds", paste0("trans_", basename(f), ".rds")))
      rm(d); gc()
    }
  }

  # CIS
  for (f in cis_files) {
    d <- load_and_filter(f, p_val_threshold = p_val_threshold,
                         threshold = threshold, n_tests_for_p_adjust = n_tests_cis)
    if (!is.null(d)) {
      saveRDS(d, file = file.path(results_path, "rds", paste0("cis_", basename(f), ".rds")))
      rm(d); gc()
    }
  }

  # Load and combine
  trans_rds_files <- list.files(file.path(results_path, "rds"), pattern = "^trans_.*\\.rds$", full.names = TRUE)
  cis_rds_files   <- list.files(file.path(results_path, "rds"), pattern = "^cis_.*\\.rds$",   full.names = TRUE)

  eqtl_data_trans <- if (length(trans_rds_files)) rbindlist(lapply(trans_rds_files, readRDS), use.names = TRUE, fill = TRUE) else data.table()
  eqtl_data_cis   <- if (length(cis_rds_files))   rbindlist(lapply(cis_rds_files,   readRDS), use.names = TRUE, fill = TRUE) else data.table()

  if (nrow(eqtl_data_trans)) eqtl_data_trans[, QTL_type := "trans"]
  if (nrow(eqtl_data_cis))   eqtl_data_cis[,   QTL_type := "cis"]

  sig_cis   <- nrow(eqtl_data_cis)
  sig_trans <- nrow(eqtl_data_trans)

  overall_numbers <- data.table(
    Type        = c("cis", "trans"),
    Significant = c(sig_cis, sig_trans),
    Total       = c(n_tests_cis, n_tests_trans)
  )
  overall_numbers[, Proportion := fifelse(Total > 0, Significant / Total, 0)]

  saveRDS(overall_numbers, file = file.path(results_path, paste(prefix, "overall_numbers.rds", sep = "_")))

  # Prepare data for non-hotspot plots
  bon_plotdata <- rbindlist(list(eqtl_data_cis, eqtl_data_trans), use.names = TRUE, fill = TRUE)
  if (!nrow(bon_plotdata)) {
    warning("No significant results after filtering/adjustment.")
    saveRDS(bon_plotdata, file = file.path(results_path, paste(prefix, "results.rds", sep = "_")))
    return(invisible(NULL))
  }

  bon_plotdata[, ppts := qunif(ppoints(.N)), by = QTL_type]
  bon_plotdata[, c("chr_location", "identity_snp") := tstrsplit(snp_id, "_", fixed = TRUE)]
  bon_plotdata[, c("chr_snp", "location_snp") := tstrsplit(chr_location, ":", fixed = TRUE)]
  bon_plotdata[, location := as.integer(location_snp)]
  bon_plotdata[, chr_snp := paste0("chr", chr_snp)]
  bon_plotdata <- bon_plotdata[nchar(chr_snp) < 6]

  check_column <- ifelse(threshold == "FDR", "FDR", "FWER")
  non_finite_rows <- bon_plotdata[!is.finite(-log10(get(check_column)))]
  if (nrow(non_finite_rows) > 0) {
    message("Non-finite adjusted p-values detected; saving debug rows.")
    saveRDS(non_finite_rows, file = file.path(directory, paste0("non_finite_values_debug_", basename(results_path), ".rds")))
  }

  bin_column <- ifelse(threshold == "FDR", "FDR_bins", "FWER_bins")
  bon_plotdata[, (bin_column) := cut_width(-log10(get(check_column)), width = 1, center = 0.5)]

  sorted_levels <- gtools::mixedsort(unique(bon_plotdata$chr_snp))
  bon_plotdata[, chr_snp := factor(chr_snp, levels = sorted_levels)]
  setorder(bon_plotdata, chr_snp)

  saveRDS(bon_plotdata, file = file.path(results_path, paste(prefix, "results.rds", sep = "_")))

  # Precompute simple summaries (no hotspots)
  plotdata_sum_chr <- bon_plotdata[, .(QTL_count = .N), by = .(chr_snp, QTL_type, bins = get(bin_column))]
  saveRDS(plotdata_sum_chr, file = file.path(results_path, paste(prefix, "plotdata_sum_chr.rds", sep = "_")))

  snp_target_stats <- bon_plotdata[, .(targets_count = uniqueN(peak_id)), by = .(snp_id, QTL_type)]
  saveRDS(snp_target_stats, file = file.path(results_path, paste(prefix, "snp_target_stats.rds", sep = "_")))

  invisible(NULL)
}

plot_results <- function(input_dir, output_dir = NULL, threshold = "FDR", prefix = "defaultprefix") {
  message("plot_results")

  overall_numbers  <- readRDS(file.path(input_dir, paste(prefix, "overall_numbers.rds", sep = "_")))
  plotdata_sum_chr <- readRDS(file.path(input_dir, paste(prefix, "plotdata_sum_chr.rds", sep = "_")))
  snp_target_stats <- readRDS(file.path(input_dir, paste(prefix, "snp_target_stats.rds", sep = "_")))

  bin_column <- ifelse(threshold == "FDR", "FDR_bins", "FWER_bins")

  common_theme <- theme(
    legend.key.size = unit(0.3, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )

  plot0 <- ggplot(overall_numbers, aes(x = Type, y = Proportion, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Significant), vjust = -0.5, size = 5) +
    labs(title = "Proportion of Significant eQTLs", x = "eQTL Type", y = "Proportion of Significant eQTLs")

  plot2 <- ggplot(plotdata_sum_chr, aes(x = chr_snp, y = QTL_count, fill = bins)) +
    geom_bar(stat = "identity") +
    labs(title = "Peak Values per Chromosome", x = "Chromosome", y = "Count") +
    guides(fill = guide_legend(
      title = bin_column,
      title.position = "top",
      label.position = "left",
      ncol = 3,
      keywidth = 0.2,
      keyheight = 0.2,
      default.unit = "cm"
    )) + theme_minimal(base_size = 10) + common_theme

  plot3 <- ggplot(snp_target_stats[snp_target_stats$targets_count > 1, ], aes(x = targets_count)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Histogram of Peak Targets Count per SNP", x = "Targets Count", y = "Frequency")

  plot4 <- plot3 + scale_y_log10()

  plotlist <- lapply(list(plot0, plot2, plot3, plot4), function(p) p + theme_minimal(base_size = 10) + common_theme)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggpubr::ggexport(plotlist, filename = file.path(output_dir, "combined_plot.pdf"))
  }

  invisible(NULL)
}

# -------------------- main --------------------
main <- function() {
  if (!interactive()) {
    args <- parse_args()
    directory <- args$directory
    prefix <- args$prefix
    p_val_threshold <- args$p_value_threshold
  } else {
    directory <- getwd()
    prefix <- "defaultprefix"
    p_val_threshold <- NULL
  }

  results_path <- file.path(directory, paste(prefix, "results", sep = "_"))

  if (!exists_and_not_empty(results_path)) {
    dir.create(file.path(results_path, "rds"), recursive = TRUE, showWarnings = FALSE)
    process_eqtl_files(
      directory = directory,
      results_path = results_path,
      trans_pattern = "eQTL_trans\\.txt\\.gz",
      cis_pattern   = "eQTL_cis\\.txt\\.gz",
      n_tests_pattern = "ntests_combined\\.txt",
      prefix = prefix,
      threshold = "FDR",
      p_val_threshold = p_val_threshold
    )
  } else {
    message("Directory with results exists. Skipping processing step.")
  }

  plots_path <- file.path(directory, paste(prefix, "plots", sep = "_"))
  if (!exists_and_not_empty(plots_path)) {
    dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)
    plot_results(input_dir = results_path, output_dir = plots_path, threshold = "FDR", prefix = prefix)
  } else {
    message("Directory with plots exists. Skipping plotting step.")
  }
}

# Run
main()
