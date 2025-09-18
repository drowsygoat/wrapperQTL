#!/usr/bin/env Rscript

# ==============================
# eQTL summarizer (no dbscan)
# ==============================

suppressPackageStartupMessages({
  require(data.table)
  require(tidyverse)
  require(gtools)   # mixedsort
  require(ggplot2)
  require(cowplot)
  require(ggpubr)
  require(argparse)
})

options(datatable.showProgress = TRUE)

# ---------- utilities ----------

vlog <- function(..., verbose = FALSE) {
  if (isTRUE(verbose)) message(format(Sys.time(), "%H:%M:%S "), paste0(...))
}

exists_and_not_empty <- function(path) {
  dir.exists(path) && length(list.files(path, all.files = FALSE, no.. = TRUE)) > 0
}

# Robust -log10 transform that avoids -Inf and NA
neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))

# ---------- CLI args ----------

parse_args <- function() {
  parser <- ArgumentParser(description = "Summarize cis/trans eQTL results in all directories matched by a glob pattern.")
  parser$add_argument("dir_pattern", nargs = "+",
      help = "One or more glob patterns (e.g., '/data/run_*' or multiple patterns).")
  parser$add_argument("--prefix", default = "defaultprefix",
                      help = "Prefix used to select files and to name outputs (default: 'defaultprefix').")
  parser$add_argument("--p_value_threshold", default = "NULL",
                      help = "Numeric p-value filter (use 'NULL' to skip).")
  parser$add_argument("--threshold", default = "FDR",
                      choices = c("FDR", "FWER"),
                      help = "Multiple-testing metric to compute and use for binning (default: FDR).")
  parser$add_argument("--verbose", action = "store_true", help = "Print progress messages.")
  parser$add_argument("--cis_pattern", default = "eQTL_cis\\.txt\\.gz", help = "Regex for cis files.")
  parser$add_argument("--trans_pattern", default = "eQTL_trans\\.txt\\.gz", help = "Regex for trans files.")
  parser$add_argument("--ntests_pattern", default = "ntests_combined\\.txt", help = "Regex for n-tests file.")
  args <- parser$parse_args()

  # normalize p-value threshold
  if (tolower(args$p_value_threshold) %in% c("null", "na", "none", "")) {
    args$p_value_threshold <- NULL
  } else {
    args$p_value_threshold <- as.numeric(args$p_value_threshold)
  }
  args
}

# ---------- n-tests reader ----------

read_ntests <- function(directory, pattern, prefix, verbose = FALSE) {
  paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  paths <- grep(prefix, paths, value = TRUE)
  if (length(paths) == 0) stop("No ntests_combined.txt file matched in: ", directory)

  vlog("Reading n-tests from: ", paste(basename(paths), collapse = ", "), verbose = verbose)
  nt <- fread(paths[1])
  # Expect 3 columns like: 1325010769940 1322420431298 2590338642
  nt <- lapply(nt, as.double)
  list(n_tests_cis = nt[[3]], n_tests_trans = nt[[2]])
}

# ---------- single file loader (no FDR/FWER here) ----------

load_and_filter <- function(file, p_val_threshold, threshold, n_tests_for_p_adjust, verbose = FALSE) {
  vlog("Loading ", basename(file), verbose = verbose)

  n_lines <- R.utils::countLines(file)[[1]]
  if (n_lines <= 2) {
    vlog("  Skipping (≤2 lines).", verbose = verbose)
    return(NULL)
  }

  # Files look like:
  # SNP  gene  beta  t-stat  p-value
  # Use header=TRUE to be robust to column order; then standardize names.
  dt <- tryCatch(
    fread(file, header = TRUE, sep = "\t", showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt)) {
    # Fallback to old positional read (skip header)
    dt <- fread(file, header = FALSE, skip = 1, select = c(1, 2, 4, 5))
    setnames(dt, c("V1", "V2", "V4", "V5"), c("snp_id", "peak_id", "beta", "p_value"))
  } else {
    # Standardize names
    setnames(dt, make.names(names(dt)))
    # Expect columns something like: SNP, gene, beta, t.stat, p.value
    required <- c("SNP", "gene", "p.value")
    if (!all(required %in% names(dt))) {
      stop("Unexpected columns in file: ", file, "\nHave: ", paste(names(dt), collapse = ", "))
    }
    # beta can be missing in some outputs; tolerate by creating if absent
    if (!"beta" %in% names(dt)) dt[, beta := NA_real_]

    setnames(dt, c("SNP", "gene", "p.value"), c("snp_id", "peak_id", "p_value"))
    dt <- dt[, .(snp_id, peak_id, beta, p_value)]
  }

  if (nrow(dt) < 1L) return(NULL)

  # Optional early p-value prefilter ONLY (to reduce memory); no FDR/FWER yet.
  if (!is.null(p_val_threshold)) {
    before <- nrow(dt)
    dt <- dt[p_value < p_val_threshold]
    vlog("  p-value filter kept ", nrow(dt), "/", before, " rows.", verbose = verbose)
  }

  # DO NOT compute FDR/FWER here anymore — this now happens AFTER merging all files.
  vlog("  Loaded ", nrow(dt), " rows (no multiple-testing adjustment yet).", verbose = verbose)
  dt
}

# ---------- batch processing for one directory ----------

process_eqtl_files <- function(directory, results_path,
                               cis_pattern, trans_pattern, ntests_pattern,
                               threshold, prefix, p_val_threshold,
                               verbose = FALSE) {

  vlog("Reading n-tests…", verbose = verbose)
  n_tests <- read_ntests(directory, ntests_pattern, prefix, verbose)
  n_trans <- n_tests$n_tests_trans
  n_cis   <- n_tests$n_tests_cis

  if (any(c(n_trans, n_cis) <= 1)) stop("Invalid n_tests value (<= 1).")

  cis_files   <- grep(prefix, list.files(directory, pattern = cis_pattern,   full.names = TRUE), value = TRUE)
  trans_files <- grep(prefix, list.files(directory, pattern = trans_pattern, full.names = TRUE), value = TRUE)

  dir.create(file.path(results_path, "rds"), recursive = TRUE, showWarnings = FALSE)

  # Load (and optional p-value prefilter) — NO FDR/FWER yet
  for (f in trans_files) {
    d <- load_and_filter(f, p_val_threshold, threshold, n_trans, verbose)
    if (!is.null(d)) {
      saveRDS(d, file = file.path(results_path, "rds", paste0("trans_", basename(f), ".rds")))
      rm(d); gc()
    }
  }

  for (f in cis_files) {
    d <- load_and_filter(f, p_val_threshold, threshold, n_cis, verbose)
    if (!is.null(d)) {
      saveRDS(d, file = file.path(results_path, "rds", paste0("cis_", basename(f), ".rds")))
      rm(d); gc()
    }
  }

  # Combine all loaded chunks per type
  trans_rds <- list.files(file.path(results_path, "rds"), pattern = "^trans_.*\\.rds$", full.names = TRUE)
  cis_rds   <- list.files(file.path(results_path, "rds"), pattern = "^cis_.*\\.rds$",   full.names = TRUE)

  eqtl_trans <- if (length(trans_rds)) rbindlist(lapply(trans_rds, readRDS), use.names = TRUE, fill = TRUE) else data.table()
  eqtl_cis   <- if (length(cis_rds))   rbindlist(lapply(cis_rds,   readRDS), use.names = TRUE, fill = TRUE) else data.table()

  # ---------- NOW compute multiple-testing adjustment on the merged sets ----------
  if (nrow(eqtl_trans)) {
    if (threshold == "FWER") {
      eqtl_trans[, FWER := p.adjust(p_value, method = "FWER", n = n_trans)]
      eqtl_trans <- eqtl_trans[FWER < 0.1]
    } else {
      eqtl_trans[, FDR := p.adjust(p_value, method = "BH", n = n_trans)]
      eqtl_trans <- eqtl_trans[FDR < 0.1]
    }
  }

  if (nrow(eqtl_cis)) {
    if (threshold == "FWER") {
      eqtl_cis[, FWER := p.adjust(p_value, method = "FWER", n = n_cis)]
      eqtl_cis <- eqtl_cis[FWER < 0.1]
    } else {
      eqtl_cis[, FDR := p.adjust(p_value, method = "BH", n = n_cis)]
      eqtl_cis <- eqtl_cis[FDR < 0.1]
    }
  }
  # -------------------------------------------------------------------------------

  if (nrow(eqtl_trans)) eqtl_trans[, QTL_type := "trans"]
  if (nrow(eqtl_cis))   eqtl_cis[,   QTL_type := "cis"]

  sig_cis   <- nrow(eqtl_cis)
  sig_trans <- nrow(eqtl_trans)

  overall_numbers <- data.table(
    Type = c("cis", "trans"),
    Significant = c(sig_cis, sig_trans),
    Total = c(n_cis, n_trans)
  )[, Proportion := fifelse(Total > 0, Significant / Total, 0)]

  saveRDS(overall_numbers, file = file.path(results_path, paste(prefix, "overall_numbers.rds", sep = "_")))

  # Merge for downstream summaries (post-filter, post-adjustment)
  bon <- rbindlist(list(eqtl_cis, eqtl_trans), use.names = TRUE, fill = TRUE)
  if (!nrow(bon)) {
    vlog("No significant results found after adjustment/filtering.", verbose = verbose)
    saveRDS(bon, file = file.path(results_path, paste(prefix, "results.rds", sep = "_")))
    return(invisible(NULL))
  }

  # Location parsing & QC
  bon[, c("chr_location", "identity_snp") := tstrsplit(snp_id, "_", fixed = TRUE)]
  bon[, c("chr_snp", "location_snp") := tstrsplit(chr_location, ":", fixed = TRUE)]
  bon[, location := as.integer(location_snp)]
  bon[, chr_snp := paste0("chr", chr_snp)]
  bon <- bon[nchar(chr_snp) < 6]  # drop weird chromosomes

  check_column <- if (threshold == "FDR") "FDR" else "FWER"

  # Save any non-finite after transform for debugging
  nf <- bon[!is.finite(neglog10(get(check_column)))]
  if (nrow(nf) > 0) {
    vlog("Found ", nrow(nf), " non-finite ", check_column, " values; saving for debug.", verbose = verbose)
    saveRDS(nf, file = file.path(directory, paste0("non_finite_values_debug_", basename(results_path), ".rds")))
  }

  # Bins of -log10(metric)
  bin_col <- if (threshold == "FDR") "FDR_bins" else "FWER_bins"
  bon[, (bin_col) := cut_width(neglog10(get(check_column)), width = 1, center = 0.5)]

  # Order chromosomes nicely
  bon[, chr_snp := factor(chr_snp, levels = gtools::mixedsort(unique(chr_snp)))]
  setorder(bon, chr_snp)

  saveRDS(bon, file = file.path(results_path, paste(prefix, "results.rds", sep = "_")))
}

# ---------- summaries for plotting (no dbscan) ----------

make_all_plots_data <- function(results_path, prefix, threshold = "FDR", verbose = FALSE) {
  vlog("Preparing plot data…", verbose = verbose)

  bin_col <- if (threshold == "FDR") "FDR_bins" else "FWER_bins"
  bon <- readRDS(file.path(results_path, paste(prefix, "results.rds", sep = "_")))
  if (!nrow(bon)) return(invisible(NULL))

  # How many QTLs per SNP (by type)
  snp_target_stats <- bon[, .(targets_count = uniqueN(peak_id)), by = .(snp_id, QTL_type)]
  saveRDS(snp_target_stats, file = file.path(results_path, paste(prefix, "snp_target_stats.rds", sep = "_")))

  # Per-chromosome counts by significance bin
  plotdata_sum_chr <- bon[, .(QTL_count = .N), by = .(chr_snp, QTL_type, bins = get(bin_col))]
  saveRDS(plotdata_sum_chr, file = file.path(results_path, paste(prefix, "plotdata_sum_chr.rds", sep = "_")))

  invisible(NULL)
}

# ---------- plotting ----------

plot_results <- function(input_dir, output_dir, threshold = "FDR", prefix, verbose = FALSE) {
  vlog("Plotting results…", verbose = verbose)

  plotdata_sum_chr <- readRDS(file.path(input_dir, paste(prefix, "plotdata_sum_chr.rds", sep = "_")))
  snp_target_stats <- readRDS(file.path(input_dir, paste(prefix, "snp_target_stats.rds", sep = "_")))
  overall_numbers  <- readRDS(file.path(input_dir, paste(prefix, "overall_numbers.rds", sep = "_")))

  bin_col <- if (threshold == "FDR") "FDR_bins" else "FWER_bins"

  common_theme <- theme(
    legend.key.size = unit(0.3, "lines"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  )

  # Proportion bar
  plot0 <- ggplot(overall_numbers, aes(x = Type, y = Proportion, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = Significant), vjust = -0.5, size = 5) +
    labs(title = "Proportion of Significant eQTLs", x = "eQTL Type", y = "Proportion of Significant eQTLs") +
    theme_minimal(base_size = 10) + common_theme

  # Per-chromosome counts by bin
  plot2 <- ggplot(plotdata_sum_chr, aes(x = chr_snp, y = QTL_count, fill = bins)) +
    geom_bar(stat = "identity") +
    labs(title = "eQTL Counts per Chromosome", x = "Chromosome", y = "Count") +
    guides(fill = guide_legend(
      title = bin_col,
      title.position = "top",
      label.position = "left",
      ncol = 3,
      keywidth = 0.2,
      keyheight = 0.2,
      default.unit = "cm"
    )) +
    theme_minimal(base_size = 10) + common_theme

  # SNP → number of targets (linear & log variants)
  plot3 <- ggplot(snp_target_stats[snp_target_stats$targets_count > 1, ], aes(x = targets_count)) +
    geom_histogram(binwidth = 1) +
    labs(title = "Targets per SNP", x = "Targets Count", y = "Frequency") +
    theme_minimal(base_size = 10) + common_theme

  plot4 <- plot3 + scale_y_log10() + labs(title = "Targets per SNP (log-scaled Y)")

  plotlist <- list(plot0, plot2, plot3, plot4)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(output_dir, "combined_plot.pdf")
  ggpubr::ggexport(plotlist, filename = outfile)

  vlog("Saved plots to: ", outfile, verbose = verbose)
}

# ---------- main ----------

main <- function() {
  args <- parse_args()
  patterns    <- args$dir_pattern
  prefix      <- args$prefix
  p_thr       <- args$p_value_threshold
  threshold   <- args$threshold
  verbose     <- args$verbose

  # Expand the glob pattern to directories
  candidates <- unique(unlist(lapply(patterns, Sys.glob)))
  candidates <- candidates[dir.exists(candidates)]

  if (length(candidates) == 0) {
    stop("No directories matched by pattern(s): ", paste(patterns, collapse = " | "))
  }

  vlog("Matched ", length(candidates), " director", ifelse(length(candidates)==1,"y","ies"),
       ":\n  - ", paste(candidates, collapse = "\n  - "), verbose = verbose)

  errors <- list()
  for (d in candidates) {
    vlog("=== Processing directory: ", d, " ===", verbose = verbose)
    # directory existence already assured
    results_path <- file.path(d, paste(prefix, "results", sep = "_"))
    plots_path   <- file.path(d, paste(prefix, "plots",   sep = "_"))
    dir.create(file.path(results_path, "rds"), recursive = TRUE, showWarnings = FALSE)

    # Try per-directory run so a failure doesn't stop the batch
    tryCatch(
      {
        vlog("Starting… prefix=", prefix, " threshold=", threshold,
             " p-value-threshold=", ifelse(is.null(p_thr), "NULL", p_thr), verbose = verbose)

        process_eqtl_files(
          directory      = d,
          results_path   = results_path,
          trans_pattern  = args$trans_pattern,
          cis_pattern    = args$cis_pattern,
          ntests_pattern = args$ntests_pattern,
          prefix         = prefix,
          threshold      = threshold,
          p_val_threshold= p_thr,
          verbose        = verbose
        )

        make_all_plots_data(
          results_path = results_path,
          prefix       = prefix,
          threshold    = threshold,
          verbose      = verbose
        )

        plot_results(
          input_dir = results_path,
          output_dir = plots_path,
          threshold = threshold,
          prefix = prefix,
          verbose = verbose
        )

        vlog("=== Done: ", d, " ===", verbose = verbose)
      },
      error = function(e) {
        msg <- paste0("Directory failed: ", d, " -> ", conditionMessage(e))
        message(msg)
        errors[[length(errors) + 1]] <<- msg
      }
    )
  }

  if (length(errors)) {
    message("\nCompleted with errors in ", length(errors), " director",
            ifelse(length(errors)==1,"y","ies"), ":\n  - ",
            paste(errors, collapse = "\n  - "))
  } else {
    message("\nCompleted successfully for all ", length(candidates), " directories.")
  }
}

# run
main()
