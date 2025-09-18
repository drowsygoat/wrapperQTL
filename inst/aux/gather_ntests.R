#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
})

# Parse CLI arguments only
parse_args <- function() {
  parser <- argparse::ArgumentParser(description = "Process eQTL files")

  parser$add_argument("directory", help = "Directory containing the files")
  parser$add_argument("--pattern",
                      default = "ntests\\.txt",
                      help = "Regex pattern to match files [default: %(default)s]")
  parser$add_argument("--prefix",
                      default = "defaultprefix",
                      help = "Prefix to filter files [default: %(default)s]")

  parser$parse_args()
}

# Gather and Summarize eQTL Test Counts
gather_ntests <- function(directory, pattern, prefix) {
  files <- list.files(path = directory, pattern = pattern, full.names = TRUE)
  files <- grep(prefix, files, value = TRUE)

  if (length(files) < 1) {
    stop("No files found to gather.")
  }

  data_list <- list()
  file_counter <- 0

  for (file in files) {
    file_counter <- file_counter + 1
    cat("Processing file", file_counter, ":", file, "\n")

    data <- read.table(file, header = FALSE, sep = " ")
    data <- data[, seq_len(min(3, ncol(data))), drop = FALSE]
    data[] <- lapply(data, as.numeric)

    data_list[[length(data_list) + 1]] <- data
  }

  combined_data <- dplyr::bind_rows(data_list)

  if (anyNA(combined_data)) {
    warning("NA values were present in the data and were removed")
  }

  summed_data <- colSums(combined_data, na.rm = TRUE)
  if (length(summed_data) < 3) summed_data <- c(summed_data, rep(0, 3 - length(summed_data)))

  result <- data.frame(all = summed_data[1], trans = summed_data[2], cis = summed_data[3])

  out_path <- file.path(directory, paste(prefix, "ntests_combined.txt", sep = "_"))
  write.table(result, file = out_path, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)

  return(result)
}

# Main CLI entry point
main <- function() {
  args <- parse_args()
  result <- gather_ntests(args$directory, args$pattern, args$prefix)
  print(result)
}

# Force CLI-only execution
if (!interactive()) {
  main()
} else {
  stop("This script is intended to be run from the command line, not interactively.")
}
