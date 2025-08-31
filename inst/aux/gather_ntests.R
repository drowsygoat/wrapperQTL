#!/usr/bin/env Rscript

parse_args <- function() {
  parser <- argparse::ArgumentParser(description = "Process eQTL files")

  # required positional
  parser$add_argument("directory", help = "Directory containing the files")

  # optional flags with defaults
  parser$add_argument("--pattern",
                      default = "ntests\\.txt",
                      help = "Regex pattern to match files [default: %(default)s]")
  parser$add_argument("--prefix",
                      default = "defaultprefix",
                      help = "Prefix to filter files [default: %(default)s]")

  args <- parser$parse_args()
  return(args)
}

#' Gather and Summarize eQTL Test Counts
#'
#' Reads multiple "ntests" files from a directory, filters by pattern and prefix,
#' and outputs the sum of test counts (all, trans, cis).
#'
#' @param directory A character string specifying the path to the directory containing files.
#' @param pattern A regex pattern to match files (default: "ntests\\\\.txt").
#' @param prefix A string to further filter matched files (default: "defaultprefix").
#'
#' @return A data.frame with one row and columns: `all`, `trans`, `cis`.
#' Also writes the result to a file named `{prefix}_ntests_combined.txt` in the given directory.
#'
#' @examples
#' \dontrun{
#' gather_ntests("results/", "ntests\\.txt", "run1")
#' }
#'
#' @export
gather_ntests <- function(directory, pattern = "ntests\\.txt", prefix = "defaultprefix") {
  # Treat empty strings as "use default"
  if (is.null(pattern) || identical(pattern, "")) pattern <- "ntests\\.txt"
  if (is.null(prefix)  || identical(prefix,  "")) prefix  <- "defaultprefix"

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
    # ensure we only keep first three columns (all, trans, cis) and coerce numeric
    data <- data[, seq_len(min(3, ncol(data))), drop = FALSE]
    data[] <- lapply(data, as.numeric)

    data_list[[length(data_list) + 1]] <- data
  }

  combined_data <- dplyr::bind_rows(data_list)

  if (anyNA(combined_data)) {
    warning("NA values were present in the data and were removed")
  }

  summed_data <- colSums(combined_data, na.rm = TRUE)
  # pad to length 3 in case some files had fewer columns
  if (length(summed_data) < 3) summed_data <- c(summed_data, rep(0, 3 - length(summed_data)))

  result <- data.frame(all = summed_data[1], trans = summed_data[2], cis = summed_data[3])

  out_path <- file.path(directory, paste(prefix, "ntests_combined.txt", sep = "_"))
  write.table(result, file = out_path, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)

  return(result)
}

# Main function to handle arguments and call the gather_ntests function
#' @keywords internal
main <- function() {
  args <- parse_args()
  result <- gather_ntests(args$directory, args$pattern, args$prefix)
  print(result)
}

# Only run main() if the script is executed directly
if (sys.nframe() == 0) {
  suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
  })
  main()
}
