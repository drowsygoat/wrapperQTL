# Function to parse command-line arguments using argparse
parse_args <- function() {
  parser <- argparse::ArgumentParser(description = "Process eQTL files")
  parser$add_argument("directory", help = "Directory containing the files")
  parser$add_argument("pattern", help = "Pattern to match files")
  parser$add_argument("prefix", help = "Prefix to filter files")
  
  args <- parser$parse_args()
  return(args)
}

#' Gather and Summarize eQTL Test Counts
#'
#' Reads multiple "ntests" files from a directory, filters by pattern and prefix,
#' and outputs the sum of test counts (all, trans, cis).
#'
#' @param directory A character string specifying the path to the directory containing files.
#' @param pattern A regex pattern to match files (e.g., "ntests").
#' @param prefix A string to further filter matched files.
#'
#' @return A data.frame with one row and columns: `all`, `trans`, `cis`.
#' Also writes the result to a file named `{prefix}_ntests_combined.txt` in the given directory.
#' 
#' @examples
#' \dontrun{
#' gather_ntests("results/", "ntests", "run1")
#' }
#'
#' @export
gather_ntests <- function(directory, pattern, prefix) {
  ...
}
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
    data_list <- append(data_list, list(data))
  }

  combined_data <- dplyr::bind_rows(data_list)

  if (anyNA(combined_data)) {
    warning("NA values were present in the data and were removed") 
  }

  summed_data <- colSums(combined_data, na.rm = TRUE)
  result <- data.frame(all = summed_data[1], trans = summed_data[2], cis = summed_data[3])

  write.table(result,
              file = file.path(directory, paste(prefix, "ntests_combined.txt", sep = "_")),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

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
