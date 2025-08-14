

#' Clean and Standardize Sample IDs
#'
#' Cleans sample IDs by optionally selecting parts of IDs, removing duplicate halves, 
#' and converting to uppercase or lowercase. Can operate directly on a character vector 
#' or on the header line of a tab-delimited file (saving a repaired version).
#'
#' @param input A character vector of IDs, or a file path to a tab-delimited text file 
#'   where the first row contains IDs (first column is assumed to be a non-ID column name).
#' @param Capitalize Logical; if `TRUE`, converts IDs to uppercase, otherwise to lowercase. 
#'   Default is `TRUE`.
#' @param pattern Optional integer vector specifying which underscore-separated parts of 
#'   the ID to keep (by position). If `NULL` (default), keeps all parts, but if the ID 
#'   consists of two identical halves, only keeps the first half.
#'
#' @return 
#' - If `input` is a vector: a cleaned character vector.
#' - If `input` is a file path: (invisibly) returns the path to the repaired file.
#'
#' @details
#' If a file path is provided, the function:
#' 1. Reads the file.
#' 2. Cleans the sample IDs in the header.
#' 3. Writes a new file with `_names_repaired.txt` appended to the name.
#'
#' When cleaning IDs:
#' - IDs are split on underscores.
#' - If `pattern` is specified, only the specified parts are kept.
#' - If `pattern` is `NULL` and the ID has two identical halves, only the first half is kept.
#' - Otherwise, all parts are concatenated.
#'
#' @examples
#' # Clean a vector of IDs
#' CleanIDs(c("S1_A", "S2_A"))
#'
#' # Keep only the first two parts of IDs
#' CleanIDs(c("S1_part1_part2_extra", "S2_part1_part2_extra"), pattern = 1:2)
#'
#' # Clean IDs in a file and save repaired file
#' \dontrun{
#' CleanIDs("matrix_eqtl_input/SNP.txt")
#' }
#'
#' @export
CleanIDs <- function(input, Capitalize = TRUE, pattern = NULL) {
  clean_vector <- function(ids) {
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

  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    lines <- readLines(input)
    header <- strsplit(lines[[1]], "\t")[[1]]
    first_id <- header[1]
    sample_ids <- header[-1]

    cleaned_ids <- clean_vector(sample_ids)
    new_header <- c(first_id, cleaned_ids)

    out_path <- file.path(
      dirname(input),
      paste0(tools::file_path_sans_ext(basename(input)), "_names_repaired.txt")
    )

    data_lines <- lines[-1]
    writeLines(
      c(paste(new_header, collapse = "\t"),
        data_lines),
      out_path
    )
    message("Cleaned file saved to: ", out_path)
    return(invisible(out_path))
  } else if (is.vector(input)) {
    return(clean_vector(input))
  } else {
    stop("Input must be a character vector or a valid file path.")
  }
}
