CleanIDs <- function(input, Capitalize = TRUE, pattern = NULL) {
  library(data.table)

  clean_vector <- function(ids) {
    ids <- as.character(ids)
    parts <- strsplit(ids, "_")

    cleaned <- vapply(parts, function(x) {
      if (!is.null(pattern) && all(pattern <= length(x))) {
        x_keep <- x[pattern]
      } else {
        # Check for repeated halves: A_B_A_B
        n <- length(x)
        if (n %% 2 == 0 && all(x[1:(n / 2)] == x[(n / 2 + 1):n])) {
          x_keep <- x[1:(n / 2)]
        } else {
          x_keep <- x
        }
      }
      paste0(x_keep, collapse = "")
    }, character(1))

    if (Capitalize) {
      cleaned <- toupper(cleaned)
    } else {
      cleaned <- tolower(cleaned)
    }

    return(cleaned)
  }

  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    # File input
    df <- fread(input, header = TRUE)
    original_ids <- colnames(df)
    first_id <- original_ids[1]
    sample_ids <- original_ids[-1]

    cleaned_ids <- clean_vector(sample_ids)
    setnames(df, c(first_id, cleaned_ids))

    out_path <- file.path(
      dirname(input),
      paste0(tools::file_path_sans_ext(basename(input)), "_names_repaired.txt")
    )
    fwrite(df, out_path, sep = "\t", quote = FALSE)
    message("Cleaned file saved to: ", out_path)
    return(invisible(out_path))
  } else if (is.vector(input)) {
    # Vector input
    return(clean_vector(input))
  } else {
    stop("Input must be a character vector or a valid file path.")
  }
}
