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
