#' Extract First Regex Match from String
#'
#' A safer wrapper around `regexpr()` + `regmatches()` to extract the first pattern match from each string.
#'
#' @param strings Character vector to search in.
#' @param pattern Regular expression pattern to match.
#'
#' @return Character vector of matched substrings. If no match is found, returns NA.
#' @export
grep_o <- function(strings, pattern) {
  matches <- regexpr(pattern, strings, perl = TRUE)
  matched <- regmatches(strings, matches)
  matched[matches == -1] <- NA_character_
  matched
}