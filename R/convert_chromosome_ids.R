#' Convert chromosome names by adding or removing "chr" prefix
#'
#' This function converts chromosome identifiers by either adding the "chr" prefix (e.g., `1` → `chr1`) 
#' or removing it (e.g., `chr1` → `1`). It also handles the mitochondrial chromosome (`MT` ↔ `chrM`).
#'
#' @param chroms A character vector of chromosome names.
#' @param direction A string specifying the conversion direction: `"add"` to add the "chr" prefix, or `"remove"` to remove it.
#'
#' @return A character vector of converted chromosome names.
#'
#' @examples
#' convert_chromosome_ids(c("1", "2", "MT"), direction = "add")
#' convert_chromosome_ids(c("chr1", "chr2", "chrM"), direction = "remove")
#'
#' @export
convert_chromosome_ids <- function(chroms, direction = c("add", "remove")) {
  direction <- match.arg(direction)
  
  if (direction == "add") {
    return(sapply(chroms, function(chr) {
      chr <- toupper(as.character(chr))
      if (chr == "MT") {
        return("chrM")
      } else {
        return(paste0("chr", chr))
      }
    }))
  } else if (direction == "remove") {
    return(sapply(chroms, function(chr) {
      chr <- as.character(chr)
      chr <- sub("^chr", "", chr, ignore.case = TRUE)
      if (toupper(chr) == "M") {
        return("MT")
      } else {
        return(chr)
      }
    }))
  }
}
