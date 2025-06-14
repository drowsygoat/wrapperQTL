% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_chromosome_ids.R
\name{convert_chromosome_ids}
\alias{convert_chromosome_ids}
\title{Convert chromosome names by adding or removing "chr" prefix}
\usage{
convert_chromosome_ids(chroms, direction = c("add", "remove"))
}
\arguments{
\item{chroms}{A character vector of chromosome names.}

\item{direction}{A string specifying the conversion direction: \code{"add"} to add the "chr" prefix, or \code{"remove"} to remove it.}
}
\value{
A character vector of converted chromosome names.
}
\description{
This function converts chromosome identifiers by either adding the "chr" prefix (e.g., \code{1} → \code{chr1})
or removing it (e.g., \code{chr1} → \code{1}). It also handles the mitochondrial chromosome (\code{MT} ↔ \code{chrM}).
}
\examples{
convert_chromosome_ids(c("1", "2", "MT"), direction = "add")
convert_chromosome_ids(c("chr1", "chr2", "chrM"), direction = "remove")

}
