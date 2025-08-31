#' Derive Significant QTL Ranges Around Local Peaks
#'
#' Iteratively identifies local peaks (smallest p-values) per chromosome and
#' grows significant ranges around each peak using one of two strategies:
#' \emph{DistanceFromPeak} (bounded by absolute genomic distance from the peak)
#' or \emph{AdjacentPoints} (bounded by the largest allowed gap between adjacent
#' significant positions). After collecting candidate ranges, overlapping/nearby
#' ranges are merged (controlled by \code{MINGAPSIZE}).
#'
#' @param qtl_data A data frame/data.table/tibble with at least the columns:
#'   \describe{
#'     \item{\code{chr_snp}}{Chromosome identifier (character/factor).}
#'     \item{\code{location}}{Genomic position (integer/numeric, increasing within chromosome).}
#'     \item{\code{FDR}}{Adjusted p-value (numeric); smaller is more significant.}
#'   }
#'   Additional columns are ignored.
#'
#' @param METHOD Character scalar. Range-building strategy:
#'   \itemize{
#'     \item \code{"DistanceFromPeak"}: Define the region as the farthest significant
#'       points within \code{SEARCHDISTANCE} bp upstream/downstream of the peak.
#'       Significance within the region is evaluated vs. a local dynamic threshold
#'       (\code{PEAKVALUE / LOCALPEAKTHRESHOLD}).
#'     \item \code{"AdjacentPoints"}: Walk upstream and downstream from the peak,
#'       extending the region while adjacent significant points are no farther apart
#'       than \code{SEARCHDISTANCE} bp. Uses the same local dynamic threshold for significance.
#'   }
#'   Default \code{"DistanceFromPeak"}.
#'
#' @param SIGTHRESHOLD Numeric in (0,1]. Global significance cutoff for peaks.
#'   Iterations stop when the current best (minimum) p-value exceeds this threshold.
#'   Default \code{0.01}.
#'
#' @param SEARCHDISTANCE Numeric (bp). For \code{"DistanceFromPeak"}, the maximum
#'   absolute distance from the peak position. For \code{"AdjacentPoints"}, the
#'   maximum gap allowed between consecutive significant positions while extending
#'   upstream/downstream. Default \code{3e6}.
#'
#' @param LOCALPEAKTHRESHOLD Numeric > 0. The dynamic threshold factor relative to the
#'   local peak p-value; a value of \code{0.01} means points are considered locally
#'   significant if \code{p <= PEAKVALUE / 0.01} (i.e., up to two orders of magnitude
#'   less significant than the peak). Default \code{0.01}.
#'
#' @param MINGAPSIZE Non-negative numeric (bp). Adjacent/overlapping ranges with a gap
#'   \eqn{\le} \code{MINGAPSIZE} are merged. Default \code{0}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{Chromosome}}{Chromosome id.}
#'     \item{\code{PeakPosition}}{Position of local peak.}
#'     \item{\code{PeakPvalue}}{Local peak p-value.}
#'     \item{\code{RegionStart}}{Start coordinate (inclusive).}
#'     \item{\code{RegionEnd}}{End coordinate (inclusive).}
#'     \item{\code{NPeaks}}{Number of peaks merged into the range.}
#'     \item{\code{NSignificant}}{Count of positions with \code{p <= SIGTHRESHOLD} inside the range.}
#'     \item{\code{RegionLength}}{Computed as \code{RegionEnd - RegionStart + 1}.}
#'   }
#'
#' @details
#' The algorithm loops per chromosome: pick the minimum p-value (the peak);
#' construct a candidate range around it using \code{METHOD} and local dynamic
#' threshold; remove covered rows; repeat until no positions with
#' \code{p <= SIGTHRESHOLD} remain. Finally, overlapping/nearby ranges are merged
#' with tolerance \code{MINGAPSIZE}.
#'
#' @section Assumptions:
#' \itemize{
#'   \item Positions are on the same coordinate system and comparable within a chromosome.
#'   \item For reproducibility and performance, positions within each chromosome should be sorted ascending.
#'   \item \code{qtl_data$FDR} contains the p-values to be used (adjusted). If you wish to use raw p-values,
#'         pass them in the same column or rename accordingly.
#' }
#'
#' @examples
#' \dontrun{
#' out <- peak_range(qtl_data, METHOD = "DistanceFromPeak",
#'                   SIGTHRESHOLD = 0.01, SEARCHDISTANCE = 1e6,
#'                   LOCALPEAKTHRESHOLD = 0.01, MINGAPSIZE = 5e3)
#' head(out)
#' }
#'
#' @export
peak_range <- function(
  qtl_data,                 # overall_results_table
  METHOD = "DistanceFromPeak",  # Method used to obtain ranges: DistanceFromPeak, AdjacentPoints.
  SIGTHRESHOLD = 0.01,          # P-value significance threshold. Stop searching for peaks if the p-value of the next peak is above SIGTHRESHOLD.
  SEARCHDISTANCE = 3e+6,        # Maximum distance from local peak (if DistanceFromPeak method), or maximum distance between adjacent significant points (if AdjacentPoints method).
  LOCALPEAKTHRESHOLD = 0.01,    # Relative threshold of local peak (0.01 represents 2 orders of magnitude).
  MINGAPSIZE = 0                # Minimum size of gaps between ranges. Adjacent ranges with a gap <= MINGAPSIZE will be merged.
) {
  # 1. Find peak position (lowest p-value).
  # 2. Set local peak p-value threshold (e.g., 0.01 = two orders of magnitude below local peak) to limit local peak height.
  # 3. DistanceFromPeak: Find the farthest significant points (below local peak p-value threshold) within the defined
  #    distance from local peak. This method limits the extension of the region, without limiting the gap between adjacent significant points.
  # 3. AdjacentPoints: Find the farthest significant points (below local peak p-value threshold) upstream/downstream
  #    the last picked point upstream/downstream the local peak. This method limits the gap between adjacent significant points, not
  #    the extension of the total regions.
  # 4. Save peak and range.
  # 5. Remove range from data.
  # 6. Go back to step 1.

  library(tidyverse)
  library(data.table)

  CHRVECTOR <- qtl_data$chr_snp   # Vector of chromosomes.
  POSVECTOR <- qtl_data$location  # Vector of chromosomal positions.
  PVALUEVECTOR <- qtl_data$FDR    # Vector of p-values (FDR).

  if (missing(CHRVECTOR)) {
    CHRVECTOR <- rep("ChromosomeC", length(POSVECTOR))
  }

  INPUTDATA <- data.frame(
    Chromosome = CHRVECTOR,
    Position = POSVECTOR,
    Pvalue = PVALUEVECTOR
  )

  # Create output data frame.
  OUTPUT <- data.frame(
    Chromosome = NA,
    PeakPosition = NA,
    PeakPvalue = NA,
    RegionStart = NA,
    RegionEnd = NA,
    NPeaks = NA,
    NSignificant = NA
  )

  for (C in unique(INPUTDATA$Chromosome)) {
    # Get data for chromosome C.
    DATA_CHRC <- INPUTDATA %>%
      dplyr::filter(Chromosome == C) %>%
      dplyr::arrange(Position)

    # Create output data frame.
    OUTPUT_CHRC <- data.frame(
      Chromosome = NA,
      PeakPosition = NA,
      PeakPvalue = NA,
      RegionStart = NA,
      RegionEnd = NA
    )

    # Continue searching while there are significant positions that haven't been assigned to a range.
    while (nrow(DATA_CHRC) > 0) {
      # Find peak value.
      PEAKVALUE <- min(DATA_CHRC$Pvalue)

      # No meaningful peaks left.
      if (PEAKVALUE > SIGTHRESHOLD) break

      # Local peak dynamic threshold.
      X <- PEAKVALUE / LOCALPEAKTHRESHOLD
      # if (X > SIGTHRESHOLD) { X <- SIGTHRESHOLD }

      # Find all peak indices.
      PEAKINDICES <- which(DATA_CHRC$Pvalue == PEAKVALUE)

      if (METHOD == "DistanceFromPeak") {
        # Find ranges around each peak.
        RANGES_CHRC <- lapply(PEAKINDICES, function(idx) {
          PEAKPOSITION <- DATA_CHRC$Position[idx]
          PEAKPVALUE <- DATA_CHRC$Pvalue[idx]
          DATASUBSET <- DATA_CHRC %>%
            dplyr::filter(
              Position >= (PEAKPOSITION - SEARCHDISTANCE),
              Position <= (PEAKPOSITION + SEARCHDISTANCE)
            )

          # Search upstream peak.
          LEFTLIMIT <- min(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
          # Search downstream peak.
          RIGHTLIMIT <- max(DATASUBSET$Position[DATASUBSET$Pvalue <= X])

          return(c(PEAKPOSITION, PEAKPVALUE, LEFTLIMIT, RIGHTLIMIT))
        })

        RANGES_CHRC <- do.call(rbind, RANGES_CHRC)
        RANGES_CHRC <- as.data.frame(RANGES_CHRC) %>%
          dplyr::mutate(Chromosome = C)
        colnames(RANGES_CHRC) <- c("PeakPosition", "PeakPvalue", "RegionStart", "RegionEnd", "Chromosome")
      }

      if (METHOD == "AdjacentPoints") {
        # Find ranges around each peak.
        RANGES_CHRC <- lapply(PEAKINDICES, function(idx) {
          PEAKPOSITION <- DATA_CHRC$Position[idx]
          PEAKPVALUE <- DATA_CHRC$Pvalue[idx]
          DATASUBSET <- DATA_CHRC %>%
            dplyr::filter(
              Position >= (PEAKPOSITION - SEARCHDISTANCE),
              Position <= (PEAKPOSITION + SEARCHDISTANCE)
            )

          # Search upstream peak.
          CONTINUE <- TRUE
          LEFTLIMIT <- PEAKPOSITION
          while (CONTINUE == TRUE) {
            DATASUBSET <- DATA_CHRC %>%
              dplyr::filter(
                Position >= (LEFTLIMIT - SEARCHDISTANCE),
                Position < LEFTLIMIT
              )
            if (nrow(DATASUBSET) > 0) {
              if (sum(DATASUBSET$Pvalue <= X) > 0) {
                LEFTLIMIT <- min(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
              } else {
                CONTINUE <- FALSE
              }
            } else {
              CONTINUE <- FALSE
            }
          }

          # Search downstream peak.
          CONTINUE <- TRUE
          RIGHTLIMIT <- PEAKPOSITION
          while (CONTINUE == TRUE) {
            DATASUBSET <- DATA_CHRC %>%
              dplyr::filter(
                Position > RIGHTLIMIT,
                Position <= (RIGHTLIMIT + SEARCHDISTANCE)
              )
            if (nrow(DATASUBSET) > 0) {
              if (sum(DATASUBSET$Pvalue <= X) > 0) {
                RIGHTLIMIT <- max(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
              } else {
                CONTINUE <- FALSE
              }
            } else {
              CONTINUE <- FALSE
            }
          }

          return(c(PEAKPOSITION, PEAKPVALUE, LEFTLIMIT, RIGHTLIMIT))
        })

        RANGES_CHRC <- do.call(rbind, RANGES_CHRC)
        RANGES_CHRC <- as.data.frame(RANGES_CHRC) %>%
          dplyr::mutate(Chromosome = C)
        colnames(RANGES_CHRC) <- c("PeakPosition", "PeakPvalue", "RegionStart", "RegionEnd", "Chromosome")
      }

      # Save ranges.
      OUTPUT_CHRC <- rbind(OUTPUT_CHRC, RANGES_CHRC)

      # Remove ranges from the data.
      for (i in 1:nrow(RANGES_CHRC)) {
        DATA_CHRC <- DATA_CHRC %>%
          dplyr::filter(!(Position >= RANGES_CHRC$RegionStart[i] & Position <= RANGES_CHRC$RegionEnd[i]))
      }
    }

    # Merge overlapping regions.
    # If several peaks with the same p-value exist within a range, repeat range for each peak.
    if (nrow(OUTPUT_CHRC[!is.na(OUTPUT_CHRC$PeakPosition), ]) > 0) {
      OUTPUT_CHRC <- OUTPUT_CHRC %>%
        dplyr::filter(!is.na(Chromosome)) %>%
        dplyr::arrange(RegionStart) %>%
        dplyr::group_by(
          Index = cumsum(
            cummax(dplyr::lag(RegionEnd + MINGAPSIZE / 2, default = data.table::first(RegionEnd + MINGAPSIZE / 2))) <
              (RegionStart - MINGAPSIZE / 2)
          )
        ) %>%
        dplyr::summarise(
          RegionStart = min(RegionStart),
          RegionEnd = max(RegionEnd),
          PeakPosition = PeakPosition[PeakPvalue == min(PeakPvalue)],
          PeakPvalue = min(PeakPvalue),
          NPeaks = dplyr::n()
        ) %>%
        dplyr::mutate(Chromosome = C) %>%
        as.data.frame()

      OUTPUT_CHRC <- OUTPUT_CHRC[, c("Chromosome", "PeakPosition", "PeakPvalue", "RegionStart", "RegionEnd", "NPeaks")]

      # Get number of significant points per region.
      DATA_CHRC <- INPUTDATA %>%
        dplyr::filter(
          Chromosome == C,
          Pvalue <= SIGTHRESHOLD
        )

      OUTPUT_CHRC$NSignificant <- NA_integer_
      for (i in 1:nrow(OUTPUT_CHRC)) {
        OUTPUT_CHRC$NSignificant[i] <- DATA_CHRC %>%
          dplyr::filter(
            Position >= OUTPUT_CHRC$RegionStart[i],
            Position <= OUTPUT_CHRC$RegionEnd[i]
          ) %>%
          nrow()
      }

      # Save ranges.
      OUTPUT <- rbind(OUTPUT, OUTPUT_CHRC)
    }
  }

  # Remove NA entries.
  OUTPUT <- OUTPUT %>%
    dplyr::filter(!is.na(PeakPosition)) %>%
    dplyr::arrange(Chromosome, PeakPvalue)

  OUTPUT$RegionLength <- OUTPUT$RegionEnd - OUTPUT$RegionStart + 1

  return(OUTPUT)
}
