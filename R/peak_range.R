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
#' @param fdr_col Character scalar. Name of the column in \code{qtl_data}
#'   containing the (adjusted) p-values used for peak finding (default \code{"FDR"}).
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
#' @param SEARCHDISTANCE Numeric (bp). For \code{"DistanceFromPeak"}, the maximum
#'   absolute distance from the peak position. For \code{"AdjacentPoints"}, the
#'   maximum gap allowed between consecutive significant positions while extending
#'   upstream/downstream. Default \code{3e6}.
#' @param LOCALPEAKTHRESHOLD Numeric > 0. The dynamic threshold factor relative to the
#'   local peak p-value; a value of \code{0.01} means points are considered locally
#'   significant if \code{p <= PEAKVALUE / 0.01}. Default \code{0.01}.
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
#' @export
peak_range <- function(
  qtl_data,                     # overall_results_table
  fdr_col = "FDR",              # NEW: which column to use for (adjusted) p-values
  METHOD = "DistanceFromPeak",  # DistanceFromPeak or AdjacentPoints
  SIGTHRESHOLD = 0.01,          # Stop if next peak's p-value > SIGTHRESHOLD
  SEARCHDISTANCE = 3e+6,        # Max distance from peak (or max gap between adjacent points)
  LOCALPEAKTHRESHOLD = 0.01,    # Local threshold factor relative to peak p
  MINGAPSIZE = 0                # Merge ranges with gaps <= MINGAPSIZE
) {
  library(tidyverse)
  library(data.table)

  # Pull required vectors
  CHRVECTOR    <- qtl_data$chr_snp
  POSVECTOR    <- qtl_data$location

  if (!fdr_col %in% names(qtl_data))
    stop("Column '", fdr_col, "' not found in qtl_data.")
  PVALUEVECTOR <- qtl_data[[fdr_col]]
  if (!is.numeric(PVALUEVECTOR))
    stop("Column '", fdr_col, "' must be numeric.")

  if (missing(CHRVECTOR)) {
    CHRVECTOR <- rep("ChromosomeC", length(POSVECTOR))
  }

  INPUTDATA <- data.frame(
    Chromosome = CHRVECTOR,
    Position   = POSVECTOR,
    Pvalue     = PVALUEVECTOR
  )

  # Create output data frame.
  OUTPUT <- data.frame(
    Chromosome   = NA,
    PeakPosition = NA,
    PeakPvalue   = NA,
    RegionStart  = NA,
    RegionEnd    = NA,
    NPeaks       = NA,
    NSignificant = NA
  )

  for (C in unique(INPUTDATA$Chromosome)) {

    DATA_CHRC <- INPUTDATA %>%
      dplyr::filter(Chromosome == C) %>%
      dplyr::arrange(Position)

    OUTPUT_CHRC <- data.frame(
      Chromosome   = NA,
      PeakPosition = NA,
      PeakPvalue   = NA,
      RegionStart  = NA,
      RegionEnd    = NA
    )

    while (nrow(DATA_CHRC) > 0) {
      PEAKVALUE <- min(DATA_CHRC$Pvalue)
      if (PEAKVALUE > SIGTHRESHOLD) break

      X <- PEAKVALUE / LOCALPEAKTHRESHOLD
      PEAKINDICES <- which(DATA_CHRC$Pvalue == PEAKVALUE)

      if (METHOD == "DistanceFromPeak") {
        RANGES_CHRC <- lapply(PEAKINDICES, function(idx) {
          PEAKPOSITION <- DATA_CHRC$Position[idx]
          PEAKPVALUE   <- DATA_CHRC$Pvalue[idx]
          DATASUBSET <- DATA_CHRC %>%
            dplyr::filter(
              Position >= (PEAKPOSITION - SEARCHDISTANCE),
              Position <= (PEAKPOSITION + SEARCHDISTANCE)
            )
          LEFTLIMIT  <- min(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
          RIGHTLIMIT <- max(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
          c(PEAKPOSITION, PEAKPVALUE, LEFTLIMIT, RIGHTLIMIT)
        }) %>% do.call(rbind, .) %>% as.data.frame() %>%
          dplyr::mutate(Chromosome = C)
        colnames(RANGES_CHRC) <- c("PeakPosition","PeakPvalue","RegionStart","RegionEnd","Chromosome")
      }

      if (METHOD == "AdjacentPoints") {
        RANGES_CHRC <- lapply(PEAKINDICES, function(idx) {
          PEAKPOSITION <- DATA_CHRC$Position[idx]
          PEAKPVALUE   <- DATA_CHRC$Pvalue[idx]

          CONTINUE  <- TRUE
          LEFTLIMIT <- PEAKPOSITION
          while (CONTINUE) {
            DATASUBSET <- DATA_CHRC %>%
              dplyr::filter(
                Position >= (LEFTLIMIT - SEARCHDISTANCE),
                Position < LEFTLIMIT
              )
            if (nrow(DATASUBSET) > 0 && any(DATASUBSET$Pvalue <= X)) {
              LEFTLIMIT <- min(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
            } else CONTINUE <- FALSE
          }

          CONTINUE   <- TRUE
          RIGHTLIMIT <- PEAKPOSITION
          while (CONTINUE) {
            DATASUBSET <- DATA_CHRC %>%
              dplyr::filter(
                Position > RIGHTLIMIT,
                Position <= (RIGHTLIMIT + SEARCHDISTANCE)
              )
            if (nrow(DATASUBSET) > 0 && any(DATASUBSET$Pvalue <= X)) {
              RIGHTLIMIT <- max(DATASUBSET$Position[DATASUBSET$Pvalue <= X])
            } else CONTINUE <- FALSE
          }

          c(PEAKPOSITION, PEAKPVALUE, LEFTLIMIT, RIGHTLIMIT)
        }) %>% do.call(rbind, .) %>% as.data.frame() %>%
          dplyr::mutate(Chromosome = C)
        colnames(RANGES_CHRC) <- c("PeakPosition","PeakPvalue","RegionStart","RegionEnd","Chromosome")
      }

      OUTPUT_CHRC <- rbind(OUTPUT_CHRC, RANGES_CHRC)

      for (i in 1:nrow(RANGES_CHRC)) {
        DATA_CHRC <- DATA_CHRC %>%
          dplyr::filter(!(Position >= RANGES_CHRC$RegionStart[i] &
                          Position <= RANGES_CHRC$RegionEnd[i]))
      }
    }

    if (nrow(OUTPUT_CHRC[!is.na(OUTPUT_CHRC$PeakPosition), ]) > 0) {
      OUTPUT_CHRC <- OUTPUT_CHRC %>%
        dplyr::filter(!is.na(Chromosome)) %>%
        dplyr::arrange(RegionStart) %>%
        dplyr::group_by(
          Index = cumsum(
            cummax(dplyr::lag(RegionEnd + MINGAPSIZE / 2,
                              default = data.table::first(RegionEnd + MINGAPSIZE / 2))) <
              (RegionStart - MINGAPSIZE / 2)
          )
        ) %>%
        dplyr::summarise(
          RegionStart  = min(RegionStart),
          RegionEnd    = max(RegionEnd),
          PeakPosition = PeakPosition[PeakPvalue == min(PeakPvalue)],
          PeakPvalue   = min(PeakPvalue),
          NPeaks       = dplyr::n()
        ) %>%
        dplyr::mutate(Chromosome = C) %>%
        as.data.frame()

      OUTPUT_CHRC <- OUTPUT_CHRC[, c("Chromosome","PeakPosition","PeakPvalue","RegionStart","RegionEnd","NPeaks")]

      DATA_CHRC <- INPUTDATA %>%
        dplyr::filter(Chromosome == C, Pvalue <= SIGTHRESHOLD)

      OUTPUT_CHRC$NSignificant <- NA_integer_
      for (i in 1:nrow(OUTPUT_CHRC)) {
        OUTPUT_CHRC$NSignificant[i] <- DATA_CHRC %>%
          dplyr::filter(Position >= OUTPUT_CHRC$RegionStart[i],
                        Position <= OUTPUT_CHRC$RegionEnd[i]) %>%
          nrow()
      }

      OUTPUT <- rbind(OUTPUT, OUTPUT_CHRC)
    }
  }

  OUTPUT <- OUTPUT %>%
    dplyr::filter(!is.na(PeakPosition)) %>%
    dplyr::arrange(Chromosome, PeakPvalue)

  OUTPUT$RegionLength <- OUTPUT$RegionEnd - OUTPUT$RegionStart + 1

  return(OUTPUT)
}