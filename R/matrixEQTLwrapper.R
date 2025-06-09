# matrixEQTLwrapper_with_logging.R

library(MatrixEQTL)
library(futile.logger)
library(R.utils)

#' Debug checkpoint logger
checkpoint <- function(msg, verbose = TRUE) {
  if (verbose) flog.info(paste0("[Checkpoint] ", msg), name = "matrixeqtl")
}

#' Initialize Logging
init_logging <- function(logfile = NULL, verbose = TRUE) {
  flog.threshold(if (verbose) INFO else WARN)
  if (!is.null(logfile)) {
    flog.appender(appender.file(logfile), name = "matrixeqtl")
  } else {
    flog.appender(appender.console(), name = "matrixeqtl")
  }
}

#' Get model constant from name
get_model_constant <- function(model) {
  switch(tolower(model),
         "linear" = modelLINEAR,
         "anova"  = modelANOVA,
         "cross"  = modelLINEAR_CROSS,
         stop("useModel must be one of 'linear', 'anova', or 'cross'")
  )
}

#' Load a SlicedData object from file
load_sliced_file <- function(file, sliceSize = 2000, first_col_is_rownames = TRUE) {
  sd <- SlicedData$new()
  sd$fileSliceSize <- sliceSize
  
  if (grepl("\\.rds$", file, ignore.case = TRUE)) {
    data <- readRDS(file)
    message("Loaded RDS preview:")
    print(head(data))

    if (is.matrix(data) || is.data.frame(data)) {
      if (first_col_is_rownames) {
        if (ncol(data) < 2) stop("RDS data has too few columns to extract rownames.")
        rownames(data) <- data[[1]]
        data <- data[, -1, drop = FALSE]
      }
      if (is.null(colnames(data))) stop("Data must have column (sample) names.")
          print(head(data))

      sd$CreateFromMatrix(as.matrix(data))
    } else {
      stop("RDS file must contain a matrix or a data.frame.")
    }

  } else {
    sd$fileDelimiter <- "\t"
    sd$fileOmitCharacters <- "NA"
    sd$fileSkipRows <- 1
    sd$fileSkipColumns <- 1
    sd$LoadFile(file)
  }

  return(sd)
}

#' Validate and reorder sample columns if needed, based on trimmed sample names
#' Trims suffixes from sample names (e.g., removes '_C1_RNA') before checking
#' @export
validate_sample_overlap <- function(feature_data, snps, covs, group_name, verbose = TRUE) {

  # Extract column names
  expr_samples_raw <- feature_data$columnNames
  snp_samples_raw <- snps$columnNames
  cov_samples_raw <- covs$columnNames

  # Trim everything after the first underscore
  trim_names <- function(x) sub("_.*$", "", x)

  expr_samples <- trim_names(expr_samples_raw)
  snp_samples <- trim_names(snp_samples_raw)
  cov_samples <- trim_names(cov_samples_raw)

  # Check for missing samples
  missing_in_snps <- setdiff(expr_samples, snp_samples)
  missing_in_covs <- setdiff(expr_samples, cov_samples)

  if (length(missing_in_snps) > 0) {
    stop(sprintf("Missing samples in SNP data: %s", paste(missing_in_snps, collapse = ", ")))
  }

  if (length(missing_in_covs) > 0) {
    stop(sprintf("Missing samples in covariates: %s", paste(missing_in_covs, collapse = ", ")))
  }

  # Reorder SNPs and covariates to match feature data
  reorder_indices_snps <- match(expr_samples, snp_samples)
  reorder_indices_covs <- match(expr_samples, cov_samples)

  if (verbose) flog.info("Reordering SNP and covariate data to match expression sample order", name = "matrixeqtl")

  snps$ColumnSubsample(reorder_indices_snps)
  covs$ColumnSubsample(reorder_indices_covs)

  # Update sample names in SlicedData objects to trimmed IDs
  snps$columnNames <- expr_samples
  covs$columnNames <- expr_samples
  feature_data$columnNames <- expr_samples
}

# #' Validate and reorder sample columns if needed
# #' @export
# validate_sample_overlap <- function(feature_data, snps, covs, group_name, verbose = TRUE) {
  
#   if (covs$nRows() > feature_data$nCols()) stop("More covariates than samples.")
#   if (length(feature_data$columnNames) > length(snps$columnNames)) stop("More samples in features than in SNPs.")

#   if (!check_complete_match(snps$columnNames, feature_data$columnNames)) {

#     if (verbose) flog.info("Reordering SNP and covariate data to match features", name = "matrixeqtl")

#     idx <- sort(which(sapply(snps$columnNames, check_match, group_name_vector = feature_data$columnNames)))

#     snps$ColumnSubsample(idx)
#     covs$ColumnSubsample(idx)
#   }
# }

#' Run MatrixEQTL
#' @export
matrixEQTLwrapper <- function(
  feature_locations_path,
  feature_data_path,
  snpFilePath,
  covFilePath,
  snpLocPath,
  group_name,
  resultsDir = getwd(),
  cisDist = 1e6,
  pvOutputThreshold = 1e-5,
  pvOutputThresholdCis = 1e-4,
  useModel = "linear",
  minPvByGeneSnp = TRUE,
  noFDRsaveMemory = FALSE,
  # SNPsInChunks = NULL,
  prefix = NULL,
  pvalueHist = NULL,
  verbose = TRUE,
  dry_run = TRUE 
) {
  checkpoint("Starting matrixEQTLwrapper", verbose)

  chunk_snp <- grep_o(snpFilePath, "chunk_[0-9]+")
  chunk_feature <- grep_o(feature_data_path, "chunk_[0-9]+")
  useModelConst <- get_model_constant(useModel)

  prefix <- if (is.null(prefix)) {
    paste(useModel, cisDist, chunk_snp, chunk_feature, sep = "_")
  } else {
    paste(useModel, cisDist, chunk_snp, chunk_feature, prefix, sep = "_")
  }

  checkpoint(paste("Prefix set to:", prefix), verbose)

  subdir <- file.path(resultsDir, paste(group_name, "MEQTL_res", sep = "_"))
    if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
    }
  base_prefix <- file.path(subdir, paste(group_name, tolower(prefix), sep = "_"))

  output_file_name_trans <- paste0(base_prefix, "_eQTL_trans.txt")
  output_file_name_cis   <- paste0(base_prefix, "_eQTL_cis.txt")
  gz_trans <- paste0(output_file_name_trans, ".gz")
  gz_cis   <- paste0(output_file_name_cis, ".gz")
  result_path <- paste0(base_prefix, "_results_MEQTL.rds")
  log_path <- paste0(base_prefix, "_MEQTL_runtime_log.txt")
  ntests_path <- paste0(base_prefix, "_ntests.txt")

  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  checkpoint(paste("Ensured log directory:", dirname(log_path)), verbose)

  init_logging(log_path, verbose)
  checkpoint(paste("Logging initialized to", log_path), verbose)

  if (dry_run) {
    flog.info("[DRY RUN] Would run eQTL with SNP: %s, Expr: %s", snpFilePath, feature_data_path, name = "matrixeqtl")
    return(invisible(NULL))
  }

  if (file.exists(gz_trans) && file.exists(gz_cis)) {
    warning(sprintf("Result files already exist: %s and/or %s. Skipping run.", gz_trans, gz_cis))
    checkpoint(sprintf("Result files already exist: %s and/or %s. Skipping run.", gz_trans, gz_cis))
    message(sprintf("Result files already exist: %s and/or %s. Skipping run.", gz_trans, gz_cis))
    return(invisible(NULL))
  }

  tryCatch({
    checkpoint("Loading feature locations", verbose)
    feature_locations <- readRDS(feature_locations_path)

    checkpoint("Loading feature data", verbose)
    feature_data <- readRDS(feature_data_path)

    checkpoint("Loading SNP locations data", verbose)
    snp_locations <- readRDS(snpLocPath)

    checkpoint("Loaded feature locations, expression data, and SNP positions", verbose)

    cov_sliced <- load_sliced_file(covFilePath)
    snp_sliced <- load_sliced_file(snpFilePath)

    checkpoint("Sliced SNP and covariate data loaded", verbose)

    validate_sample_overlap(feature_data, snp_sliced, cov_sliced, group_name, verbose)

    checkpoint("Sample overlap validated", verbose)

    start_time <- Sys.time()

    checkpoint("Running Matrix_eQTL_main", verbose)

    result <- Matrix_eQTL_main(
      snps = snp_sliced,
      gene = feature_data,
      cvrt = cov_sliced,
      output_file_name = output_file_name_trans,
      pvOutputThreshold = pvOutputThreshold,
      useModel = useModelConst,
      errorCovariance = numeric(0),
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThresholdCis,
      snpspos = snp_locations,
      genepos = feature_locations,
      cisDist = cisDist,
      pvalue.hist = pvalueHist,
      min.pv.by.genesnp = minPvByGeneSnp,
      noFDRsaveMemory = noFDRsaveMemory,
      verbose = verbose
    )

    checkpoint("Matrix_eQTL_main completed", verbose)

    write.table(
      data.frame(result$all$ntests, result$trans$ntests, result$cis$ntests),
      file = ntests_path, row.names = FALSE, col.names = FALSE
    )
    checkpoint(paste("ntests written to:", ntests_path), verbose)

    result$trans$eqtls <- NULL
    result$cis$eqtls <- NULL

    saveRDS(list(
      result = result,
      feature_data_path = feature_data_path,
      feature_locations_path = feature_locations_path,
      snpFilePath = snpFilePath,
      covFilePath = covFilePath,
      snpLocPath = snpLocPath
    ), file = result_path)

    checkpoint(paste("RDS result saved to:", result_path), verbose)

    gzip(output_file_name_trans, destname = gz_trans, remove = TRUE)
    gzip(output_file_name_cis, destname = gz_cis, remove = TRUE)
    
    checkpoint("Compressed trans and cis output files", verbose)

    flog.info("MatrixEQTL complete in %s", round(Sys.time() - start_time, 2), name = "matrixeqtl")

    checkpoint("matrixEQTLwrapper finished", verbose)
    invisible(NULL)
  }, error = function(e) {
    flog.error("Error in matrixEQTLwrapper: %s", e$message, name = "matrixeqtl")
    stop(e)
  })
}

#' Run test harness for MatrixEQTL
#' @export
test_matrixEQTL_run <- function() {
  flog.threshold(DEBUG)
  flog.info("Running test harness...", name = "matrixeqtl")

  matrixEQTLwrapper(
    feature_locations_path = "test/features_locations_chunk_1.rds",
    feature_data_path = "test/expr_chunk_1.rds",
    snpFilePath = "test/snps_chunk_1.txt",
    covFilePath = "test/covariates.txt",
    snpLocPath = "test/snps_locations_chunk_1.txt",
    group_name = "group_test",
    resultsDir = "test/results",
    dry_run = TRUE,
    verbose = TRUE
  )

  flog.info("✅ Dry run test completed.", name = "matrixeqtl")
}

check_match <- function(sample_name, group_name_vector) {
  any(sapply(group_name_vector, function(x) grepl(sample_name, x)))
}

check_complete_match <- function(sample_name_vector, group_name_vector) {
  all(sapply(sample_name_vector, check_match, group_name_vector = group_name_vector))
}
#replace these two fucntion with the one below



# #' Validate and align sample columns between feature, SNP, and covariate data
# #'
# #' Ensures that the sample columns are correctly matched and ordered.
# #' Attempts to reorder SNPs and covariates to match expression samples if necessary.
# #'
# #' @param feature_data `SlicedData` object for expression data.
# #' @param snps `SlicedData` object for SNP data.
# #' @param covs `SlicedData` object for covariate data.
# #' @param group_name Group identifier (not used here but preserved for logging compatibility).
# #' @param verbose Logical. If `TRUE`, log reordering steps.
# #'
# #' @return None. Modifies `snps` and `covs` in-place.
# #' @export
# validate_sample_overlap <- function(feature_data, snps, covs, group_name, verbose = TRUE) {
#   if (covs$nRows() > feature_data$nCols()) stop("More covariates than samples.")
#   if (feature_data$nCols() > snps$nCols()) stop("More samples in features than in SNPs.")

#   feature_names <- feature_data$columnNames
#   snp_names <- snps$columnNames

#   if (!all(feature_names %in% snp_names)) {
#     if (verbose) flog.info("Reordering SNP and covariate data to match features", name = "matrixeqtl")
#     matched_indices <- match(feature_names, snp_names)
    
#     if (any(is.na(matched_indices))) {
#       stop("Some feature sample names not found in SNP data.")
#     }

#     snps$ColumnSubsample(matched_indices)
#     covs$ColumnSubsample(matched_indices)
#   }
# }



# #' Validate and align sample columns between feature, SNP, and covariate data
# #'
# #' Ensures that the sample columns are correctly matched and ordered.
# #'
# #' @param feature_data `SlicedData` object for expression data.
# #' @param snps `SlicedData` object for SNP data.
# #' @param covs `SlicedData` object for covariate data.
# #' @param group_name A group identifier for the sample group.
# #' @param verbose Logical. If `TRUE`, log reordering operations.
# #'
# #' @return None. Alters `snps` and `covs` by reference to match `feature_data`.
# #' @export
# validate_sample_overlap <- function(feature_data, snps, covs, group_name, verbose = TRUE) {
  
#   feature_samples <- feature_data$columnNames
#   snp_samples <- snps$columnNames
#   cov_samples <- covs$columnNames

#   # Check sizes
#   if (length(cov_samples) > length(feature_samples)) {
#     stop("More covariates than expression samples.")
#   }
#   if (length(feature_samples) > length(snp_samples)) {
#     stop("More expression samples than SNP samples.")
#   }

#   # If sample names don't match exactly, try reordering SNPs and covariates
#   if (!all(feature_samples %in% snp_samples)) {
#     if (verbose) flog.info("Reordering SNP and covariate data to match expression samples", name = "matrixeqtl")
    
#     match_idx <- match(feature_samples, snp_samples)
    
#     if (any(is.na(match_idx))) {
#       stop("Some expression sample names are missing in SNP data.")
#     }

#     snps$ColumnSubsample(match_idx)
    
#     if (!all(feature_samples %in% cov_samples)) {
#       stop("Some expression sample names are missing in covariates.")
#     }

#     covs$ColumnSubsample(match(feature_samples, cov_samples))
#   }
# }
