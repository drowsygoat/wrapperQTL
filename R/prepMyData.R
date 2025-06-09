library(SummarizedExperiment)
library(MatrixEQTL)
library(Matrix)
library(tibble)
library(dplyr)
library(parallel)

# Prepare feature locations function; helper function
prepareSlicedSubsets_helper <- function(group, group_data, feature_matrix, se_sample_data, matrixName, topNrows, common_samples, resultsDir, nChunks) {

    dir_name <- file.path(resultsDir, paste0("group_", group, "_results"))
    
    if (!dir.exists(dir_name)) {
      dir.create(dir_name, recursive = TRUE)
    }

    log_path <- file.path(dir_name, paste0(group, "_input_MEQTL_log.txt"))

    # Start time for execution timing

    # Execute  with error and warning handling
    sink(log_path, append = TRUE)
    on.exit(sink()) 

    tryCatch(withCallingHandlers({      

        group_indices <- which(group_data == group & se_sample_data %in% common_samples)
    
        feature_matrix <- feature_matrix[, group_indices, drop = FALSE]

        # Clone and update SNPs_sliced to include only common samples

        if (!is.null(topNrows)) {
            rowSums <- Matrix::rowSums(feature_matrix)
            top_indices <- order(rowSums, decreasing = TRUE)[1:min(topNrows, length(rowSums))]
            feature_matrix <- feature_matrix[top_indices, , drop = FALSE]
        }

        # Splitting the matrix into chunks if nChunks is not NULL
        feature_matrices <- list()
        if (!is.null(nChunks) && nChunks > 1) {
            rows_per_chunk <- ceiling(nrow(feature_matrix) / nChunks)
            for (i in 1:nChunks) {
                result_path <- file.path(dir_name, paste("group", group, "chunk", i, "input_MEQTL.rds", sep = "_"))
                start_row <- (i - 1) * rows_per_chunk + 1
                end_row <- min(i * rows_per_chunk, nrow(feature_matrix))
                feature_matrix_subset <- feature_matrix[start_row:end_row, , drop = FALSE]
                features_data <- SlicedData$new()
                features_data$CreateFromMatrix(as(feature_matrix_subset, "matrix"))
                features_data$ResliceCombined(sliceSize = 2000)
                saveRDS(features_data, file = result_path)
            }
        } else {
            result_path <- file.path(dir_name, paste("group", group, "chunk_1_input_MEQTL.rds", sep = "_"))
            features_data <- SlicedData$new()
            features_data$ResliceCombined(sliceSize = 2000)
            features_data$CreateFromMatrix(as(feature_matrix, "matrix"))
            saveRDS(features_data, file = result_path)

        }
        
# the original posityions on the list were retained in case they are needed, but they are

    }, warning = function(w) {
        # Handle warnings
        message(paste("Warning in group", group, ":", w$message))
        cat(paste(Sys.time(), "Warning:", w$message))
        invokeRestart("muffleWarning")
    }, error = function(e) {
        # Handle errors
        message(paste("Error in group", group, ":", e$message))
        cat(paste(Sys.time(), "Error:", e$message))
    }))

    # Save results and log execution time
    return(TRUE)
}

# Prepare feature locations function; helper function
prepareLocationsForMatrixEQTL <- function(data) {
  data <- tibble::rownames_to_column(data, var = "feature_id")
  result <- data %>%
    dplyr::select(feature_id, seqnames, start, end)
    return(result)
}


prepareSlicedSubsets <- function(pathToSE, resultsDir, matrixName, covFilePath = NULL, snpFilePath = NULL, check = FALSE, groupColName, sampleColName = "Sample", topNrows = NULL, groupSubset = NULL, nChunks = 1) {

  summarizedExp <- readRDS(pathToSE) 

  # Load necessary packages
  if (!requireNamespace("MatrixEQTL", quietly = TRUE)) {
    stop("Package 'MatrixEQTL' is needed but not installed.")
  }

  # Argument checks
  if (!inherits(summarizedExp, "SummarizedExperiment")) {
    stop("summarizedExp must be a SummarizedExperiment object.")
  }

  if (!matrixName %in% names(assays(summarizedExp))) {
    stop(paste("Matrix", matrixName, "not found in SummarizedExperiment object."))
  }

  if (!groupColName %in% names(colData(summarizedExp))) {
    stop(paste("groupColName", groupColName, "not found in colData of the SummarizedExperiment."))
  }

  if (!sampleColName %in% names(colData(summarizedExp))) {
    stop(paste("groupColName", groupColName, "not found in colData of the SummarizedExperiment."))
  }

  if (!is.null(groupSubset)) {
    cat("Subsetting SE.\n")
    summarizedExp <- summarizedExp[ , colData(summarizedExp)[[groupColName]] %in% groupSubset]
  }
    
  if (!dir.exists(resultsDir)) {
    dir.create(resultsDir, recursive = TRUE)
  }

  se_sample_ids <- unique(colData(summarizedExp)[[sampleColName]])
  se_sample_data <- colData(summarizedExp)[[sampleColName]] # sample column from colData

  feature_matrix <- as(assays(summarizedExp)[[matrixName]], "sparseMatrix")

  if (!is.null(nChunks) && nChunks > 1) {
      rows_per_chunk <- ceiling(nrow(feature_matrix) / nChunks)
      for (i in 1:nChunks) {

          feature_locations_path <- file.path(resultsDir, paste("chunk", i, "loc_input_MEQTL.rds", sep = "_"))
          start_row <- (i - 1) * rows_per_chunk + 1
          end_row <- min(i * rows_per_chunk, nrow(feature_matrix))

          feature_locations_subsetted <- prepareLocationsForMatrixEQTL(as.data.frame(rowData(summarizedExp))[start_row:end_row, , drop = FALSE])
          saveRDS(feature_locations_subsetted, feature_locations_path)

      }

  } else {
    
      feature_locations <- prepareLocationsForMatrixEQTL(as.data.frame(rowData(summarizedExp)))
      feature_locations_path <- file.path(resultsDir, "chunk_1_loc_input_MEQTL.rds")
      saveRDS(feature_locations, feature_locations_path)
      
  }
        
  group_data <- colData(summarizedExp)[[groupColName]]
  groups_to_process <- unique(group_data)

  rm(summarizedExp)

  cat("Sample IDs in SE:", se_sample_ids, "\n")

    # Covariate Data
    COV_sliced <- SlicedData$new()
    COV_sliced$fileDelimiter <- "\t"
    COV_sliced$fileOmitCharacters <- "NA"
    COV_sliced$fileSkipRows <- 1
    COV_sliced$fileSkipColumns <- 1
    COV_sliced$fileSliceSize <- 2000
    COV_sliced$LoadFile(covFilePath)
    
    cov_sample_ids <- COV_sliced$columnNames

  if (check) {
    if (is.null(covFilePath) | is.null(snpFilePath)) {
      stop("'check = TRUE' but 'covFilePath' and/or 'snpFilePath' were not provided")
    }

    SNPs_sliced <- SlicedData$new()
    SNPs_sliced$fileDelimiter <- "\t"
    SNPs_sliced$fileOmitCharacters <- "NA"
    SNPs_sliced$fileSkipRows <- 1
    SNPs_sliced$fileSkipColumns <- 1
    SNPs_sliced$fileSliceSize <- 2000
    SNPs_sliced$LoadFile(snpFilePath)

    snp_sample_ids <- SNPs_sliced$columnNames

    if (!all(snp_sample_ids == cov_sample_ids)) {
      stop("SNPs and covariances sample IDs are not matching.")
    }
  }

  cat("Found ", length(cov_sample_ids), " sample IDs in the provided covariance matrix.\n")

  # common samples between groups
  common_samples <- intersect(se_sample_ids, cov_sample_ids)

  if (length(common_samples) / length(se_sample_ids) < 0.5) {
    warning("Less than 50% of the sample IDs in the SE object match those in the SNP file.")
  } 

  cat("Matched sample IDs:", common_samples, "(", length(common_samples), " out of ", length(se_sample_ids), ")\n")

  # num_cores <- max(1, detectCores() - 1)

  # Run for each group
  lapply(groups_to_process, prepareSlicedSubsets_helper, feature_matrix =     feature_matrix, group_data = group_data, se_sample_data = se_sample_data,topNrows = topNrows, resultsDir = resultsDir, common_samples = common_samples, matrixName = matrixName, nChunks = nChunks)

  return(TRUE)
}
 
#  mc.preschedule = FALSE, mc.cores = num_cores
# add argument seurat_obj that could be used instead of SE argument. Then matrix name will be the surta object slot (e.g., RNA, SCT, etc.)

# seurat_cl_v3_group_SE_scale_1000 <- readRDS("/proj/sllstore2017078/private/lech_rackham/scAnalysis_rackham/milestones/peaks_for_QTL/seurat_cl_v3_group_SE_scale_1000.rds") 

# res_pilot3 <- prepareSlicedSubsets(pathToSE = "/proj/sllstore2017078/private/lech_rackham/scAnalysis_rackham/milestones/peaks_for_QTL/seurat_cl_v3_group_SE_scale_1000.rds",

# resultsDir = "seurat_cl_v3_group_SE_scale_1000_MEQTL_one_chunk", 

# matrixName = "PeakMatrix", 

# snpFilePath = "/proj/sllstore2017078/private/lech_rackham/scAnalysis_rackham/QTL_analysis/atac_QTL/data/SNPs_VK_fixed_colnames.txt", 

# covFilePath = "/proj/sllstore2017078/private/lech_rackham/scAnalysis_rackham/QTL_analysis/atac_QTL/data/COV_VK_fixed_colnames.txt", 

# groupColName = "Cluster",

# nChunks = 1)
