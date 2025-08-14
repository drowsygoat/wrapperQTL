#' Prepare Matrix eQTL Input Files
#'
#' Prepares and filters expression, SNP, and covariate matrices for Matrix eQTL analysis.
#' Supports SummarizedExperiment or raw matrices as input, optional GTF-based annotation,
#' and automatic chunking for large datasets. Saves intermediate results as RDS files.
#'
#' @param verbose Logical; whether to print progress messages.
#' @param features SummarizedExperiment object or expression matrix (genes x samples).
#' @param sampleMetadata Optional data.frame with sample annotations (required if `features` is not a SE).
#' @param groupCol Column in metadata indicating group assignment.
#' @param sampleCol Column in metadata indicating sample identifiers (default: "sample_id").
#' @param snpPaths Character vector of paths to SNP matrices or directories.
#' @param covPaths Character vector of paths to covariate matrices or directories.
#' @param resultsDir Directory to save all output files.
#' @param minFeatureFrac Minimum fraction of samples a gene must be detected in (default: 0.8).
#' @param minFeatureMean Minimum mean expression for features to be retained (default: 0).
#' @param matrixName Name of assay to extract from SummarizedExperiment (optional).
#' @param nChunks Number of chunks to split SNP and expression matrices into (default: 5).
#' @param topNrows Optional limit on top expressed features to retain.
#' @param topNSNPs Optional limit on top SNPs by variance to retain.
#' @param groupSubset Optional vector of group names to include.
#' @param sliceSize Slice size for SlicedData export (default: 2000).
#' @param gtfFile Optional GTF file to retrieve genomic coordinates for features.
#' @param useGeneName Logical; whether to prioritize gene name over gene ID (default: FALSE).
#'
#' @return Invisible TRUE; writes RDS and text files to `resultsDir`.
#' @export
prepareMatrixEQTLInputs <- function(
  verbose = TRUE,
  features,
  sampleMetadata = NULL,
  groupCol = NULL,
  sampleCol = "sample_id",
  snpPaths,
  covPaths,
  resultsDir,
  minFeatureFrac = 0.8,
  minFeatureMean = 0,
  matrixName = NULL,
  nChunks = 1,
  topNrows = NULL,
  topNSNPs = NULL,
  groupSubset = NULL,
  sliceSize = 2000,
  gtfFile = NULL,
  useGeneName = FALSE
) {
  if (!dir.exists(resultsDir)) {
    if (verbose) message("Creating results directory: ", resultsDir)
    dir.create(resultsDir, recursive = TRUE)
  }

  filter_features_by_expression <- function(mat, min_frac = 0.1, min_mean = 0) {
    stopifnot(is(mat, "sparseMatrix"))
    n_samples <- ncol(mat)
    detect_frac <- Matrix::rowSums(mat > 0) / n_samples
    mean_expr <- Matrix::rowMeans(mat)
    keep <- detect_frac >= min_frac & mean_expr >= min_mean
    mat[keep, , drop = FALSE]
  }

  extract_snp_locations <- function(snp_input) {
    if (is.character(snp_input) && file.exists(snp_input)) {
      snp_df <- read.table(snp_input, header = TRUE, check.names = FALSE)
    } else if (is.data.frame(snp_input)) {
      snp_df <- snp_input
    } else {
      stop("snp_input must be a valid file path or a data frame.")
    }

    snp_ids <- snp_df[[1]]
    parsed_df <- do.call(rbind, lapply(snp_ids, function(x) {
      parts <- strsplit(x, "[:_]")[[1]]
      if (length(parts) >= 2) {
        chr <- paste0("chr", parts[1])
        pos <- as.integer(parts[2])
      } else {
        chr <- NA
        pos <- NA
      }
      return(c(chr, pos))
    }))

    locs <- data.frame(
      snpid = snp_ids,
      chr = parsed_df[, 1],
      pos = as.integer(parsed_df[, 2]),
      stringsAsFactors = FALSE
    )

    return(list(data = snp_df, locs = locs))
  }

  read_and_merge_files <- function(paths, pattern = NULL) {
    resolved_paths <- lapply(paths, function(p) {
      if (!is.null(pattern) && file.info(p)$isdir) {
        files <- list.files(p, pattern = pattern, full.names = TRUE)
        if (length(files) == 0) stop("No files matching pattern '", pattern, "' found in directory: ", p)
        if (length(files) > 1) stop("Multiple files matching pattern '", pattern, "' found in directory: ", p)
        return(files[1])
      } else {
        if (!file.exists(p)) stop("File does not exist: ", p)
        return(p)
      }
    })

    resolved_paths <- unlist(resolved_paths)

    file_tables <- lapply(resolved_paths, function(p) {
      df <- read.table(p, header = TRUE, check.names = FALSE)
      if (nrow(df) == 0 || ncol(df) <= 1) stop("File has insufficient content: ", p)
      df
    })

    col_sets <- lapply(file_tables, function(df) colnames(df)[-1])
    common_cols <- Reduce(intersect, col_sets)
    if (length(common_cols) == 0) stop("No common sample columns found across input files.")

    file_tables <- lapply(file_tables, function(df) {
      df <- df[, c(colnames(df)[1], common_cols), drop = FALSE]
      colnames(df)[1] <- "row_id"  # ensure the first column is uniformly named
      return(df)
    })

    merged <- do.call(rbind, file_tables)
    return(merged)
  }

  is_SE <- inherits(features, "SummarizedExperiment")
  if (is_SE) {
    if (is.null(matrixName)) matrixName <- names(assays(features))[1]
    stopifnot(matrixName %in% names(assays(features)))
    feature_matrix <- as(assays(features)[[matrixName]], "sparseMatrix")
    if (verbose) message("Filtering features based on expression thresholds...")
    feature_matrix <- filter_features_by_expression(feature_matrix, min_frac = minFeatureFrac, min_mean = minFeatureMean)
    sample_info <- as.data.frame(colData(features))

    if (ncol(rowData(features)) >= 3 && all(c("seqnames", "start", "end") %in% colnames(rowData(features)))) {
      feature_locs <- as.data.frame(rowData(features)) |>
        tibble::rownames_to_column("feature_id") |>
        dplyr::select(feature_id, seqnames, start, end)
    } else if (!is.null(gtfFile)) {
      if (verbose) message("Retrieving feature locations from GTF file...")
      feature_locs <- rtracklayer::import(gtfFile) |>
        as.data.frame() |>
        dplyr::filter(type == "gene") |>
        dplyr::select(gene_id, gene_name, seqnames, start, end) |>
        dplyr::distinct() |>
        dplyr::mutate(
          chr = as.character(seqnames),
          geneid = dplyr::case_when(
            useGeneName ~ gene_name,
            is.na(gene_name) ~ gene_id,
            TRUE ~ gene_name
          )
        ) |>
        dplyr::select(feature_id = geneid, seqnames = chr, start, end) |>
        dplyr::filter(feature_id %in% rownames(features)) |>
        dplyr::mutate(seqnames = ifelse(tolower(seqnames) == "mt", "chrM", paste0("chr", seqnames)))
    } else {
      warning("Feature locations not found in rowData and no GTF file provided.")
      feature_locs <- NULL
    }
  } else {
    stopifnot(is.matrix(features) || is.data.frame(features))
    stopifnot(!is.null(sampleMetadata))
    feature_matrix <- as(features, "sparseMatrix")
    sample_info <- sampleMetadata
    feature_locs <- NULL
  }

  # Ensure sampleCol is valid, or fall back to rownames
  if (is.null(sampleCol) || !sampleCol %in% colnames(sample_info)) {
    if (verbose) message("sampleCol is NULL or not found in colData. Using rownames(colData) as sample IDs.")
    sampleCol <- "__sample__"
    sample_info[[sampleCol]] <- rownames(sample_info)
  }

  # Ensure groupCol is valid, or treat all samples as one group
  if (is.null(groupCol) || !groupCol %in% colnames(sample_info)) {
    if (verbose) message("groupCol is NULL or not found in colData. Treating all samples as one group.")
    sample_info$`__group__` <- "all_samples"
    groupCol <- "__group__"
  }

  if (!is.null(groupSubset)) {
    keep <- sample_info[[groupCol]] %in% groupSubset
    feature_matrix <- feature_matrix[, keep, drop = FALSE]
    sample_info <- sample_info[keep, , drop = FALSE]
  }

  sample_ids <- sample_info[[sampleCol]]
  group_data <- sample_info[[groupCol]]

  groups <- unique(group_data)
  snp_df_merged <- read_and_merge_files(snpPaths)
  snp_samples <- colnames(snp_df_merged)[-1]
  common_samples <- intersect(sample_ids, snp_samples)
print(sample_ids)
print(snp_samples)


  if (verbose) message("Matched ", length(common_samples), " samples across SNP and features.")

  parsed_snp <- extract_snp_locations(snp_df_merged)
  snp_matrix <- parsed_snp$data
  snp_locs <- parsed_snp$locs

  if (!is.null(topNSNPs)) {
    snp_mat <- as.matrix(snp_matrix[, -1])
    snp_var <- apply(snp_mat, 1, function(x) var(as.numeric(x), na.rm = TRUE))
    top_idx <- order(snp_var, decreasing = TRUE)[1:min(topNSNPs, nrow(snp_matrix))]
    snp_matrix <- snp_matrix[top_idx, , drop = FALSE]
    snp_locs <- snp_locs[top_idx, , drop = FALSE]
  }

  if (nChunks <= 0) {
    saveRDS(snp_matrix, file.path(resultsDir, "merged_SNPs.rds"))
    saveRDS(snp_locs, file.path(resultsDir, "merged_SNP_locations.rds"))
  } else {
    snp_chunks_dir <- file.path(resultsDir, "snp_chunks_rds")
    dir.create(snp_chunks_dir, showWarnings = FALSE, recursive = TRUE)
    rows_per_chunk <- ceiling(nrow(snp_matrix) / nChunks)

    for (i in seq_len(nChunks)) {
      start <- (i - 1) * rows_per_chunk + 1
      end <- min(i * rows_per_chunk, nrow(snp_matrix))
      chunk <- snp_matrix[start:end, , drop = FALSE]
      locs_chunk <- snp_locs[start:end, , drop = FALSE]
      saveRDS(chunk, file.path(snp_chunks_dir, paste0("chunk_", i, "_SNPs.rds")))
      saveRDS(locs_chunk, file.path(snp_chunks_dir, paste0("chunk_", i, "_SNP_locations.rds")))
    }
  }

  for (group in groups) {
    if (verbose) message("📦 Processing group: ", group)
    group_idx <- which(group_data == group & sample_ids %in% common_samples)
    if (length(group_idx) == 0) next

    mat <- feature_matrix[, group_idx, drop = FALSE]
    if (!is.null(topNrows)) {
      top_idx <- order(Matrix::rowSums(mat), decreasing = TRUE)[1:min(topNrows, nrow(mat))]
      mat <- mat[top_idx, , drop = FALSE]
    }

    group_dir <- file.path(resultsDir, paste0("group_", group, "_results"))
    dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)

    if (!is.null(feature_locs)) {
      loc_subset <- feature_locs[feature_locs$feature_id %in% rownames(mat), , drop = FALSE]
      saveRDS(loc_subset, file.path(group_dir, paste0("group_", group, "_feature_locations.rds"))) # saving all locations, is that redundant?

      rows_per_chunk <- ceiling(nrow(loc_subset) / nChunks)
      for (i in seq_len(nChunks)) {
        chunk <- loc_subset[((i - 1) * rows_per_chunk + 1):min(i * rows_per_chunk, nrow(loc_subset)), , drop = FALSE]
        saveRDS(chunk, file.path(group_dir, paste0("group_", group, "_chunk_", i, "_loc_input_MEQTL.rds")))
      }
    }

    rows_per_chunk <- ceiling(nrow(mat) / nChunks)
    for (i in seq_len(nChunks)) {
      chunk <- mat[((i - 1) * rows_per_chunk + 1):min(i * rows_per_chunk, nrow(mat)), , drop = FALSE]
      sliced <- SlicedData$new()
      sliced$CreateFromMatrix(as(chunk, "matrix"))
      sliced$ResliceCombined(sliceSize)
      saveRDS(sliced, file.path(group_dir, paste0("group_", group, "_chunk_", i, "_input_MEQTL.rds")))
    }

    cov_df_merged <- read_and_merge_files(covPaths, pattern = paste0(group, "_"))
    cov_samples <- colnames(cov_df_merged)[-1]
    common_samples <- intersect(intersect(sample_ids, snp_samples), cov_samples)
    if (verbose) message("Matched ", length(common_samples), " samples across SNP, covariates, and features.")
    cov_export_path <- file.path(group_dir, "merged_covariates.txt")
    write.table(cov_df_merged, file = cov_export_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(TRUE)
}
