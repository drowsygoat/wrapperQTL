#' Prepare MatrixEQTL Input Files
#'
#' This function prepares MatrixEQTL input files (features, feature locations, SNPs, and covariates)
#' from a SummarizedExperiment or data.frame object. SNP and covariate data can be given as multiple files.
#'
#' @param features Either a SummarizedExperiment object or a numeric data.frame/matrix with features x samples.
#' @param sampleMetadata Optional: A data.frame with sample annotations (required if `features` is not an SE).
#' @param groupCol String; name of the column in sample metadata for groupings (e.g., clusters or conditions).
#' @param sampleCol String; name of the sample ID column.
#' @param snpPaths Character vector; paths to one or more SNP data files (MatrixEQTL format).
#' @param covPaths Character vector; paths to one or more covariate files (MatrixEQTL format).
#' @param resultsDir String; directory to store MatrixEQTL input files.
#' @param matrixName Optional; assay name if using SummarizedExperiment (default = first assay).
#' @param nChunks Integer; number of chunks to split the feature matrix into.
#' @param topNrows Optional integer; select top N most expressed features per group.
#' @param topNSNPs Optional integer; select top N most variable SNPs across samples.
#' @param groupSubset Optional; vector of group labels to process.
#' @param sliceSize Integer; slice size for SlicedData objects.
#' @param gtfFile Optional; GTF file used to retrieve feature locations if not present in the SE.
#' @param useGeneName Logical; whether to use gene_name instead of gene_id from GTF file.
#'
#' @return Logical TRUE (invisible) if execution succeeds.
prepareMatrixEQTLInputs <- function(
  verbose = TRUE,
  features,
  sampleMetadata = NULL,
  groupCol,
  sampleCol = "Sample",
  snpPaths,
  covPaths,
  resultsDir,
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

  if (verbose) message("Extracting and validating input features...")
  # Input: extract matrices
  is_SE <- inherits(features, "SummarizedExperiment")
  if (is_SE) {
    if (is.null(matrixName)) matrixName <- names(assays(features))[1]
    stopifnot(matrixName %in% names(assays(features)))
    feature_matrix <- as(assays(features)[[matrixName]], "sparseMatrix")
    sample_info <- as.data.frame(colData(features))

    # Attempt to retrieve feature locations from rowData or GTF
    if (ncol(rowData(features)) >= 3 && all(c("seqnames", "start", "end") %in% colnames(rowData(features)))) {
      feature_locs <- as.data.frame(rowData(features)) %>%
        tibble::rownames_to_column("feature_id") %>%
        dplyr::select(feature_id, seqnames, start, end)
    } else if (!is.null(gtfFile)) {
      message("Retrieving feature locations from GTF file...")
      feature_locs <- rtracklayer::import(gtfFile) %>%
        as.data.frame() %>%
        dplyr::filter(type == "gene") %>%
        dplyr::select(gene_id, gene_name, seqnames, start, end) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
          chr = as.character(seqnames),
          geneid = case_when(useGeneName ~ gene_name,
                             is.na(gene_name) ~ gene_id,
                             TRUE ~ gene_name
          ),
          s1 = start,
          s2 = end
        ) %>% 
        dplyr::select(feature_id = geneid, seqnames = chr, start = s1, end = s2) %>% 
        dplyr::filter(feature_id %in% rownames(features)
        )
    } else {
      warning("Feature locations not found in rowData and no GTF file provided. Skipping feature location export.")
      feature_locs <- NULL
    }
  } else {
    stopifnot(is.matrix(features) || is.data.frame(features))
    feature_matrix <- as(features, "sparseMatrix")
    stopifnot(!is.null(sampleMetadata))
    sample_info <- sampleMetadata
    feature_locs <- NULL
  }

  stopifnot(groupCol %in% colnames(sample_info), sampleCol %in% colnames(sample_info))
  if (!is.null(groupSubset)) {
    keep <- sample_info[[groupCol]] %in% groupSubset
    feature_matrix <- feature_matrix[, keep, drop = FALSE]
    sample_info <- sample_info[keep, , drop = FALSE]
  }

  sample_ids <- sample_info[[sampleCol]]
  group_data <- sample_info[[groupCol]]
  groups <- unique(group_data)

  # Merge SNP and Covariate sample names
  getSampleNamesFromFile <- function(path) {
    if (file.info(path)$isdir) {
      files <- list.files(path, pattern = "(?i)(cov|snp)", full.names = TRUE)
      if (length(files) == 0) stop("No SNP or covariate files found in ", path)
      path <- files
    }

    file_tables <- lapply(path, function(p) {
      if (!file.exists(p)) stop("File does not exist: ", p)
      if (verbose) {
        message("
📄 Reading file: ", p)
        header <- readLines(p, n = 2)
        message("Header: ", header[1])
        message("First row: ", header[2])
      }
      df <- read.table(p, header = TRUE, check.names = FALSE)
      if (nrow(df) == 0 || ncol(df) <= 1) stop("File has insufficient content: ", p)
      df
    })

    merged_df <- Reduce(function(x, y) {
      x_names <- colnames(x)[-1]
      y_names <- colnames(y)[-1]
      if (!identical(sort(x_names), sort(y_names))) {
        stop("Column names (sample IDs) do not match between covariate files.")
      }
      colnames(x)[1] <- "covariate"
      colnames(y)[1] <- "covariate"
      rbind(x, y)
    }, file_tables)

    colnames(merged_df)[-1]  # return sample names (skip ID column)
  }
  snp_samples <- Reduce(intersect, lapply(snpPaths, getSampleNamesFromFile))
  cov_samples <- Reduce(intersect, lapply(covPaths, getSampleNamesFromFile))
  common_samples <- intersect(intersect(sample_ids, snp_samples), cov_samples)

  if (verbose) message("Matched ", length(common_samples), " samples across SNP, covariates, and features.")

  # Save feature locations to RDS and TXT
  if (!is.null(feature_locs)) {
    if (verbose) message("Writing feature_location.txt to ", file.path(resultsDir, "feature_location.txt"))
    write.table(
      feature_locs,
      file = file.path(resultsDir, "feature_location.txt"),
      sep = "	",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    );
    rows_per_chunk <- ceiling(nrow(feature_locs) / nChunks)
    for (i in seq_len(nChunks)) {
      start <- (i - 1) * rows_per_chunk + 1
      end <- min(i * rows_per_chunk, nrow(feature_locs))
      chunk <- feature_locs[start:end, , drop = FALSE]
      saveRDS(chunk, file.path(resultsDir, paste0("chunk_", i, "_loc_input_MEQTL.rds")))
    }
  }

  # Save expression matrix chunks per group
  for (group in groups) {
    if (verbose) message("Processing group: ", group)
    group_idx <- which(group_data == group & sample_ids %in% common_samples)
    if (length(group_idx) == 0) next

    mat <- feature_matrix[, group_idx, drop = FALSE]
    if (!is.null(topNrows)) {
      top_idx <- order(Matrix::rowSums(mat), decreasing = TRUE)[1:min(topNrows, nrow(mat))]
      mat <- mat[top_idx, , drop = FALSE]
    }

    group_dir <- file.path(resultsDir, paste0("group_", group, "_results"))
    if (!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)
    if (verbose) message("Creating directory: ", group_dir)
    dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)

    rows_per_chunk <- ceiling(nrow(mat) / nChunks)
    for (i in seq_len(nChunks)) {
      start <- (i - 1) * rows_per_chunk + 1
      end <- min(i * rows_per_chunk, nrow(mat))
      chunk <- mat[start:end, , drop = FALSE]
      if (!is.null(topNSNPs)) {
        snp_var <- Matrix::colSums(chunk != 0)
        top_snps <- order(snp_var, decreasing = TRUE)[1:min(topNSNPs, ncol(chunk))]
        chunk <- chunk[, top_snps, drop = FALSE]
      }
      sliced <- SlicedData$new()
      sliced$CreateFromMatrix(as(chunk, "matrix"))
      sliced$ResliceCombined(sliceSize)
      output_path <- file.path(group_dir, paste0("group_", group, "_chunk_", i, "_input_MEQTL.rds"))
      if (verbose) message("Saving feature chunk to ", output_path)
      saveRDS(sliced, output_path)
    }
  }
  # return(feature_locs)
  invisible(TRUE)
}
