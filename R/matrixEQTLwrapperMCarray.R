matrixEQTLwrapperMCarray <- function(
  feature_locations_path,
  feature_data_path,
  snpFilePath,
  snpLocPath,
  covFilePath,
  group_name,
  chunk_id,
  resultsDir = getwd(),
  cisDist = 1e6,
  pvOutputThreshold = 1e-5,
  pvOutputThresholdCis = 1e-4,
  useModel = "linear",
  minPvByGeneSnp = TRUE,
  noFDRsaveMemory = FALSE,
  pvalueHist = NULL,
  SNPsInChunks = NULL,
  prefix = NULL,
  threads = 1,
  dry_run = FALSE,
  verbose = TRUE
) {
  message("[INFO] Checking input directories...")

  prepare_iteration_df <- function(feature_locations_path, feature_data_path, snpFilePath, snpLocPath, group_name, threads = NULL) {

  check_directories(feature_locations_path, feature_data_path, snpFilePath, snpLocPath)

  DATA <- list.files(file.path(feature_data_path, group_name), pattern = "chunk_[0-9]+_input", full.names = TRUE)

# print(DATA)

    SNP  <- list.files(snpFilePath, pattern = "chunk_[0-9]+_SNPs", full.names = TRUE)

# print(SNP)

    iteration_df <- expand.grid(
      DATA = DATA,
      SNP = SNP,
      stringsAsFactors = FALSE
    )
    
print(iteration_df)

    chunk_snp     <- grep_o(iteration_df$SNP, "chunk_[0-9]+")
    chunk_feature <- grep_o(iteration_df$DATA, "chunk_[0-9]+")

print(chunk_snp)
print(chunk_feature)

    iteration_df$SNP_LOC <- sapply(chunk_snp, function(x) {
      list.files(snpLocPath, pattern = paste0(x, "_SNP_loc"), full.names = TRUE)[[1]]
    }, USE.NAMES = FALSE)

    iteration_df$DATA_LOC <- sapply(chunk_feature, function(x) {
      list.files(file.path(feature_locations_path, group_name), pattern = paste0(x, "_loc"), full.names = TRUE)[[1]]
    }, USE.NAMES = FALSE)

    return(iteration_df)
  }

  iteration_df <- prepare_iteration_df(feature_locations_path, feature_data_path, snpFilePath, snpLocPath, group_name, threads)

  stopifnot(nrow(iteration_df) == length(iteration_df$SNP_LOC),
            nrow(iteration_df) == length(iteration_df$DATA_LOC))

  if (dry_run) {
    message("[DRY RUN] The following MatrixEQTL jobs would be executed:")
    for (i in seq_len(nrow(iteration_df))) {
      message(sprintf("  [%d/%d] Expression: %s | SNP: %s | SNP loc: %s | Expr loc: %s",
                      i, nrow(iteration_df),
                      basename(iteration_df$DATA[[chunk_id]]),
                      basename(iteration_df$SNP[[chunk_id]]),
                      basename(iteration_df$SNP_LOC[[chunk_id]]),
                      basename(iteration_df$DATA_LOC[[chunk_id]])))
    }
    return(invisible(NULL))
  }
  message(sprintf("[INFO] Starting Matrix eQTL chunked analysis for %s", chunk_id))

  run_job <- function(i) {
    matrixEQTLwrapper(
      feature_locations_path = iteration_df$DATA_LOC[[i]],
      feature_data_path      = iteration_df$DATA[[i]],
      snpFilePath            = iteration_df$SNP[[i]],
      covFilePath            = covFilePath,
      snpLocPath             = iteration_df$SNP_LOC[[i]],
      group_name             = group_name,
      resultsDir             = resultsDir,
      cisDist                = cisDist,
      pvOutputThreshold      = pvOutputThreshold,
      pvOutputThresholdCis   = pvOutputThresholdCis,
      useModel               = useModel,
      minPvByGeneSnp         = minPvByGeneSnp,
      noFDRsaveMemory        = noFDRsaveMemory,
      pvalueHist             = pvalueHist,
      # SNPsInChunks           = SNPsInChunks,
      prefix                 = prefix,
      dry_run                = dry_run,
      verbose                = verbose
    )
  }

  run_job(chunk_id)

  message(sprintf("[INFO] Array element %s processed.", chunk_id))

}

check_directories <- function(...) {
  paths <- list(...)
  if (length(paths) == 0) stop("No paths provided.")
  for (path in paths) {
    if (!dir.exists(path)) stop(sprintf("Directory does not exist: %s", path))
  }
}
