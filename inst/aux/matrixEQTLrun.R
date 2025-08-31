#!/usr/bin/env Rscript

library(argparse)
library(MatrixEQTL)
library(parallel)
library(R.utils)
library(futile.logger)
library(wrapperQTL)

# Define command-line arguments
parser <- ArgumentParser(description = 'Process QTL analysis arguments')

parser$add_argument("group_name", help = "Name of the group")
parser$add_argument("outputDir", help = "Output directory")
parser$add_argument("snpFilePath", help = "Path to the directory containing SNP chunks")
parser$add_argument("snpLocPath", help = "Path to the directory containing SNP location files")
parser$add_argument("feature_data_path", help = "Path to the directory containing expression data chunks")
parser$add_argument("feature_locations_path", help = "Path to the directory containing expression location files")
parser$add_argument("covFilePath", help = "Path to the covariates file")
parser$add_argument("--threads", type = "integer", default = 10, help = "Number of threads to use (default: 10)")

parser$add_argument("--pvOutputThresholdTrans", type = "double", default = 1e-5)

parser$add_argument("--pvOutputThresholdCis", type = "double", default = 1e-4)

parser$add_argument("--dry_run", action = "store_true", help = "Dry run mode — don't run analysis, just simulate")
# parser$add_argument("--chunk_id", action = "integer", help = "Chunk id array")  # Not used currently
parser$add_argument("--verbose", action = "store_true", help = "Enable verbose logging")

parser$add_argument(
  "--colNames_convention",
  type = "character",
  nargs = 2,
  help = "Two-element character vector specifying column name conventions, e.g., gene and sample naming pattern"
)
# Parse command-line arguments
args <- parser$parse_args()

# Optional debugging output
print("----- ARGUMENTS -----")
print(args)
print("---------------------")

# Run the MatrixEQTL multi-chunk wrapper
matrixEQTLwrapperMC(
  feature_locations_path = args$feature_locations_path,
  feature_data_path = args$feature_data_path,
  snpFilePath = args$snpFilePath,
  snpLocPath = args$snpLocPath,
  covFilePath = args$covFilePath,
  group_name = args$group_name,
  resultsDir = args$outputDir,
  cisDist = 1e6,
  pvOutputThreshold = args$pvOutputThresholdTrans,
  pvOutputThresholdCis = args$pvOutputThresholdCis,
  useModel = "linear",
  minPvByGeneSnp = FALSE,
  noFDRsaveMemory = TRUE,
  SNPsInChunks = FALSE,
  prefix = "defaultprefix",
  pvalueHist = "qqplot",
  threads = args$threads,
  dry_run = args$dry_run,
  colNames_convention = args$colNames_convention,
  verbose = args$verbose
)
