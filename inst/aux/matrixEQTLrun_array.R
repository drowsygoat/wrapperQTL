#!/usr/bin/env Rscript

# Load necessary libraries
devtools::load_all("/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/wrapperQTL")

library(argparse)
library(MatrixEQTL)
library(parallel)
library(R.utils)
library(futile.logger)

# Define command-line arguments
parser <- ArgumentParser(description = 'Process QTL analysis arguments')

parser$add_argument("group_name", help = "Name of the group")
parser$add_argument("outputDir", help = "Output directory")
parser$add_argument("snpFilePath", help = "Path to the directory containing SNP chunks")
parser$add_argument("snpLocPath", help = "Path to the directory containing SNP location files")
parser$add_argument("feature_data_path", help = "Path to the directory containing expression data chunks")
parser$add_argument("feature_locations_path", help = "Path to the directory containing expression location files")
parser$add_argument("covFilePath", help = "Path to the covariates file")
parser$add_argument("--threads", type = "integer", default = 10,
                    help = "Number of threads to use (default: 10)")
parser$add_argument("--dry_run", action = "store_true",
                    help = "Dry run mode — don't run analysis, just simulate")
parser$add_argument("--jobArrayID", type = "integer", default = NULL,
                    help = "Chunk ID (e.g. from SLURM_ARRAY_TASK_ID)")
parser$add_argument("--verbose", action = "store_true",
                    help = "Enable verbose logging")

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
  pvOutputThreshold = 1e-6,
  pvOutputThresholdCis = 1e-5,
  useModel = "linear",
  minPvByGeneSnp = FALSE,
  noFDRsaveMemory = TRUE,
  SNPsInChunks = FALSE,
  prefix = "defaultprefix",
  pvalueHist = "qqplot",
  threads = args$threads,
  dry_run = args$dry_run,
  chunk_id = args$jobArrayID

)
