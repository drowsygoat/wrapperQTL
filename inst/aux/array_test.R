#!/usr/bin/env Rscript

# Load necessary libraries
message("--------START--")

devtools::load_all("/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/wrapperQTL")

library(argparse)
library(parallel)
library(R.utils)
library(futile.logger)

# Define command-line arguments
parser <- ArgumentParser(description = 'Process QTL analysis arguments')

parser$add_argument("--threads", type = "integer", default = 10, help = "Number of threads to use (default: 10)")

parser$add_argument("--jobArrayID", type = "integer", nargs = "?", default = 1, help = "SlurmTaskID")

# Parse command-line arguments
args <- parser$parse_args()

# Optional debugging output
print("----- ARGUMENTS -----")
print(args)
print("---------------------")

# Run the MatrixEQTL multi-chunk wrapper
message("Hello World ", Sys.getenv("SLURM_ARRAY_TASK_ID"))
message("Hello World ", args$jobArrayID)
message("----------")

