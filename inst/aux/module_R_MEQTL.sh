#!/bin/bash

# MODULE PARAMETERS
RUN_COMMAND="run_shell_command.sh"
JOB_NAME="5atac_${1}"
PARTITION="main"
NODES=1
TIME="23:59:00"
TASKS=1
# MEMORY="512GB"
# JOB_ARRAY="1-4"
CPUS=8
DRY="slurm"

threads=1

# feature locations is the same as features data, the directory changes is resolved in the R script

# snp_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/toy_data/results_toy/snp_chunks_rds"

# features_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/toy_data/results_toy"

# snp_loc_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/toy_data/results_toy/snp_chunks_rds"

# features_loc_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/toy_data/results_toy"

# covFilePath="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/toy_data/results_toy/merged_covariates.txt"

# snp_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_08_0_4chk/snp_chunks_rds"

# features_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_08_0_4chk"

# snp_loc_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_08_0_4chk/snp_chunks_rds"

# features_loc_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_08_0_4chk"

snp_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_05_002_20chk_ATAC_all_clusters/snp_chunks_rds"

features_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_05_002_20chk_ATAC_all_clusters"

snp_loc_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_05_002_20chk_ATAC_all_clusters/snp_chunks_rds"

features_loc_path="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_05_002_20chk_ATAC_all_clusters"

# Check if directories exist
if [ ! -d "$snp_path" ]; then
    echo "Error: Directory '$snp_path' does not exist."
    exit 1
fi

if [ ! -d "$features_path" ]; then
    echo "Error: Directory '$features_path' does not exist."
    exit 1
fi

if [ ! -d "$snp_loc_path" ]; then
    echo "Error: Directory '$snp_loc_path' does not exist."
    exit 1
fi

if [ ! -d "$features_loc_path" ]; then
    echo "Error: Directory '$features_loc_path' does not exist."
    exit 1
fi

echo "All paths and file exist. Proceeding with the script."

R_SCRIPT="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/wrapperQTL/inst/aux/matrixEQTLrun.R"

process_file() {

    local input=$1
    local output=$1 # fix!!!!!!!!

    echo "$input"
    echo "$output"

    covFilePath="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_05_002_20chk_ATAC_all_clusters/${input}/merged_covariates.txt"

    # covFilePath="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/real_data/results_log_norm_08_0_4chk/${input}/merged_covariates.txt"

    # Check if file exists
    if [ ! -f "$covFilePath" ]; then
        echo "Error: File '$covFilePath' does not exist."
        exit 1
    fi

    export input output R_SCRIPT snp_path features_path snp_loc_path features_loc_path covFilePath threads  # Default to 1 threads if not set

    export R_LIBS_USER="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/seurat_gal7/lib"

    $RUN_COMMAND -J "$JOB_NAME" -p "$PARTITION" -n "$TASKS" -t "$TIME" -N "$NODES" -c "$CPUS" -d "$DRY" -- 'docker_sing.sh -B /cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/ -H /cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/singularity_sandboxes/home_rsc -s -c -d drowsygoat/css Rscript $R_SCRIPT $input $output $snp_path $snp_loc_path $features_path $features_loc_path $covFilePath --threads $threads --verbose --colNames_convention "2" "_x_"'
}

input_directory=${1:-.}  # Defaults to current directory if no argument is given

output_directory="${features_path}/${input_directory}"

process_file $input_directory $output_directory

# example usage: >module_R_MEQTL.sh group_C7_results (in this case the input directory is group_C7_results (must be in current directory) and the output will be in features_path/group_C7_results (subdirectory of group_C7_results in features_p