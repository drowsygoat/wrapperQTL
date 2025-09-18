#!/bin/bash

# MODULE PARAMETERS
RUN_COMMAND="run_shell_command.sh"
JOB_NAME="downstream_MEQTL"
PARTITION="main"
NODES=1
TIME="23:05:00"
TASKS=1
CPUS=32
DRY="with_eval"

PREFIX="defaultprefix"
p_value_threshold="NULL"

R_SCRIPT="/cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/QTL_analysis/wrapperQTL/inst/aux/processMEQTLdata.R"

process_file() {
    local input=$1
    local output=$2

    # Example modification: copy the file content (you can replace this with actual processing)

    echo "$input"
    echo "$output"

    export input output PREFIX p_value_threshold R_SCRIPT

    $RUN_COMMAND -J "$JOB_NAME" -p "$PARTITION" -n "$TASKS" -t "$TIME" -N "$NODES" -c "$CPUS" -d "$DRY" -- 'docker_sing.sh -B /cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/RR/scAnalysis/ -H /cfs/klemming/projects/supr/sllstore2017078/kaczma-workingdir/singularity_sandboxes/home_rsc -s -c -d drowsygoat/rsc Rscript $R_SCRIPT $input $PREFIX $p_value_threshold'
}

input_directory=${1:-.}

output_directory=${2:-$input_directory}

gather_ntests.R $input_directory ntests\\.txt $PREFIX

process_file $input_directory