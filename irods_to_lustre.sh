#!/bin/bash

module load HGI/pipelines/irods_to_lustre
WORK_DIR=/lustre/scratch126/humgen/teams/hurles/users/vb9/Results_UTR_novaseq/irods_to_lustre
irods_to_lustre \
    -w $WORK_DIR/work \
    --run_mode "csv_samples_id" \
    --input_samples_csv $WORK_DIR/input/samples.tsv \
    --samples_to_process -1        \
    --outdir $WORK_DIR/output 
