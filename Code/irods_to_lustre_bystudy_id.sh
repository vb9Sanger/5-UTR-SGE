#!/bin/bash

module load HGI/pipelines/irods_to_lustre
WORK_DIR=/lustre/scratch126/humgen/teams/hurles/users/vb9/Results_UTR_novaseq/irods_to_lustre
irods_to_lustre \
    -w $WORK_DIR/work \
    --run_mode "study_id" \
    --input_studies 7885 \
    --samples_to_process -1        \
    --crams_to_fastq_min_reads 10 \
    --run_imeta_study true        \
    --run_merge_crams false \
    --run_crams_to_fastq false     \
    --outdir $WORK_DIR/output 
