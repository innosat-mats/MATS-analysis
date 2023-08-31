#!/bin/bash

#### DATA AND HEALTH MONITORING ####

# activate conda environment ?

source /Users/lindamegner/miniconda3/etc/profile.d/conda.sh
conda activate mats-analysis


#### PARAMETERS ####

# main path
MATS_dir='/Users/lindamegner/MATS/MATS-retrieval/'

# current date
day_start=$(date -v -1d +'%Y_%m_%d')
# output folder for images
outdir=${MATS_dir}'MATS-data/daily_'${day_start}

# monitoring parameters (False by default)
sampling_period=60 # sampling period in seconds
mode=daily
show_plots=False # shows interactive matplotlib plots (might be buggy)

# start and stop of the monitoring
start_time=$(date -v -1d +"%Y:%m:%d_00:00:00")
stop_time=$(date +"%Y:%m:%d_00:00:00")


#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid

echo 'End of program .....'
