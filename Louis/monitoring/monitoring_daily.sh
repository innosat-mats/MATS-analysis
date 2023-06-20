#!/bin/bash

#### DATA AND HEALTH MONITORING ####

# activate conda environment ?

source /home/louis/miniconda3/etc/profile.d/conda.sh
conda activate mats-analysis


#### PARAMETERS ####

# main path
MATS_dir='/home/louis/MATS/'

# monitoring parameters (False by default)
sampling_period=60 # sampling period in seconds
data_processing=True # data generation and processing success rate between levels
mode=daily # monitoring mode
show_plots=False # shows interactive matplotlib plots (might be buggy)



#### MONITORING ROUTINE ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

# monitoring the day before
start_time=$(date -d "$date -1 days" +"%Y:%m:%d_00:00:00")
stop_time=$(date +"%Y:%m:%d_00:00:00")
day_start=$(date -d "$date -1 days" +'%Y_%m_%d')

# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/daily_'${day_start}

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode} --show_plots ${show_plots}; } &
pid=$!
wait $pid

# monitoring 2 days before
start_time=$(date -d "$date -2 days" +"%Y:%m:%d_00:00:00")
stop_time=$(date +"%Y:%m:%d_00:00:00")
day_start=$(date -d "$date -2 days" +'%Y_%m_%d')

# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/daily_'${day_start}

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode} --show_plots ${show_plots}; } &
pid=$!
wait $pid

echo 'End of program .....'
