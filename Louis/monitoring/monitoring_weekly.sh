#!/bin/bash

#### DATA AND HEALTH MONITORING ####

# activate conda environment ?

source /home/louis/miniconda3/etc/profile.d/conda.sh
conda activate mats-analysis


#### PARAMETERS ####

# main path
MATS_dir='/home/louis/MATS/'

# current date
week_start=$(date -d "$date -"$(date +'%w')" days - 6 days" +'%Y_%m_%d')

# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/weekly_'${week_start}

# monitoring parameters (False by default)
sampling_period=600 # sampling period in seconds
data_processing=True # data generation and processing success rate between levels
CRBD=True # CRBD temperature data
HTR=True # Heat sensors data
PWR=True # Temperature, current and voltage data from the PWR module
CPRU=True # CPRU voltage data including power and overvoltage
schedule=True # Mode plot
show_plots=False # shows interactive matplotlib plots (might be buggy)

# start and stop of the monitoring
start_time=$(date -d "$date -"$(date +'%w')" days - 6 days" +'%Y:%m:%d_00:00:00')
stop_time=$(date -d "$date -"$(date +'%w')" days + 1 days" +'%Y:%m:%d_00:00:00')


#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --data_processing ${data_processing} --CRBD ${CRBD} --HTR ${HTR} --PWR ${PWR} --CPRU ${CPRU} --schedule ${schedule} --show_plots ${show_plots}; } &
pid=$!
wait $pid

echo 'End of program .....'
