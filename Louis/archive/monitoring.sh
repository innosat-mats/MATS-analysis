#!/bin/bash

#### DATA AND HEALTH MONITORING ####

# activate conda environment ?

source /home/louis/miniconda3/etc/profile.d/conda.sh
conda activate mats-analysis


#### PARAMETERS ####

# main path
MATS_dir='/home/louis/MATS/'

# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/daily_monitoring_2023_06_17'

# monitoring parameters (False by default)
sampling_period=60 # sampling period in seconds
data_processing=Tue # data generation and processing success rate between levels
CRBD=True # CRBD temperature data
HTR=Tre # Heat sensors data
PWR=Tre # Temperature, current and voltage data from the PWR module
CPRU=Tre # CPRU voltage data including power and overvoltage
schedule=Tue # Mode plot
show_plots=False # shows interactive matplotlib plots (might be buggy)

# start and stop of the monitoring
start_time='2023:06:17_0:0:0'
stop_time='2023:06:18_0:0:0'


#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --data_processing ${data_processing} --CRBD ${CRBD} --HTR ${HTR} --PWR ${PWR} --CPRU ${CPRU} --schedule ${schedule} --show_plots ${show_plots}; } &
pid=$!
wait $pid

echo 'End of program .....'
