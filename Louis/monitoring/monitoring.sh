#!/bin/bash

#### DATA AND HEALTH MONITORING ####

# activate conda environment ?

source /home/louis/miniconda3/etc/profile.d/conda.sh
conda activate mats-analysis


#### PARAMETERS ####

# main path
MATS_dir='/home/louis/MATS/'

# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/test_all'

# monitoring parameters (False by default)
sampling_period=60 # sampling period in seconds
show_plots=False # shows interactive matplotlib plots (might be buggy)
mode=all # the modes are the following :
# daily : daily summary
# temp : all temperature data
# all : all accessible monitoring data
# power : current and voltage data


# start and stop of the monitoring
start_time='2023:03:26_13:0:0'
stop_time='2023:03:26_13:10:0'


#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid

echo 'End of program .....'
