#!/bin/bash

#### DATA AND HEALTH MONITORING ####
# The following script is used to monitor the health and data of the MATS instrument.
# It is run on a specific date range and generates plots as .png images



# LINES TO BE CHANGED BY USER #

# activate conda environment
source /home/louis/miniconda3/etc/profile.d/conda.sh # path to conda installation
conda activate mats-analysis

# path to the MATS folder containing the repositories (MATS-analysis, calibration_data, MATS-data, MATS-L1-processing ...)
MATS_dir='/home/louis/MATS/'

# output folder for the generated images
outdir=${MATS_dir}'MATS-Data/Monitoring/schedule'

# monitoring parameters
sampling_period=120 # sampling period in seconds
show_plots=False # shows interactive matplotlib plots (might be buggy)
mode=daily # the modes are the following :
# daily : daily summary
# temp : all temperature data
# all : all accessible monitoring data
# power : current and voltage data
# custom : custom mode defined directly in the monitoring_routine.py file

# start and stop of the monitoring
start_time='2023:07:26_00:00:00'
stop_time='2023:07:29_00:00:00'


#### RUNNING SCRIPT ####

# a python script is run in the background, if your MATS-analysis repository is up to date, nothing has to be changed here

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid

echo 'End of program .....'
