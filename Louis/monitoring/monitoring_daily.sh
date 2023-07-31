#!/bin/bash

#### DATA AND HEALTH MONITORING ####
# The following script is used to monitor the health and data of the MATS instrument.
# It analyzes the last day of data and generates plots as .png images

# LINES TO BE CHANGED BY USER #

# activate conda environment
source /home/louis/miniconda3/etc/profile.d/conda.sh # path to conda installation
conda activate mats-analysis

# path to the MATS folder containing the repositories (MATS-analysis, calibration_data, MATS-data, MATS-L1-processing ...)
MATS_dir='/home/louis/MATS/'

#### RUNNING SCRIPT #### (no need to change anything except if you use MAC OS or Windows)
# current date
day_start=$(date -d "$date -1 days" +'%Y_%m_%d') #For Linux and Windows users
#day_start=$(date -v -1d +'%Y_%m_%d') #For MAC OS users

# output folder for the generated images
outdir=${MATS_dir}'MATS-Data/Monitoring/daily_'${day_start}

# monitoring parameters
sampling_period=60 # sampling period in seconds
mode=daily
show_plots=False # shows interactive matplotlib plots (might be buggy)

# start and stop of the monitoring
start_time=$(date -d "$date -1 days" +"%Y:%m:%d_00:00:00") #For Linux and Windows users
#start_time=$(date -v -1d +"%Y:%m:%d_00:00:00") #For MAC OS users
stop_time=$(date +"%Y:%m:%d_00:00:00")


#### RUNNING  PYTHON SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid

echo 'End of program .....'
