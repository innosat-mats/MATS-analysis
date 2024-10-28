#!/bin/bash

#### DATA AND HEALTH MONITORING ####
# The following script is used to monitor the health and data of the MATS instrument.
# It analyzes the last week of data and generates plots as .png images

# LINES TO BE CHANGED BY USER #

# activate conda environment
source /Users/lindamegner/miniconda3/etc/profile.d/conda.sh
conda activate MatsAnalysis

# path to the MATS folder containing the repositories (MATS-analysis, calibration_data, MATS-data, MATS-L1-processing ...)

MATS_dir='/Users/lindamegner/MATS/MATS-retrieval/'
#### RUNNING SCRIPT #### (no need to change anything except if you use MAC OS or Windows)
# current date
week_start=$(date -d "$date -"$(date +'%w')" days - 6 days" +'%Y_%m_%d') #For Linux and Windows users
#week_start=$(date -v-$(($(date +'%w') + 6))d +'%Y_%m_%d') #For MAC OS users

# output folder for images
outdir=${MATS_dir}'data/Monitoring/weekly_'${week_start}

# monitoring parameters
sampling_period=600 # sampling period in seconds
mode=daily
show_plots=False # shows interactive matplotlib plots (might be buggy)

# start and stop of the monitoring
#start_time=$(date -d "$date -"$(date +'%w')" days - 6 days" +'%Y:%m:%d_00:00:00') #For Linux and Windows users
#stop_time=$(date -d "$date -"$(date +'%w')" days + 1 days" +'%Y:%m:%d_00:00:00') #For Linux and Windows users
start_time=$(date -v-$(($(date +'%w') + 6))d +"%Y:%m:%d_00:00:00") #For MAC OS users
stop_time=$(date -v-$(($(date +'%w') - 1))d +"%Y:%m:%d_00:00:00") #For MAC OS users


#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid

echo 'End of program .....'
