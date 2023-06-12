#!/bin/bash


#### DATA AND HEALTH MONITORING ####

# activate conda environment ?
# conda activate mats_analysis


#### PARAMETERS ####
 
# main path
MATS_dir='/home/louis/MATS/'

# output folder
outdir=${MATS_dir}'MATS-Data/Monitoring/test'

# monitoring parameters
sampling_period=3600 # sampling period in seconds
data_processing=true
CRBD=true
HTR=true
PWR=true
CPRU=true

# start and stop of the monitoring
start_time='2023:05:1_0:0:0'
stop_time='2023:06:1_0:0:0'


#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --data_processing=${data_processing} --CRBD=${CRBD} --HTR=${HTR} --PWR=${PWR} --CPRU=${CPRU} ; } &
pid=$!
wait $pid

echo 'End of program .....'
