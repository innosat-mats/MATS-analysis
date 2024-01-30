#!/bin/bash

#### DATA AND HEALTH MONITORING ####

# shell script to chack MATS health data during interesting events

# activate conda environment


source /Users/lindamegner/miniconda3/etc/profile.d/conda.sh
conda activate mats-analysis

#### PARAMETERS ####

# main path

MATS_dir='/Users/lindamegner/MATS/MATS-retrieval/'


# monitoring parameters (False by default)
sampling_period=600 # sampling period in seconds
show_plots=False # shows interactive matplotlib plots (might be buggy)
mode=all # the modes are the following :
# daily : daily summary
# temp : all temperature data
# all : all accessible monitoring data
# power : current and voltage data





#### RUNNING SCRIPT ####

# initiate monitoring
echo -e "Data and health monitoring: Initiating ..."

'''
# start and stop of the monitoring
start_time='2023:02:17_0:0:0'
stop_time='2023:02:19_0:0:0'
mode=power
# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/event_analysis/02_17_02_19/'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid
'''

# start and stop of the monitoring
start_time='2023:12:21_0:0:0'
stop_time='2023:12:24_0:0:0'
mode=all
# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/event_analysis/03_25_03_29/'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid

'''

# start and stop of the monitoring
start_time='2023:04:06_0:0:0'
stop_time='2023:04:07_0:0:0'
mode=power
# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/event_analysis/04_06_04_07/'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid



# start and stop of the monitoring
start_time='2023:04:19_0:0:0'
stop_time='2023:04:23_0:0:0'
# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/event_analysis/04_19_04_23/'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid



# start and stop of the monitoring
start_time='2023:05:03_0:0:0'
stop_time='2023:05:18_0:0:0'
mode=all
sampling_period=600
# output folder for images
outdir=${MATS_dir}'MATS-Data/Monitoring/event_analysis/05_03_05_18/'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time} --sampling_period ${sampling_period} --mode ${mode}; } &
pid=$!
wait $pid



'''

echo 'End of program .....'
