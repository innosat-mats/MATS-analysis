#!/bin/bash

# activate conda environment ?
conda activate MATS_analysis

# main path
MATS_dir='/home/louis/MATS/'


# output
outdir=${MATS_dir}'MATS-Data/Monitoring'

#### Daily data and health monitoring

# daily monitoring
daily_start='2023:04:24_0:0:0'
daily_stop='2023:04:24_0:10:0'

# initiate daily monitoring
echo -e "Daily data and health monitoring: Initiating ..."

# generate figures
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid


#### Weekly summary generation

# weekly summary 
weekly_start='2023:04:01_0:0:0'
weekly_stop='2023:05:01_0:0:0'

# initiate weekly monitoring
echo -e "Weekly data and health monitoring summary: Initiating ..."

# generate figures
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_weekly.py --outdir ${outdir} --start_time ${weekly_start} --stop_time ${weekly_stop}; } &
pid=$!
wait $pid

echo 'End of program .....'