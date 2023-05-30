#!/bin/bash

# activate conda environment ?
conda activate MATS_analysis

# main path
MATS_dir='/home/louis/MATS/'


# output
outdir=${MATS_dir}'MATS-Data/Monitoring'

#### Daily data and health monitoring

# initiate daily monitoring
echo -e "Daily data and health monitoring: Initiating ..."



daily_start='2023:05:07_0:0:0'
daily_stop='2023:05:08_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid


daily_start='2023:05:09_0:0:0'
daily_stop='2023:05:10_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid

daily_start='2023:05:10_0:0:0'
daily_stop='2023:05:11_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid

daily_start='2023:05:12_0:0:0'
daily_stop='2023:05:13_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid

daily_start='2023:05:13_0:0:0'
daily_stop='2023:05:14_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid

daily_start='2023:05:14_0:0:0'
daily_stop='2023:05:15_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid

daily_start='2023:05:15_0:0:0'
daily_stop='2023:05:16_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid

daily_start='2023:05:16_0:0:0'
daily_stop='2023:05:17_0:0:0'
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_daily.py --outdir ${outdir} --start_time ${daily_start} --stop_time ${daily_stop}; } &
pid=$!
wait $pid
# #### Weekly summary generation

# # weekly summary 
# weekly_start='2023:05:01_0:0:0'
# weekly_stop='2023:06:01_0:0:0'

# # initiate weekly monitoring
# echo -e "Weekly data and health monitoring summary: Initiating ..."

# # generate figures (fixed time range)
# # { python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_weekly.py --outdir ${outdir} --start_time ${weekly_start} --stop_time ${weekly_stop}; } &

# # generate figures (current week)
# #{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_weekly.py --outdir ${outdir}; } &

# # pid=$!
# # wait $pid

echo 'End of program .....'