#!/bin/bash

# activate conda environment ?
conda activate MATS_analysis

# main path
MATS_dir='/home/louis/projects/MATS/'


# output
outdir=${MATS_dir}'animations/daily/'$dates'/'

# measurement date
start_time='2022:4:24_0:0:0'
stop_time='2022:4:27_0:0:0

# initiate
echo -e "DAILY ANIMATION (${dates}): Initiating ..."
echo 'MATS DATA: L'${level}' data (v'${version}')..'
echo 'Generating figures in: '${outdir}

# generate figures
{ python ${MATS_dir}MATS-analysis/Louis/monitoring/monitoring_routine_weekly.py --outdir ${outdir} --start_time ${start_time} --stop_time ${stop_time}; } &
pid=$!
wait $pid

echo 'End of program .....'