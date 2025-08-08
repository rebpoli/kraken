#!/bin/bash
#$ -m ea
#$ -pe orte 8   # orte="fill-first", orte_rr="round-robin"
#$ -cwd
#$ -j y
#$ -N ipars
#$ -S /bin/bash

msg="`date`: START  job=${JOB_ID} pwd=`pwd` system=${SYSTEM} cores=${NSLOTS}"
echo "$msg"
echo "$msg" >> ~/${SYSTEM}.log

STARTTIME=$(date +%s)

mpirun -np $NSLOTS ./ipars
STATUS=$?

ENDTIME=$(date +%s)
elap=$((ENDTIME - STARTTIME))
ELAPSED=`printf "$((elap / 86400))d";date -d "0 $elap sec" +"%Hh%Mm"`

msg="`date`: FINISH job=${JOB_ID} pwd=`pwd` system=${SYSTEM} cores=${NSLOTS} elapsed=${ELAPSED} status=${STATUS}"
echo "$msg"
echo "$msg" >> ~/${SYSTEM}.log

