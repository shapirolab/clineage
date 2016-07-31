#!/bin/bash

#$ -S /bin/bash
#$ -o $HOME/cluster_out/o.$JOB_ID.$TASK_ID.$USER.stdout
#$ -e $HOME/cluster_out/o.$JOB_ID.$TASK_ID.$USER.stderr
#$ -cwd
echo $HOSTNAME
export LD_LIBRARY_PATH=/usr/lib64:/usr/lib64/atlas/:$LD_LIBRARY_PATH
/usr/wisdom/python/bin/python2.7 dyn_calibration_fc.py -i $JOB_ID -n $SGE_TASK_ID -u $USER --c1_min 5 --c1_max 15 --c2_min 20 --c2_max 60