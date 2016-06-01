#!/bin/bash

#$ -S /bin/bash
#$ -o $HOME/cluster_output/o.$JOB_ID.$TASK_ID.$USER.$HOSTNAME.stdout
#$ -e $HOME/cluster_output/o.$JOB_ID.$TASK_ID.$USER.$HOSTNAME.stderr
#$ -cwd
#$ -N dworker.$JOB_ID
echo $JOB_ID
echo $TASK_ID
echo $USER
echo $HOSTNAME
source /home/dcsoft/.venvs/cl/bin/activate
PYTHONPATH='/home/dcsoft/clineage' DJANGO_SETTINGS_MODULE='clineage.settings' ./clworker.py --nthreads 1 132.76.81.158:8786
