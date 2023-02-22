#!/bin/bash
set -e

. $(which activate) py36
PYTHONPATH='/clineage' DJANGO_SETTINGS_MODULE='clineage.settings' dask-worker --nthreads 1 --memory-limit 300000000000 --no-bokeh --preload /clineage/worker-setup.py task-scheduler:8786
