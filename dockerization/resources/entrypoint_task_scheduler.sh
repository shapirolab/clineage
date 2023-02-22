#!/bin/bash
set -e

. $(which activate) py36
dask-scheduler
