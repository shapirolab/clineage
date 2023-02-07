#!/bin/bash
set -e

. $(which activate) py36
jupyter notebook --allow-root --ip=0.0.0.0
