#!/bin/bash
set -e

bash /init_once.sh

. $(which activate) py36

cd /clineage
./manage.py runserver 0.0.0.0:8000
