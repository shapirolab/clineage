#!/bin/bash
source ~/.virtualenvs/cl/bin/activate
cd ~/CLineage
exec gunicorn -p /home/ofirr/clineage.pid -b 127.0.0.1:9000 -w 2 clineage.wsgi:application
