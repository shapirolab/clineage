#!/usr/bin/env python3
import django
django.setup()
from distributed.cli.dworker import go
go()
