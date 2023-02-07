#!/bin/bash
set -euo pipefail

DOCKER_BUILDKIT=1 docker build -t clineage:0.0.3 .
