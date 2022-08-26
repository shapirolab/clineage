#!/bin/bash
set -euo pipefail

docker run --rm -u root --entrypoint /usr/bin/fish -it clineage 
