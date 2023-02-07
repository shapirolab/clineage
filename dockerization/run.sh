#!/bin/bash
set -euo pipefail

docker run --rm -p 8888:80 -u root --entrypoint /bin/bash -it --network=host clineage:0.0.3
