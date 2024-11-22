#! /usr/bin/env bash
source /opt/conda/etc/profile.d/conda.sh
conda activate starscope_scATAC_env
exec "$@"
