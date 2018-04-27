#!/usr/bin/env bash

set -e
set -xv

# wrapper for msi_measure to find matching normal or write dummy output for normals
REGIONS="$1"
TUMOUR="$2"
SUMMARY="$3"

filename=$(basename $TUMOUR)
TUMOUR_ID=${filename/\.*/}
NORMAL_ID=$(src/find_normal.py $TUMOUR_ID < cfg/sample-metadata.csv)
NORMAL=${TUMOUR/$TUMOUR_ID/$NORMAL_ID}

if [ "$NORMAL_ID" = "NOTFOUND" ]; then
  # problem
  exit 1
else
  python src/msi_measure.py --tumour $TUMOUR --normal $NORMAL --verbose --regions $REGIONS --summary $SUMMARY
fi
