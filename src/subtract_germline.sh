#!/usr/bin/env bash

set -e
set -xv

# wrapper for msi_measure to find matching normal or write dummy output for normals
TUMOUR="$1"

tfn=$(basename $TUMOUR)
TUMOUR_ID=${tfn/\.*/}
NORMAL_ID=$(src/find_normal.py $TUMOUR_ID < cfg/sample-metadata.csv)
NORMAL=${TUMOUR/$TUMOUR_ID/$NORMAL_ID}

if [ "$NORMAL_ID" = "NOTFOUND" ]; then
  # problem
  exit 1
else
  src/subtract_vcf.py --subtract $NORMAL < $TUMOUR
fi

