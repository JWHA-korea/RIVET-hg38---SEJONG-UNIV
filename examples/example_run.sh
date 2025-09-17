#!/usr/bin/env bash
set -euo pipefail
rivet-gs --disease "Hereditary breast and/or ovarian cancer syndrome"
--disease-label "breast cancer"
--outdir runs
echo "Done. Check ./runs/"
