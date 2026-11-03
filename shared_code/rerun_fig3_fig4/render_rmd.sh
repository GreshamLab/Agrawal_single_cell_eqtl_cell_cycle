#!/bin/bash
# Knit one Rmd in its own directory; log next to it as <Rmd>.render_<TAG>.log.
# Usage: render_rmd.sh <path/to/file.Rmd>   (returns the render exit code)
TAG=2026-11-01_arial
rmd_path="$1"
dir=$(dirname "$rmd_path"); rmd=$(basename "$rmd_path")
cd "$dir" || exit 1
log="${rmd%.Rmd}.render_${TAG}.log"
echo "$(date '+%F %T')  START  $rmd_path"
Rscript -e "rmarkdown::render('${rmd}', knit_root_dir = getwd())" > "$log" 2>&1
rc=$?
echo "$(date '+%F %T')  $([ $rc -eq 0 ] && echo DONE || echo "FAILED (exit $rc, see $dir/$log)")  $rmd_path"
exit $rc
