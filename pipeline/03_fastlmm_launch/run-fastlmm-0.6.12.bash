#!/bin/bash

__base_dir=$(realpath -- $(dirname -- "${BASH_SOURCE[0]}"))

args=
for i in "$@"; do 
  i="${i//\\/\\\\}"
  args="${args} \"${i//\"/\\\"}\""
done

if [[ "${args}" == "" ]]; then args="/bin/bash"; fi

if [[ "${SINGULARITY_CONTAINER}" != "" ]]; then 
  export PATH=/share/apps/apptainer/bin:${PATH}
fi

singularity exec ${nv} \
--overlay ${__base_dir}/fastlmm-0.6.12.sqf:ro \
${__base_dir}/ubuntu-24.04.2.sif \
/bin/bash -c "
unset -f which
if [[ -e /ext3/env.sh ]]; then source /ext3/env.sh; fi
${args}
"
