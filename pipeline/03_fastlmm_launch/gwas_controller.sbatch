#!/bin/bash
#SBATCH --job-name=gwas_controller
#SBATCH --account=torch_pr_130_general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all/gwas_queue/logs/controller_%j.log

# GWAS queue controller.
#
# Keeps MAX_ACTIVE variant arrays (condition x variant from queue.tsv) running at
# all times. Every INTERVAL it:
#   1. resubmits itself (so it keeps going without an open session),
#   2. marks items done when every gene in their batches has a .threshold.txt,
#   3. submits the next not-done items, in queue.tsv order, until MAX_ACTIVE run.
# An item is "running" while any job named GWAS_<variant> whose working
# directory is that condition's batch_scripts dir is in the queue (this includes
# tasks the timeout handler resubmitted). An item that stops with genes missing
# is resubmitted (finished genes are skipped) up to MAX_SUBMITS times in total.
#
# Stop it:   touch  .../gwas_queue/STOP        (remove the file and sbatch again to restart)
# Progress:  cat    .../gwas_queue/status.txt  and  .../gwas_queue/controller.log
# Dry run:   DRY_RUN=1 bash gwas_controller.sbatch   (prints actions, submits nothing)

set -uo pipefail

ROOT=/projects/rps/cgsb/gresham/Akriti/wildyeast1/common_analysis_for_all
Q=$ROOT/gwas_queue
QUEUE=$Q/queue.tsv
STATE=$Q/state            # <condition>__<variant>.submits  /  .done  /  .gave_up
LOG=$Q/controller.log
STATUS=$Q/status.txt
SELF=$Q/gwas_controller.sbatch
MAX_ACTIVE=${MAX_ACTIVE:-4}
MAX_SUBMITS=${MAX_SUBMITS:-3}
INTERVAL=${INTERVAL:-10}  # minutes
DRY_RUN=${DRY_RUN:-0}

mkdir -p "$STATE" "$Q/logs"
log() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }

# only one controller at a time
exec 9>"$Q/.lock"
flock -n 9 || { log "another controller run holds the lock, exiting"; exit 0; }

if [ -f "$Q/STOP" ]; then
    log "STOP file present, not running and not rescheduling"
    exit 0
fi

# 1. reschedule first, so a failure below never breaks the chain
if [ "$DRY_RUN" != 1 ]; then
    others=$(squeue -h -u "$USER" -n gwas_controller -o "%i" | grep -v -x "${SLURM_JOB_ID:-none}" | wc -l)
    if [ "$others" -eq 0 ]; then
        sbatch --begin=now+${INTERVAL}minutes "$SELF" > /dev/null \
            || log "WARNING: could not reschedule controller"
    fi
fi

# running (job name, working dir) pairs, one per line: "GWAS_SNPs<TAB>/path"
running=$(squeue -h -u "$USER" -o "%j	%Z" | sort -u)

n_active=0
n_done=0
declare -a todo=()
status_lines=()

while IFS=$'\t' read -r cond var rel; do
    [[ -z "$cond" || "$cond" == \#* ]] && continue
    dir=$ROOT/$rel
    key=${cond}__${var}
    submits=$(cat "$STATE/$key.submits" 2>/dev/null || echo 0)

    if [ -f "$STATE/$key.done" ]; then
        n_done=$((n_done + 1)); status_lines+=("$cond	$var	done"); continue
    fi
    if [ -f "$STATE/$key.gave_up" ]; then
        status_lines+=("$cond	$var	gave_up ($(cat "$STATE/$key.gave_up"))"); continue
    fi

    if grep -q -x -F "GWAS_${var}	${dir}" <<< "$running"; then
        n_active=$((n_active + 1)); status_lines+=("$cond	$var	running (submits: $submits)"); continue
    fi

    # not running: count finished genes
    outdir=$(grep -m1 '^OUTDIR=' "$dir/run_task_${var}.bash" | cut -d'"' -f2)
    batchdir=$(grep -m1 '^BATCH_DIR=' "$dir/run_task_${var}.bash" | cut -d'"' -f2)
    expected=$(cat "$batchdir"/phenotype_batch_*.txt | grep -c .)
    finished=0
    [ -d "$outdir" ] && finished=$(find "$outdir" -mindepth 2 -maxdepth 2 -name '*.threshold.txt' | wc -l)

    if [ "$finished" -ge "$expected" ]; then
        [ "$DRY_RUN" != 1 ] && touch "$STATE/$key.done"
        log "DONE: $cond $var ($finished/$expected genes)"
        n_done=$((n_done + 1)); status_lines+=("$cond	$var	done"); continue
    fi

    if [ "$submits" -ge "$MAX_SUBMITS" ]; then
        [ "$DRY_RUN" != 1 ] && echo "$finished/$expected genes after $submits submissions" > "$STATE/$key.gave_up"
        log "GAVE UP: $cond $var ($finished/$expected genes after $submits submissions)"
        status_lines+=("$cond	$var	gave_up"); continue
    fi

    todo+=("$cond	$var	$dir	$finished	$expected")
    status_lines+=("$cond	$var	waiting ($finished/$expected genes, submits: $submits)")
done < "$QUEUE"

# 2. fill free slots in queue order
slots=$((MAX_ACTIVE - n_active))
for item in "${todo[@]}"; do
    [ "$slots" -le 0 ] && break
    IFS=$'\t' read -r cond var dir finished expected <<< "$item"
    key=${cond}__${var}
    if [ "$DRY_RUN" = 1 ]; then
        log "DRY RUN would submit: $cond $var ($finished/$expected genes done)"
        slots=$((slots - 1)); continue
    fi
    if out=$(cd "$dir" && sbatch "submit_gwas_${var}.sh" 2>&1); then
        submits=$(( $(cat "$STATE/$key.submits" 2>/dev/null || echo 0) + 1 ))
        echo "$submits" > "$STATE/$key.submits"
        log "SUBMITTED: $cond $var -> ${out##* } (submission $submits, $finished/$expected genes already done)"
        slots=$((slots - 1))
    else
        log "sbatch failed for $cond $var, will retry next run: $out"
        break
    fi
done

{
    echo "Updated: $(date '+%F %T')   active: $n_active   done: $n_done / ${#status_lines[@]}"
    printf '%s\n' "${status_lines[@]}" | column -t -s $'\t'
} > "$STATUS"

if [ "$n_active" -eq 0 ] && [ "${#todo[@]}" -eq 0 ] && [ "$DRY_RUN" != 1 ]; then
    log "QUEUE FINISHED ($n_done done, $(( ${#status_lines[@]} - n_done )) gave up) - stopping controller"
    touch "$Q/STOP"
    scancel -u "$USER" -n gwas_controller --state=PENDING 2>/dev/null
fi
