#!/bin/zsh
# Konvergierte (ii-a)-Messung lambda=11,13,15 mit Fortschritts-Logging.
# Parallelitaet harter Cap via PID-Tracking, damit ABC-HCT-Jobs nicht clobbern.
set -u
cd ~/compute/rh_l2strip
source ~/.venvs/science/bin/activate

LOG=progress.log
SLOTS=5

JOBS=(
  "j01_l11_d440_sup:python xihat_strip_audit.py  --lambda 11 --N 145 --dps 440 --R 20 --h 0.02"
  "j02_l11_d440_l2:python xihat_l2strip_audit.py --lambda 11 --N 145 --dps 440"
  "j03_l11_d560_sup:python xihat_strip_audit.py  --lambda 11 --N 145 --dps 560 --R 20 --h 0.02"
  "j04_l11_d560_l2:python xihat_l2strip_audit.py --lambda 11 --N 145 --dps 560"
  "j05_l13_d540_sup:python xihat_strip_audit.py  --lambda 13 --N 170 --dps 540 --R 20 --h 0.02"
  "j06_l13_d540_l2:python xihat_l2strip_audit.py --lambda 13 --N 170 --dps 540"
  "j07_l13_d680_sup:python xihat_strip_audit.py  --lambda 13 --N 170 --dps 680 --R 20 --h 0.02"
  "j08_l13_d680_l2:python xihat_l2strip_audit.py --lambda 13 --N 170 --dps 680"
  "j09_l15_d640_sup:python xihat_strip_audit.py  --lambda 15 --N 195 --dps 640 --R 20 --h 0.02"
  "j10_l15_d640_l2:python xihat_l2strip_audit.py --lambda 15 --N 195 --dps 640"
  "j11_l15_d800_sup:python xihat_strip_audit.py  --lambda 15 --N 195 --dps 800 --R 20 --h 0.02"
  "j12_l15_d800_l2:python xihat_l2strip_audit.py --lambda 15 --N 195 --dps 800"
)
TOTAL=${#JOBS[@]}

ts() { date '+%Y-%m-%d %H:%M:%S' ; }
log() { printf '[%s] %s\n' "$(ts)" "$1" >> "$LOG" ; }

printf '[%s] [INIT] %d jobs, parallelism cap %d, on %s\n' \
  "$(ts)" "$TOTAL" "$SLOTS" "$(hostname)" > "$LOG"
ps -axo pid,%cpu,%mem,etime,command | awk '$1==70420 || $1==60896 {printf "[INIT]    pid=%s cpu=%s mem=%s etime=%s\n", $1,$2,$3,$4}' >> "$LOG"
echo >> "$LOG"

run_one() {
  local tag="$1"; local cmd="$2"
  local out="${tag}.txt"
  local t0=$(date +%s)
  log "[START] $tag"
  eval "$cmd" > "$out" 2>&1
  local rc=$?
  local dt=$(($(date +%s) - t0))
  if [ $rc -eq 0 ]; then
    log "[DONE]  $tag elapsed=${dt}s"
  else
    log "[FAIL]  $tag rc=$rc elapsed=${dt}s"
  fi
}

# PID-Cap: hartes Limit ueber kill -0 Polling
declare -a PIDS=()
prune() {
  local new=()
  for p in "${PIDS[@]}"; do
    if kill -0 "$p" 2>/dev/null; then new+=("$p"); fi
  done
  PIDS=("${new[@]}")
}

idx=0
for jobspec in "${JOBS[@]}"; do
  prune
  while [ ${#PIDS[@]} -ge $SLOTS ]; do
    sleep 5
    prune
  done
  idx=$((idx + 1))
  tag="${jobspec%%:*}"
  cmd="${jobspec#*:}"
  log "[QUEUE $idx/$TOTAL] $tag  (running: ${#PIDS[@]})"
  run_one "$tag" "$cmd" &
  PIDS+=("$!")
done
wait
log "[ALL_DONE] $TOTAL jobs finished"
