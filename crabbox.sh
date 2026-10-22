#!/usr/bin/env bash
#
# crabbox.sh — islo.dev sandbox harness for the T2T-CHM13 vs GRCh38 comparison.
#
#   ./crabbox.sh t2t      warm one sandbox with T2T-CHM13v2.0 (UCSC hs1.fa.gz), snapshot it,
#                         run scripts/compute.py on 24 forks of the snapshot (MAXP at a time),
#                         and reduce against the GRCh38 receipts -> data/compare.new.json
#   ./crabbox.sh grch38   recompute the 24 GRCh38 receipts (data/wg_warm_<chr>.json) the same
#                         way. Optional: without them, `t2t` fetches the published receipts
#                         from zozo123/genomics-sandboxes at a pinned commit.
#
# Both modes exit non-zero unless all 24 shard outputs exist. Neither writes data/compare.json.
#
# Provenance
#   The published data/compare.json came from the direct islo CLI path (commit 5357a12).
#   The OpenClaw Crabbox shard path (used when crabbox is installed and ISLO_API_KEY or
#   CRABBOX_ISLO_API_KEY is set) was added in 61d1488 and has not been re-run.
#
# Requirements
#   - islo (https://islo.dev), logged in: sandboxes, snapshot save/list
#   - optional: crabbox (https://github.com/openclaw/crabbox) and ISLO_API_KEY
#
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"
export ISLO_OUTPUT_FORMAT=plain
if [ -z "${ISLO_API_KEY:-}" ] && [ -n "${CRABBOX_ISLO_API_KEY:-}" ]; then
  export ISLO_API_KEY="$CRABBOX_ISLO_API_KEY"
fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
DATA="$ROOT/data"; mkdir -p "$DATA"
LOG="$ROOT/scripts/wg.log"

# Pins. The original run recorded no versions; the image is what python:3.12-slim
# resolved to on the run date (2026-06-14). numpy==2.4.6 is pinned in scripts/*warmup.sh.
IMG="docker.io/library/python:3.12.13-slim"
CRABBOX_VERSION="0.36.0"   # the crabbox flags below were checked against this version's --help
GRCH38_RECEIPTS="https://raw.githubusercontent.com/zozo123/genomics-sandboxes/79ce3e0671f8c176c0983c20829768e8c95054fd/data"

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
MAXP="${MAXP:-8}"                       # concurrent forks
FORKMEM="${FORKMEM:-8192}"              # chr1 (~250 Mb) needs headroom
TS="$(date +%s)"
SNAP=""                                 # set per mode

now(){ python3 -c 'import time;print("%.3f"%time.time())'; }
since(){ python3 -c "print(round($(now)-$1,2))"; }
say(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOG"; }
b64(){ base64 < "$1" | tr -d '\n'; }
have_crabbox_islo(){ command -v crabbox >/dev/null 2>&1 && [ -n "${ISLO_API_KEY:-}" ]; }
check_crabbox_version(){
  local v; v="$(crabbox --version 2>/dev/null | head -1 || true)"
  if [ "$v" != "$CRABBOX_VERSION" ]; then
    say "warning: crabbox ${v:-unknown}; flags were checked against $CRABBOX_VERSION"
  fi
}
lease_sandbox_name(){
  python3 - "$1" <<'PY'
import json, sys
doc=json.load(open(sys.argv[1]))
vals=[]
def walk(x):
    if isinstance(x, dict):
        for k,v in x.items():
            if k in ("sandbox","sandboxName","name","id","leaseID","lease_id") and isinstance(v,str):
                vals.append(v)
            walk(v)
    elif isinstance(x, list):
        for v in x: walk(v)
walk(doc)
for v in vals:
    if v.startswith("crabbox-"):
        print(v); raise SystemExit
for v in vals:
    if v.startswith("isb_"):
        print(v[4:]); raise SystemExit
for v in vals:
    print(v); raise SystemExit
raise SystemExit("could not find sandbox name in lease output")
PY
}

# Copy the kernel and warm-up scripts into a fresh sandbox (islo CLI path).
inject(){
  local w c t
  w="$(b64 scripts/warmup.sh)"; c="$(b64 scripts/compute.py)"; t="$(b64 scripts/t2t_warmup.sh)"
  islo use "$1" --image "$IMG" --cpu 2 --memory 4096 --disk 40 -- bash -lc \
    "mkdir -p \$HOME/gx && echo $w | base64 -d > \$HOME/gx/warmup.sh && echo $c | base64 -d > \$HOME/gx/compute.py && echo $t | base64 -d > \$HOME/gx/t2t_warmup.sh && chmod +x \$HOME/gx/warmup.sh \$HOME/gx/t2t_warmup.sh && echo ok" >>"$LOG" 2>&1
}

# Keep the kernel's JSON line from a shard's stdout. $1 raw output, $2 destination.
keep_json(){
  if grep -q '^{"chrom"' "$1"; then grep '^{"chrom"' "$1" > "$2"; else return 1; fi
}

# One fork per chromosome from snapshot $SNAP. $1 box name, $2 chromosome.
run_t2t_shard(){
  local box="$1" chr="$2" out="$DATA/_t2t_$2.out" try ts
  for try in 1 2; do
    islo rm "$box" >>"$LOG" 2>&1 || true
    ts=$(now)
    if have_crabbox_islo; then
      crabbox run --provider islo \
        --islo-image "$IMG" \
        --islo-snapshot-name "$SNAP" \
        --islo-vcpus 2 \
        --islo-memory-mb "$FORKMEM" \
        --islo-disk-gb 40 \
        --label "t2t-$chr" \
        -- bash -lc "GX_REF=\$HOME/gx/ref_t2t python3 scripts/compute.py $chr" \
        > "$out" 2>>"$LOG" || true
    else
      islo use "$box" --snapshot "$SNAP" --cpu 2 --memory "$FORKMEM" --disk 40 \
        -- bash -lc "GX_REF=\$HOME/gx/ref_t2t python3 \$HOME/gx/compute.py $chr" \
        > "$out" 2>>"$LOG" || true
    fi
    if keep_json "$out" "$DATA/t2t_$chr.json"; then
      since "$ts" > "$DATA/t2t_$chr.time"
      say "  [t2t] $chr $(cat "$DATA/t2t_$chr.time")s (try $try)"; islo rm "$box" >>"$LOG" 2>&1 || true; return 0
    fi
    say "  [t2t] $chr try $try: no JSON, retry"
  done
  say "  [t2t] $chr FAILED"; islo rm "$box" >>"$LOG" 2>&1 || true; return 1
}

run_grch38_shard(){ # $1 box name, $2 chromosome
  local box="$1" chr="$2" out="$DATA/_wg_warm_$2.out" try ts
  for try in 1 2; do
    islo rm "$box" >>"$LOG" 2>&1 || true
    ts=$(now)
    islo use "$box" --snapshot "$SNAP" --cpu 2 --memory "$FORKMEM" \
      -- bash -lc "python3 \$HOME/gx/compute.py $chr" > "$out" 2>>"$LOG" || true
    if keep_json "$out" "$DATA/wg_warm_$chr.json"; then
      since "$ts" > "$DATA/wg_warm_$chr.time"
      say "  [grch38] $chr $(cat "$DATA/wg_warm_$chr.time")s (try $try)"; islo rm "$box" >>"$LOG" 2>&1 || true; return 0
    fi
    say "  [grch38] $chr try $try: no JSON, retry"
  done
  say "  [grch38] $chr FAILED"; islo rm "$box" >>"$LOG" 2>&1 || true; return 1
}

# Run $1 (a shard function) for every chromosome, MAXP at a time. $2 box-name prefix.
fanout(){
  local fn="$1" prefix="$2" i=0 j p chr pids
  while [ "$i" -lt "${#CHROMS[@]}" ]; do
    pids=(); j="$i"
    while [ "$j" -lt "${#CHROMS[@]}" ] && [ "$j" -lt $((i + MAXP)) ]; do
      chr="${CHROMS[$j]}"
      "$fn" "$prefix-$chr" "$chr" & pids+=("$!")
      j=$((j + 1))
    done
    for p in "${pids[@]}"; do wait "$p" || true; done
    i=$((i + MAXP))
  done
}

# Exit 1 unless data/<prefix>_<chr>.json exists for all 24 chromosomes.
require_outputs(){
  local c missing=""
  for c in "${CHROMS[@]}"; do
    if [ ! -s "$DATA/$1_$c.json" ]; then missing="$missing $c"; fi
  done
  if [ -n "$missing" ]; then
    say "ERROR: no $1 output for:$missing. Not reducing."
    exit 1
  fi
}

# Use local GRCh38 receipts only if all 24 exist; otherwise fetch all 24 at the pinned commit.
grch38_receipts(){
  local c n=0
  for c in "${CHROMS[@]}"; do
    if [ -s "$DATA/wg_warm_$c.json" ]; then n=$((n + 1)); fi
  done
  if [ "$n" -eq "${#CHROMS[@]}" ]; then
    say "using local GRCh38 receipts (data/wg_warm_*.json)"; return 0
  fi
  say "fetching GRCh38 receipts from $GRCH38_RECEIPTS"
  for c in "${CHROMS[@]}"; do
    curl -fsSL "$GRCH38_RECEIPTS/wg_warm_$c.json" -o "$DATA/wg_warm_$c.json"
  done
}

cmd_t2t(){
  local warm_box="gx-t2t-warm" via t warm snaps snapsz map
  : > "$LOG"
  SNAP="genomics-t2t-$TS"
  if have_crabbox_islo; then
    via="OpenClaw Crabbox (islo provider)"; check_crabbox_version
  else
    via="islo CLI"
    say "crabbox or ISLO_API_KEY not available; using the islo CLI path (the one that produced the published numbers)"
  fi
  say "=== T2T-CHM13v2.0 vs GRCh38 · ${#CHROMS[@]} chromosomes · shards via $via ==="
  rm -f "$DATA"/t2t_chr*.json "$DATA"/t2t_chr*.time "$DATA"/_t2t_*.out "$DATA/compare.new.json"

  say "=== 1/4 warm base: download hs1, md5-check, split per chromosome ==="
  islo rm "$warm_box" >>"$LOG" 2>&1 || true
  t=$(now)
  if have_crabbox_islo; then
    crabbox run --provider islo \
      --islo-image "$IMG" \
      --islo-vcpus 2 \
      --islo-memory-mb 4096 \
      --islo-disk-gb 40 \
      --keep \
      --lease-output "$DATA/t2t_warm_lease.json" \
      --label t2t-warm \
      -- bash -lc "bash scripts/t2t_warmup.sh && mkdir -p \$HOME/gx && cp scripts/compute.py \$HOME/gx/compute.py" \
      >>"$LOG" 2>&1
    warm_box="$(lease_sandbox_name "$DATA/t2t_warm_lease.json")"
    say "warm base kept by Crabbox as $warm_box"
  else
    inject "$warm_box"
    islo use "$warm_box" --cpu 2 --memory 4096 --disk 40 -- bash -lc "\$HOME/gx/t2t_warmup.sh" >>"$LOG" 2>&1
  fi
  warm=$(since "$t"); say "T2T warm base built in ${warm}s"

  say "=== 2/4 snapshot -> $SNAP ==="
  t=$(now); islo snapshot save "$warm_box" --name "$SNAP" >>"$LOG" 2>&1; snaps=$(since "$t")
  snapsz=$(islo snapshot ls 2>/dev/null | grep -F "$SNAP" | awk '{print $2" "$3}' || true)
  say "snapshot saved in ${snaps}s (size ${snapsz:-n/a})"
  if have_crabbox_islo; then
    crabbox stop --provider islo "$warm_box" >>"$LOG" 2>&1 || true
  else
    islo rm "$warm_box" >>"$LOG" 2>&1 || true
  fi

  say "=== 3/4 fan-out: ${#CHROMS[@]} shards, $MAXP concurrent ==="
  t=$(now); fanout run_t2t_shard gx-t2t; map=$(since "$t")
  say "T2T fan-out wall-clock ${map}s"
  require_outputs t2t

  say "=== 4/4 reduce against GRCh38 -> data/compare.new.json ==="
  grch38_receipts
  require_outputs wg_warm
  python3 scripts/reduce_compare.py "$DATA" "$SNAP" "$snapsz" "$warm" "$snaps" "$map" "${CHROMS[@]}" \
    --harness "crabbox.sh t2t, shards via $via" --out "$DATA/compare.new.json" | tee -a "$LOG"
  say "=== DONE (snapshot $SNAP retained). Diff data/compare.new.json against data/compare.json. ==="
}

cmd_grch38(){
  local warm_box="gx-wg-warm" t warm snaps map
  : > "$LOG"
  SNAP="genomics-wg-$TS"
  say "=== GRCh38 receipts · ${#CHROMS[@]} chromosomes · islo CLI ==="
  rm -f "$DATA"/wg_warm_chr*.json "$DATA"/wg_warm_chr*.time "$DATA"/_wg_warm_*.out

  say "=== 1/3 warm base: download ${#CHROMS[@]} hg38 chromosomes, md5-check ==="
  islo rm "$warm_box" >>"$LOG" 2>&1 || true
  t=$(now); inject "$warm_box"
  islo use "$warm_box" --cpu 2 --memory 4096 --disk 25 -- bash -lc "\$HOME/gx/warmup.sh ${CHROMS[*]}" >>"$LOG" 2>&1
  warm=$(since "$t"); say "warm base built in ${warm}s"

  say "=== 2/3 snapshot -> $SNAP ==="
  t=$(now); islo snapshot save "$warm_box" --name "$SNAP" >>"$LOG" 2>&1; snaps=$(since "$t")
  say "snapshot saved in ${snaps}s"
  islo rm "$warm_box" >>"$LOG" 2>&1 || true

  say "=== 3/3 fan-out: ${#CHROMS[@]} shards, $MAXP concurrent ==="
  t=$(now); fanout run_grch38_shard gx-wg; map=$(since "$t")
  say "GRCh38 fan-out wall-clock ${map}s"
  require_outputs wg_warm
  say "=== DONE: data/wg_warm_*.json written (snapshot $SNAP retained); ./crabbox.sh t2t reduces against them ==="
}

case "${1:-}" in
  t2t) cmd_t2t ;;
  grch38) cmd_grch38 ;;
  *) echo "usage: $0 {t2t|grch38}" >&2; exit 1 ;;
esac
