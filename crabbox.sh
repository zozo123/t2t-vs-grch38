#!/usr/bin/env bash
#
# crabbox.sh — genome-wide reference-broadcast fan-out with OpenClaw Crabbox.
#
#   ./crabbox.sh run        # warm one box (all 24 GRCh38 chromosomes + index),
#                           # snapshot it, fork one box per chromosome, reduce, time it, tear down.
#
# Why this exists
#   The reference genome + indices are large, read-only, and shared by every shard
#   of a scatter-gather. Instead of re-staging them per worker, we warm one box,
#   snapshot its initialized state once, and restore-and-fork it copy-on-write per
#   chromosome. This is the genome-wide version of the 5-chromosome demo.
#
# How it works
#   OpenClaw Crabbox (https://github.com/openclaw/crabbox) is the ephemeral "warm a box,
#   sync the working tree, run the suite, tear down" path; it leases boxes from a
#   broker or from a provider. For islo.dev, use `crabbox run --provider islo`
#   with ISLO_API_KEY exported. Snapshot creation is still done with the islo CLI
#   because it is the snapshot primitive the paper is testing; shard execution uses
#   Crabbox when the islo provider is authenticated.
#
# Requirements
#   - crabbox         (https://github.com/openclaw/crabbox)
#   - ISLO_API_KEY    required for `crabbox run --provider islo`
#   - islo            (https://islo.dev), logged in, for snapshot save/list
#
set -uo pipefail
export PATH="$HOME/.local/bin:$PATH"
export ISLO_OUTPUT_FORMAT=plain
if [ -z "${ISLO_API_KEY:-}" ] && [ -n "${CRABBOX_ISLO_API_KEY:-}" ]; then
  export ISLO_API_KEY="$CRABBOX_ISLO_API_KEY"
fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; cd "$ROOT"
DATA="$ROOT/data"; mkdir -p "$DATA"
LOG="$ROOT/scripts/wg.log"; : > "$LOG"

# all 24 GRCh38 chromosomes
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
        chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
MAXP="${MAXP:-8}"                       # concurrent forks (fabric-friendly waves)
TS="$(python3 -c 'import time;print(int(time.time()))')"
SNAP="genomics-wg-$TS"
IMG="docker.io/library/python:3.12-slim"
FORKMEM="${FORKMEM:-8192}"              # chr1 (~250 Mb) needs headroom
COLD_REF="chr1"                         # cold baseline shard (largest)

now(){ python3 -c 'import time;print("%.3f"%time.time())'; }
say(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOG"; }
have_crabbox_islo(){ command -v crabbox >/dev/null 2>&1 && [ -n "${ISLO_API_KEY:-}" ]; }
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
ensure_grch38_receipts(){
  local c url
  for c in "${CHROMS[@]}"; do
    [ -s "$DATA/wg_warm_$c.json" ] && continue
    url="https://raw.githubusercontent.com/zozo123/genomics-sandboxes/main/data/wg_warm_$c.json"
    say "fetching GRCh38 receipt for $c"
    curl -fsSL "$url" -o "$DATA/wg_warm_$c.json"
  done
}
WARMB64="$(base64 < "$ROOT/scripts/warmup.sh")"
COMPB64="$(base64 < "$ROOT/scripts/compute.py")"
T2TB64="$(base64 < "$ROOT/scripts/t2t_warmup.sh" 2>/dev/null || true)"
inject(){ islo use "$1" --image "$IMG" --cpu 2 --memory 4096 --disk 40 -- bash -lc \
  "mkdir -p \$HOME/gx && echo $WARMB64 | base64 -d > \$HOME/gx/warmup.sh && echo $COMPB64 | base64 -d > \$HOME/gx/compute.py && echo $T2TB64 | base64 -d > \$HOME/gx/t2t_warmup.sh && chmod +x \$HOME/gx/warmup.sh \$HOME/gx/t2t_warmup.sh && echo ok" >>"$LOG" 2>&1; }

# T2T fork: compute one chromosome from the T2T snapshot (GX_REF points at the split T2T ref)
run_t2t_shard(){ # $1 box $2 chr
  local box="$1" chr="$2" try ts el
  for try in 1 2; do
    islo rm "$box" >>"$LOG" 2>&1 || true
    ts=$(now)
    if command -v crabbox >/dev/null 2>&1 && [ -n "${ISLO_API_KEY:-}" ]; then
      crabbox run --provider islo \
        --islo-image "$IMG" \
        --islo-snapshot-name "$SNAP" \
        --islo-vcpus 2 \
        --islo-memory-mb "$FORKMEM" \
        --islo-disk-gb 40 \
        --label "t2t-$chr" \
        -- bash -lc "GX_REF=\$HOME/gx/ref_t2t python3 scripts/compute.py $chr" \
        > "$DATA/_t2t_$chr.out" 2>>"$LOG"
    else
      [ -z "${ISLO_API_KEY:-}" ] && say "  [t2t] OpenClaw crabbox islo provider needs ISLO_API_KEY; falling back to direct islo CLI for $chr"
      islo use "$box" --snapshot "$SNAP" --cpu 2 --memory "$FORKMEM" --disk 40 -- bash -lc "GX_REF=\$HOME/gx/ref_t2t python3 \$HOME/gx/compute.py $chr" > "$DATA/_t2t_$chr.out" 2>>"$LOG"
    fi
    el=$(python3 -c "print(round($(now)-$ts,2))")
    if grep -q '^{"chrom"' "$DATA/_t2t_$chr.out"; then
      grep '^{"chrom"' "$DATA/_t2t_$chr.out" > "$DATA/t2t_$chr.json"; echo "$el" > "$DATA/t2t_$chr.time"
      say "  [t2t] $chr ${el}s (try $try)"; islo rm "$box" >>"$LOG" 2>&1 || true; return 0
    fi
    say "  [t2t] $chr try $try: no JSON, retry"
  done
  say "  [t2t] $chr FAILED"; islo rm "$box" >>"$LOG" 2>&1 || true; return 1
}

run_shard(){ # $1 box $2 chr $3 mode(cold|warm)
  local box="$1" chr="$2" mode="$3" try ts el
  for try in 1 2; do
    islo rm "$box" >>"$LOG" 2>&1 || true
    ts=$(now)
    if [ "$mode" = cold ]; then
      inject "$box"
      islo use "$box" --cpu 2 --memory "$FORKMEM" --disk 25 -- bash -lc "\$HOME/gx/warmup.sh $chr >/dev/null 2>&1 && python3 \$HOME/gx/compute.py $chr" > "$DATA/_wg_${mode}_$chr.out" 2>>"$LOG"
    else
      islo use "$box" --snapshot "$SNAP" --cpu 2 --memory "$FORKMEM" -- bash -lc "python3 \$HOME/gx/compute.py $chr" > "$DATA/_wg_${mode}_$chr.out" 2>>"$LOG"
    fi
    el=$(python3 -c "print(round($(now)-$ts,2))")
    if grep -q '^{"chrom"' "$DATA/_wg_${mode}_$chr.out"; then
      grep '^{"chrom"' "$DATA/_wg_${mode}_$chr.out" > "$DATA/wg_${mode}_$chr.json"
      echo "$el" > "$DATA/wg_${mode}_$chr.time"
      say "  [$mode] $chr ${el}s (try $try)"; islo rm "$box" >>"$LOG" 2>&1 || true; return 0
    fi
    say "  [$mode] $chr try $try: no JSON, retry"
  done
  say "  [$mode] $chr FAILED"; islo rm "$box" >>"$LOG" 2>&1 || true; return 1
}

cmd_run(){
  say "=== genome-wide reference-broadcast fan-out · ${#CHROMS[@]} chromosomes · harness=Claude Code ==="
  command -v crabbox >/dev/null && say "crabbox present; coordinator: $(crabbox whoami 2>&1 | head -1)"

  # 1 · WARM BASE (download all 24 chromosomes + index, once)
  say "=== 1/5 warm base: toolchain + ${#CHROMS[@]} GRCh38 chromosomes + index ==="
  islo rm gx-wg-warm >>"$LOG" 2>&1 || true
  t=$(now); inject gx-wg-warm
  islo use gx-wg-warm --cpu 2 --memory 4096 --disk 25 -- bash -lc "\$HOME/gx/warmup.sh ${CHROMS[*]}" >>"$LOG" 2>&1
  WARM=$(python3 -c "print(round($(now)-$t,2))"); say "warm base built in ${WARM}s"

  # 2 · SNAPSHOT (the broadcast)
  say "=== 2/5 snapshot warm base -> $SNAP ==="
  t=$(now); islo snapshot save gx-wg-warm --name "$SNAP" >>"$LOG" 2>&1
  SNAPS=$(python3 -c "print(round($(now)-$t,2))")
  SNAPSZ=$(islo snapshot ls 2>/dev/null | grep -F "$SNAP" | awk '{print $2" "$3}')
  say "snapshot saved in ${SNAPS}s (size ${SNAPSZ:-n/a})"

  # 3 · COLD BASELINE (one large shard from scratch, for honest per-shard cold cost)
  say "=== 3/5 cold baseline: $COLD_REF from scratch ==="
  t=$(now); run_shard gx-wg-cold "$COLD_REF" cold; COLD=$(python3 -c "print(round($(now)-$t,2))")
  say "cold baseline ($COLD_REF) ${COLD}s"

  # 4 · WARM FAN-OUT (fork per chromosome, in waves of MAXP)
  say "=== 4/5 warm fan-out: fork ${#CHROMS[@]} shards from snapshot (MAXP=$MAXP) ==="
  t_map=$(now); i=0
  while [ $i -lt ${#CHROMS[@]} ]; do
    pids=()
    for j in $(seq 0 $((MAXP-1))); do
      idx=$((i+j)); [ $idx -ge ${#CHROMS[@]} ] && break
      run_shard "gx-wg-${CHROMS[$idx]}" "${CHROMS[$idx]}" warm & pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p"; done
    i=$((i+MAXP))
  done
  MAP=$(python3 -c "print(round($(now)-$t_map,2))"); say "warm fan-out wall-clock ${MAP}s"

  # 5 · REDUCE
  say "=== 5/5 reduce -> data/receipts.json (genome-wide) ==="
  python3 scripts/reduce_wg.py "$DATA" "$SNAP" "$SNAPSZ" "$WARM" "$SNAPS" "$COLD" "$COLD_REF" "$MAP" "$MAXP" "${CHROMS[@]}" | tee -a "$LOG"
  say "=== DONE (snapshot $SNAP retained) ==="
}

# ---- T2T follow-up: same harness, pointed at the complete T2T-CHM13 assembly ----
cmd_t2t(){
  SNAP="genomics-t2t-$TS"
  say "=== T2T-CHM13 vs GRCh38 delta · ${#CHROMS[@]} chromosomes · harness=Claude Code + OpenClaw Crabbox ==="
  # 1 · warm base: download + split the complete T2T assembly (paid once)
  say "=== 1/4 warm base: OpenClaw crabbox run --provider islo --keep ==="
  WARM_BOX="gx-t2t-warm"
  islo rm "$WARM_BOX" >>"$LOG" 2>&1 || true
  t=$(now)
  if have_crabbox_islo; then
    LEASE_JSON="$DATA/t2t_warm_lease.json"
    crabbox run --provider islo \
      --islo-image "$IMG" \
      --islo-vcpus 2 \
      --islo-memory-mb 4096 \
      --islo-disk-gb 40 \
      --keep \
      --lease-output "$LEASE_JSON" \
      --label t2t-warm \
      -- bash -lc "bash scripts/t2t_warmup.sh && mkdir -p \$HOME/gx && cp scripts/compute.py \$HOME/gx/compute.py" \
      >>"$LOG" 2>&1
    WARM_BOX="$(lease_sandbox_name "$LEASE_JSON")"
    say "warm base kept by Crabbox as $WARM_BOX"
  else
    say "ISLO_API_KEY/CRABBOX_ISLO_API_KEY missing; falling back to direct islo CLI for warm base"
    inject "$WARM_BOX"
    islo use "$WARM_BOX" --cpu 2 --memory 4096 --disk 40 -- bash -lc "\$HOME/gx/t2t_warmup.sh" >>"$LOG" 2>&1
  fi
  WARM=$(python3 -c "print(round($(now)-$t,2))"); say "T2T warm base built in ${WARM}s"
  # 2 · snapshot
  say "=== 2/4 snapshot -> $SNAP ==="
  t=$(now); islo snapshot save "$WARM_BOX" --name "$SNAP" >>"$LOG" 2>&1
  SNAPS=$(python3 -c "print(round($(now)-$t,2))")
  SNAPSZ=$(islo snapshot ls 2>/dev/null | grep -F "$SNAP" | awk '{print $2" "$3}')
  say "snapshot saved in ${SNAPS}s (size ${SNAPSZ:-n/a})"
  if have_crabbox_islo; then
    crabbox stop --provider islo "$WARM_BOX" >>"$LOG" 2>&1 || true
  else
    islo rm "$WARM_BOX" >>"$LOG" 2>&1 || true
  fi
  # 3 · fan-out (waves of MAXP)
  say "=== 3/4 T2T fan-out: fork ${#CHROMS[@]} shards (MAXP=$MAXP) ==="
  t_map=$(now); i=0
  while [ $i -lt ${#CHROMS[@]} ]; do
    pids=()
    for j in $(seq 0 $((MAXP-1))); do
      idx=$((i+j)); [ $idx -ge ${#CHROMS[@]} ] && break
      run_t2t_shard "gx-t2t-${CHROMS[$idx]}" "${CHROMS[$idx]}" & pids+=($!)
    done
    for p in "${pids[@]}"; do wait "$p"; done
    i=$((i+MAXP))
  done
  MAP=$(python3 -c "print(round($(now)-$t_map,2))"); say "T2T fan-out wall-clock ${MAP}s"
  # 4 · reduce vs GRCh38 (reuses data/wg_warm_*.json from the genome-wide run)
  say "=== 4/4 reduce: T2T vs GRCh38 -> data/compare.json ==="
  ensure_grch38_receipts
  python3 scripts/reduce_compare.py "$DATA" "$SNAP" "$SNAPSZ" "$WARM" "$SNAPS" "$MAP" "${CHROMS[@]}" | tee -a "$LOG"
  say "=== DONE (snapshot $SNAP retained) ==="
}

case "${1:-run}" in
  run) cmd_run ;;
  t2t) cmd_t2t ;;
  *) echo "usage: crabbox.sh {run|t2t}"; exit 1 ;;
esac
