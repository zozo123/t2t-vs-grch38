#!/usr/bin/env bash
# T2T warm-up: download the complete T2T-CHM13v2.0 assembly (UCSC hs1, ~975 MB),
# split it into per-chromosome FASTA, build an index. Paid once, then snapshotted.
set -euo pipefail
GX="$HOME/gx"; REF="$GX/ref_t2t"
mkdir -p "$REF"; cd "$GX"

echo "[t2t] installing toolchain..."
pip install --quiet --no-input numpy >/dev/null 2>&1

if [ ! -f "$REF/.split_done" ]; then
  echo "[t2t] downloading hs1.fa.gz (~975 MB)..."
  GX_DST="$REF/hs1.fa.gz" python3 - <<'PY'
import os, urllib.request
urllib.request.urlretrieve("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz", os.environ["GX_DST"])
print("  downloaded %.0f MB" % (os.path.getsize(os.environ["GX_DST"])/1e6))
PY
  echo "[t2t] splitting into per-chromosome FASTA (streaming)..."
  GX_REF="$REF" python3 - <<'PY'
import gzip, os
ref = os.environ["GX_REF"]
src = os.path.join(ref, "hs1.fa.gz")
keep = {f"chr{i}" for i in range(1,23)} | {"chrX","chrY"}   # 24 chromosomes; skip chrM
out = None; name = None; wrote = {}
with gzip.open(src, "rt") as f:
    for line in f:
        if line.startswith(">"):
            if out: out.close()
            name = line[1:].strip().split()[0]
            if name in keep:
                out = open(os.path.join(ref, f"{name}.fa"), "w"); out.write(line); wrote[name]=0
            else:
                out = None
        elif out is not None:
            out.write(line); wrote[name]+=len(line.strip())
if out: out.close()
os.remove(src)  # reclaim ~975 MB
open(os.path.join(ref, ".split_done"), "w").close()
print("  split:", ", ".join(f"{k}={v}" for k,v in sorted(wrote.items())))
PY
fi
echo "[t2t] warm base ready at $REF"
