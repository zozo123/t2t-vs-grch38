#!/usr/bin/env bash
# WARM-UP: the expensive setup that we pay ONCE and bake into a snapshot.
#   1. install the genomics toolchain
#   2. download the reference chromosomes (the data broadcast)
#   3. stage a faidx-style index
# A cold map shard has to run all of this before it can compute.
# A warm fork inherits all of it for free.
set -euo pipefail
GX="$HOME/gx"
mkdir -p "$GX/ref"
cd "$GX"

echo "[warmup] installing toolchain..."
pip install --quiet --no-input numpy biopython pyfaidx >/dev/null 2>&1

echo "[warmup] downloading reference chromosomes..."
BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
for chr in "$@"; do
  dst="$GX/ref/${chr}.fa.gz"
  if [ ! -f "$dst" ]; then
    GX_DST="$dst" GX_URL="$BASE/${chr}.fa.gz" GX_CHR="$chr" python3 - <<'PY'
import os, urllib.request
dst, url, chr = os.environ["GX_DST"], os.environ["GX_URL"], os.environ["GX_CHR"]
urllib.request.urlretrieve(url, dst)
print(f"  {chr}: {os.path.getsize(dst)/1e6:.1f} MB")
PY
  fi
done

echo "[warmup] building faidx index..."
GX_REF="$GX/ref" python3 - <<'PY'
import gzip, glob, os, json
ref = os.environ["GX_REF"]
idx = {}
for p in sorted(glob.glob(os.path.join(ref, "*.fa.gz"))):
    name = os.path.basename(p).replace(".fa.gz", "")
    n = 0
    with gzip.open(p, "rt") as f:
        for line in f:
            if line[0] != ">":
                n += len(line.strip())
    idx[name] = n
json.dump(idx, open(os.path.join(ref, "index.json"), "w"))
print("  indexed:", ", ".join(f"{k}={v}" for k, v in idx.items()))
PY

echo "[warmup] done. warm base ready at $GX"
