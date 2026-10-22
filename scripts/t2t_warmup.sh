#!/usr/bin/env bash
# T2T warm-up, run once in the base sandbox before it is snapshotted: download the
# complete T2T-CHM13v2.0 assembly (UCSC hs1.fa.gz, ~975 MB), check its md5, and split
# it into per-chromosome FASTA (chr1-22, X, Y; chrM skipped) under ~/gx/ref_t2t.
set -euo pipefail
GX="$HOME/gx"; REF="$GX/ref_t2t"
mkdir -p "$REF"; cd "$GX"

echo "[t2t] installing numpy..."
pip install --quiet --no-input numpy==2.4.6 >/dev/null

if [ ! -f "$REF/.split_done" ]; then
  echo "[t2t] downloading hs1.fa.gz (~975 MB)..."
  GX_DST="$REF/hs1.fa.gz" python3 - <<'PY'
import hashlib, os, sys, urllib.request
dst = os.environ["GX_DST"]
# From https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/md5sum.txt (file dated 2022-07-18).
MD5 = "a493d5402cc86ecc3f54f6346d980036"
if not os.path.exists(dst):
    urllib.request.urlretrieve("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz", dst)
h = hashlib.md5()
with open(dst, "rb") as f:
    for block in iter(lambda: f.read(1 << 20), b""):
        h.update(block)
if h.hexdigest() != MD5:
    os.remove(dst)
    sys.exit(f"[t2t] md5 mismatch for hs1.fa.gz: got {h.hexdigest()}, expected {MD5}")
print("  downloaded %.0f MB, md5 ok" % (os.path.getsize(dst) / 1e6))
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
