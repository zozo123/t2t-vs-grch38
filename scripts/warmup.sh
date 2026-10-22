#!/usr/bin/env bash
# GRCh38 warm-up, run once in the base sandbox before it is snapshotted:
#   1. install numpy (the kernel's only dependency)
#   2. download the requested UCSC hg38 chromosomes and check their md5
#   3. record chromosome lengths in ref/index.json
#   usage: warmup.sh chr1 chr2 ...
set -euo pipefail
GX="$HOME/gx"
mkdir -p "$GX/ref"
cd "$GX"

echo "[warmup] installing numpy..."
pip install --quiet --no-input numpy==2.4.6 >/dev/null

echo "[warmup] downloading reference chromosomes..."
GX_REF="$GX/ref" python3 - "$@" <<'PY'
import hashlib, os, sys, urllib.request
BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
# From $BASE/md5sum.txt (files dated 2014-01-24).
MD5 = {
    "chr1": "f069c41e7cc8c2d3a7655cbb2d4186b8", "chr2": "609cc41e8a44eb99f87839456c6ed333",
    "chr3": "410294b87e71dbd40cf1b98fbadb443b", "chr4": "2a54872d6f884b8a993a7aa16db288df",
    "chr5": "69ee0b7e266a45d1b21fd665f917b009", "chr6": "fe629afc52a1a06468b60fde0109e32a",
    "chr7": "9c1df22d2a5ce712a410c65fb6e4aeb6", "chr8": "44c1556e73628c8ccbaa62ec39bae86a",
    "chr9": "f7e98217b35f5a7451f1166196b4b33c", "chr10": "b7d478cb06a84400f74a8932cd5fad0d",
    "chr11": "76542ddafa00440fb67606649e387d7a", "chr12": "ea464740f7571bd35d8bc4776be874fa",
    "chr13": "faddb3ad7e9c468ade70b1c88171f1eb", "chr14": "c970b8826b6434daa4292d95a300bfd9",
    "chr15": "fd152c2f9d14aaa1b9267e8a1d7561ff", "chr16": "6add68e486a69be4ae8cbbc9453cc173",
    "chr17": "023ccefd298db339159a446e88d742ac", "chr18": "f3c68bdafea015e6f724c3d646fa2282",
    "chr19": "5668c1a7edca32cc5d68bbeb26061d81", "chr20": "1d42fbb5014b165373aca053f4f24061",
    "chr21": "184df2bd9b812b6e6b6da16c6021369e", "chr22": "41b47ce1cc21b558409c19b892e1c0d1",
    "chrX": "58259f5bd35fe7e16aeec4c24ee02206", "chrY": "d26954bb71ff5e3b8f794fa4c61bfb15",
}
def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
for chrom in sys.argv[1:]:
    dst = os.path.join(os.environ["GX_REF"], f"{chrom}.fa.gz")
    if not os.path.exists(dst):
        urllib.request.urlretrieve(f"{BASE}/{chrom}.fa.gz", dst)
    got = md5(dst)
    if got != MD5[chrom]:
        os.remove(dst)
        sys.exit(f"[warmup] md5 mismatch for {chrom}: got {got}, expected {MD5[chrom]}")
    print(f"  {chrom}: {os.path.getsize(dst)/1e6:.1f} MB, md5 ok")
PY

echo "[warmup] recording chromosome lengths..."
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
print("  lengths:", ", ".join(f"{k}={v}" for k, v in idx.items()))
PY

echo "[warmup] done. warm base ready at $GX"
