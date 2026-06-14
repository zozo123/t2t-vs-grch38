#!/usr/bin/env python3
"""MAP kernel: per-chromosome genomic landscape.

Given a chromosome already present in the warm base (/opt/gx/ref/<chr>.fa.gz),
compute the signals ENCODE-style regulatory genomics cares about:

  - GC% in 1 Mb bins  (the isochore / GC landscape)
  - CpG observed/expected ratio per bin
  - CpG-island candidates (Gardiner-Garden & Frommer 1987:
        200 bp window, GC > 50%, obs/exp CpG > 0.6)
  - N-gap fraction per bin (assembly gaps / centromere)

Pure-stdlib + numpy so it runs identically in a cold box and a warm fork.
Emits one JSON object on stdout. This is the unit of work a single sandbox does.
"""
import sys, time, json, gzip, os
import numpy as np

BIN = 1_000_000          # 1 Mb landscape bins
ISL_WIN = 200            # CpG-island detection window (bp)
ISL_STEP = 200           # non-overlapping for speed; honest, slightly conservative

def load_seq(path):
    opener = gzip.open if path.endswith(".gz") else open
    parts = []
    with opener(path, "rt") as f:
        for line in f:
            if line and line[0] != ">":
                parts.append(line.strip())
    return "".join(parts).upper()

def main():
    chrom = sys.argv[1]
    ref = os.environ.get("GX_REF", os.path.expanduser("~/gx/ref"))
    path = os.path.join(ref, f"{chrom}.fa.gz")
    if not os.path.exists(path):
        alt = os.path.join(ref, f"{chrom}.fa")
        if os.path.exists(alt):
            path = alt
    t0 = time.time()
    s = load_seq(path)
    n = len(s)

    # Encode to bytes once; vectorize with numpy. Keep only the masks we reuse
    # (C, G, N, CpG) as bool to stay memory-frugal on the largest chromosomes.
    b = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    del s
    counts = np.bincount(b, minlength=256)
    C = (b == 67); G = (b == 71); N = (b == 78)
    cg = C[:-1] & G[1:]   # CpG: base==C followed by base==G

    a_tot = int(counts[65]); t_tot = int(counts[84])
    c_tot = int(counts[67]); g_tot = int(counts[71])
    gc_tot = c_tot + g_tot
    at_tot = a_tot + t_tot
    n_tot = int(counts[78])
    cpg_tot = int(cg.sum())
    # genome-wide obs/exp CpG = (CpG/L) / ((C/L)*(G/L)) = CpG*L / (C*G)
    usable = gc_tot + at_tot
    oe = (cpg_tot * usable / (c_tot * g_tot)) if (c_tot and g_tot) else 0.0

    # 1 Mb bins (counts via bincount on each slice — cheap, no full A/T arrays)
    bins = []
    for start in range(0, n, BIN):
        end = min(start + BIN, n)
        cnt = np.bincount(b[start:end], minlength=256)
        ci = int(cnt[67]); gi = int(cnt[71]); nb = int(cnt[78])
        atb = int(cnt[65] + cnt[84])
        gcb = ci + gi
        cpgi = int(cg[start:min(end, n - 1)].sum())
        denom = gcb + atb
        gc_pct = (100.0 * gcb / denom) if denom else 0.0
        oeb = (cpgi * denom / (ci * gi)) if (ci and gi) else 0.0
        bins.append({
            "start": start,
            "gc": round(gc_pct, 2),
            "cpg_oe": round(oeb, 3),
            "n_frac": round(nb / (end - start), 3),
        })

    # CpG-island candidates (Gardiner-Garden & Frommer), non-overlapping windows.
    # reduceat directly on bool masks -> small int outputs (no big int32 copies).
    islands = 0
    nwin = n // ISL_STEP
    if nwin:
        idx = np.arange(0, nwin * ISL_STEP, ISL_STEP)
        idxcg = np.arange(0, min(nwin * ISL_STEP, n - 1), ISL_STEP)
        cC = np.add.reduceat(C, idx)
        cG = np.add.reduceat(G, idx)
        cN = np.add.reduceat(N, idx)
        cCG = np.add.reduceat(cg, idxcg)
        m = min(len(cC), len(cCG))
        cC, cG, cN, cCG = cC[:m], cG[:m], cN[:m], cCG[:m]
        win_len = ISL_WIN - cN  # usable bases
        with np.errstate(divide="ignore", invalid="ignore"):
            gcw = np.where(win_len > 0, (cC + cG) / win_len, 0)
            oew = np.where((cC > 0) & (cG > 0), cCG * win_len / (cC * cG), 0)
        islands = int(np.sum((win_len >= 150) & (gcw > 0.5) & (oew > 0.6)))

    out = {
        "chrom": chrom,
        "length": n,
        "n_bases": n_tot,
        "gc_pct": round(100.0 * gc_tot / usable, 2) if usable else 0.0,
        "cpg_count": cpg_tot,
        "cpg_oe": round(oe, 3),
        "cpg_islands": islands,
        "bins": bins,
        "compute_sec": round(time.time() - t0, 2),
    }
    print(json.dumps(out))

if __name__ == "__main__":
    main()
