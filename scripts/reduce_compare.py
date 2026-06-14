#!/usr/bin/env python3
"""Compare T2T-CHM13 (hs1) vs GRCh38 (hg38) per chromosome → data/compare.json.

GRCh38 per-chrom results are reused from data/wg_warm_<chr>.json (the genome-wide run).
T2T per-chrom results are read from data/t2t_<chr>.json (this run)."""
import sys, json, os

data, snap, snapsz, warm, snaps, mapw = sys.argv[1:7]
chroms = sys.argv[6:]
warm = float(warm); snaps = float(snaps); mapw = float(mapw)

def load(prefix, chrom):
    p = os.path.join(data, f"{prefix}_{chrom}.json")
    return json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) else None

def keyf(c):
    s = c.replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, {"X": 23, "Y": 24}.get(s, 99))

order = sorted(chroms, key=keyf)
per = {}
tot = {"grch38": {"usable": 0, "n": 0, "islands": 0, "len": 0},
       "t2t": {"usable": 0, "n": 0, "islands": 0, "len": 0}}
for c in order:
    h = load("wg_warm", c)      # GRCh38
    t = load("t2t", c)          # T2T
    if not h or not t:
        continue
    h_us = h["length"] - h["n_bases"]
    t_us = t["length"] - t["n_bases"]
    per[c] = {
        "grch38": {"len": h["length"], "usable": h_us, "n": h["n_bases"],
                   "islands": h["cpg_islands"], "gc": h["gc_pct"],
                   "isl_per_mb": round(h["cpg_islands"]/(h_us/1e6), 1) if h_us else 0},
        "t2t": {"len": t["length"], "usable": t_us, "n": t["n_bases"],
                "islands": t["cpg_islands"], "gc": t["gc_pct"],
                "isl_per_mb": round(t["cpg_islands"]/(t_us/1e6), 1) if t_us else 0},
        "delta_bp": t_us - h_us,                 # newly resolved sequence (net)
        "delta_islands": t["cpg_islands"] - h["cpg_islands"],
    }
    for k, src in (("grch38", h), ("t2t", t)):
        tot[k]["usable"] += (src["length"] - src["n_bases"])
        tot[k]["n"] += src["n_bases"]; tot[k]["islands"] += src["cpg_islands"]; tot[k]["len"] += src["length"]

net_new = tot["t2t"]["usable"] - tot["grch38"]["usable"]
# island density of the *newly resolved* sequence (approx): extra islands / extra Mb
extra_isl = tot["t2t"]["islands"] - tot["grch38"]["islands"]
new_isl_per_mb = round(extra_isl / (net_new/1e6), 1) if net_new else 0
# biggest gainers
gainers = sorted(per.items(), key=lambda kv: kv[1]["delta_bp"], reverse=True)

out = {
    "assemblies": {"grch38": "GRCh38 (hg38)", "t2t": "T2T-CHM13v2.0 (hs1)"},
    "harness": "Claude Code agent",
    "snapshot": {"name": snap, "size": (snapsz.replace("ready", "").strip() if snapsz else "n/a"), "save_sec": snaps},
    "timings": {"warmup_sec": warm, "snapshot_save_sec": snaps, "map_wallclock_sec": mapw, "shards": len(per)},
    "totals": {
        "grch38_usable_bp": tot["grch38"]["usable"], "grch38_n_bp": tot["grch38"]["n"],
        "t2t_usable_bp": tot["t2t"]["usable"], "t2t_n_bp": tot["t2t"]["n"],
        "net_new_bp": net_new,
        "grch38_islands": tot["grch38"]["islands"], "t2t_islands": tot["t2t"]["islands"],
        "extra_islands": extra_isl, "new_seq_isl_per_mb": new_isl_per_mb,
    },
    "top_gainers": [{"chrom": c, "delta_bp": v["delta_bp"], "delta_mb": round(v["delta_bp"]/1e6, 1)} for c, v in gainers[:6]],
    "per_chrom": {c: per[c] for c in order if c in per},
}
json.dump(out, open(os.path.join(data, "compare.json"), "w"), indent=2)
print(json.dumps({k: out[k] for k in ("assemblies", "totals", "top_gainers", "timings")}, indent=2))
