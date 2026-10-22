#!/usr/bin/env python3
"""Reduce per-chromosome kernel outputs into the schema of data/compare.json.

Inputs (in DATA): wg_warm_<chr>.json (GRCh38) and t2t_<chr>.json (T2T), one per
chromosome, as printed by scripts/compute.py. Every requested chromosome must
have both files, or the reducer exits non-zero without writing anything.

Writes DATA/compare.new.json by default, so a rerun never overwrites the
published data/compare.json; pass --out to choose another path.

  reduce_compare.py DATA SNAPSHOT SNAPSHOT_SIZE WARMUP_SEC SNAPSHOT_SEC MAP_SEC CHROM... [--out PATH]
"""
import argparse, datetime, json, os, sys

# Shared with data/compare.json. `islands` keeps its historical key name.
DEFINITIONS = {
    "usable": "non-N bases",
    "n": "N bases (assembly gaps)",
    "gc": "GC % of usable bases",
    "islands": "CpG-rule windows: non-overlapping 200-bp windows (>= 150 non-N bases) with GC > 50% and CpG observed/expected > 0.6 (Gardiner-Garden & Frommer 1987 thresholds); not merged into islands, not annotated CpG islands",
    "density": "islands per Mb of usable sequence",
    "delta_bp": "t2t.usable - grch38.usable (net)",
    "delta_islands": "t2t.islands - grch38.islands",
    "new_seq_isl_per_mb": "extra_islands / (net_new_bp / 1e6): ratio of chromosome-level deltas, not a per-region density",
}


def keyf(c):
    s = c.replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, {"X": 23, "Y": 24}.get(s, 99))


def side(rec):
    usable = rec["length"] - rec["n_bases"]
    return {
        "usable": usable,
        "n": rec["n_bases"],
        "islands": rec["cpg_islands"],
        "gc": rec["gc_pct"],
        "density": round(rec["cpg_islands"] / (usable / 1e6), 1) if usable else 0,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("data")
    p.add_argument("snapshot")
    p.add_argument("snapshot_size")
    p.add_argument("warmup_sec", type=float)
    p.add_argument("snapshot_sec", type=float)
    p.add_argument("map_sec", type=float)
    p.add_argument("chroms", nargs="+")
    p.add_argument("--out", help="output path (default: DATA/compare.new.json)")
    p.add_argument("--harness", default="crabbox.sh t2t", help="how the shards were run")
    a = p.parse_args(argv)
    out_path = a.out or os.path.join(a.data, "compare.new.json")

    order = sorted(set(a.chroms), key=keyf)
    recs, missing = {}, []
    for c in order:
        for prefix in ("wg_warm", "t2t"):
            path = os.path.join(a.data, f"{prefix}_{c}.json")
            if not (os.path.exists(path) and os.path.getsize(path)):
                missing.append(os.path.basename(path))
                continue
            rec = json.load(open(path))
            if rec.get("chrom") != c:
                sys.exit(f"{path}: chrom is {rec.get('chrom')!r}, expected {c!r}")
            recs[prefix, c] = rec
    if missing:
        sys.exit(f"missing {len(missing)} shard output(s), nothing written: {', '.join(missing)}")

    rows = []
    for c in order:
        h, t = side(recs["wg_warm", c]), side(recs["t2t", c])
        rows.append({"chrom": c, "grch38": h, "t2t": t,
                     "delta_bp": t["usable"] - h["usable"],
                     "delta_islands": t["islands"] - h["islands"]})

    def total(k, f):
        return sum(r[k][f] for r in rows)

    net_new = total("t2t", "usable") - total("grch38", "usable")
    extra = total("t2t", "islands") - total("grch38", "islands")
    gainers = sorted(rows, key=lambda r: r["delta_bp"], reverse=True)[:6]

    out = {
        "assemblies": {"grch38": "GRCh38 (hg38)", "t2t": "T2T-CHM13v2.0 (hs1)"},
        "harness": a.harness,
        "run_date": datetime.datetime.now(datetime.timezone.utc).date().isoformat(),
        "definitions": DEFINITIONS,
        "snapshot": {"name": a.snapshot,
                     "size": a.snapshot_size.replace("ready", "").strip() or "n/a",
                     "save_sec": a.snapshot_sec},
        "timings": {"warmup_sec": a.warmup_sec, "snapshot_save_sec": a.snapshot_sec,
                    "map_wallclock_sec": a.map_sec, "shards": len(rows)},
        "totals": {
            "grch38_usable_bp": total("grch38", "usable"), "grch38_n_bp": total("grch38", "n"),
            "t2t_usable_bp": total("t2t", "usable"), "t2t_n_bp": total("t2t", "n"),
            "net_new_bp": net_new,
            "grch38_islands": total("grch38", "islands"), "t2t_islands": total("t2t", "islands"),
            "extra_islands": extra,
            "new_seq_isl_per_mb": round(extra / (net_new / 1e6), 1) if net_new else 0,
        },
        "top_gainers": [{"chrom": r["chrom"], "delta_bp": r["delta_bp"], "delta_mb": round(r["delta_bp"] / 1e6, 1)}
                        for r in gainers],
        "chromosomes": rows,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print(json.dumps({"wrote": out_path, "totals": out["totals"]}, indent=2))


if __name__ == "__main__":
    main()
