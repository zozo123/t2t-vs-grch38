# What T2T changed in the human reference

Per-chromosome comparison of GRCh38 and T2T-CHM13v2.0: non-N sequence and counts of 200-bp windows passing the Gardiner-Garden & Frommer CpG rule.

**Result:** +179.6 Mb usable sequence (net), +40,492 windows (+15%). Heterochromatic additions (chrY, chr9, chr16) add few windows; the acrocentric chromosomes (chr13–15, 21, 22) add 80% of them.

**Caveat:** the windows are not annotated CpG islands; this is sequence composition only.

[Reproduce](#reproduce) · [Live page](https://zozo123.github.io/t2t-vs-grch38/) · run 2026-06-14 (commit `5357a12`)

## What is measured

For chr1–22, X and Y of UCSC `hg38` (GRCh38) and UCSC `hs1` (T2T-CHM13v2.0), `scripts/compute.py` reports:

- **usable sequence**: non-N bases.
- **CpG-rule windows** (JSON key `islands`, kept for compatibility): non-overlapping 200-bp windows passing the Gardiner-Garden & Frommer rule (GC > 50%, observed/expected CpG > 0.6, at least 150 non-N bases) [3], tiled from the chromosome start. They are not merged into islands and are not annotated CpG islands.

## Result

| metric | value |
|---|---:|
| Usable sequence, T2T − GRCh38 (net) | +179.6 Mb |
| GRCh38 N (gap) sequence | 150.6 Mb |
| CpG-rule windows, GRCh38 | 264,816 |
| CpG-rule windows, T2T | 305,308 |
| Δ windows | +40,492 (+15%) |
| Net Δwindows / net ΔMb | 225.4 (ratio of chromosome-level deltas, not a per-region density) |

Notable chromosomes (T2T − GRCh38; density is windows per Mb of usable sequence):

| chrom | Δ usable | Δ windows | density, GRCh38 → T2T |
|---|---:|---:|---:|
| chrY | +36.0 Mb | +263 | 79.7 → 37.9 |
| chr9 | +28.8 Mb | +324 | 94.7 → 78.7 |
| chr16 | +14.5 Mb | +232 | 147.8 → 128.0 |
| chr13 | +15.6 Mb | +10,487 | 68.5 → 151.4 |
| chr14 | +10.6 Mb | +4,422 | 85.1 → 119.9 |
| chr15 | +15.1 Mb | +6,777 | 94.1 → 147.8 |
| chr21 | +5.0 Mb | +6,925 | 111.9 → 253.0 |
| chr22 | +12.2 Mb | +3,728 | 174.0 → 205.4 |
| chrX | −0.6 Mb | +235 | 66.6 → 68.4 |

- chrY, chr9 and chr16 add large blocks of sequence but few windows, so their density falls. Most of what T2T adds on these chromosomes is satellite heterochromatin [1, 2].
- The five acrocentric chromosomes add 32,339 of the 40,492 extra windows. GC- and CpG-rich repeats can pass a composition rule without being promoter-type CpG islands [4]; which of these windows are repeats was not measured.
- chrX is the only chromosome with less non-N sequence in T2T.

All per-chromosome values are in [`data/compare.json`](data/compare.json); its `definitions` block documents each field.

## Limits

- **Windows are not CpG islands.** On chr21 the kernel gives 4,484 windows on GRCh38 (2,781 runs after merging adjacent windows) and 11,409 on T2T; UCSC's CpG-island track lists 446 (hg38) and 864 (hs1). Stricter, repeat-excluding criteria exist [4] and were not applied.
- **Composition only.** No repeat-annotation overlap was performed; no functional, promoter or methylation claim is made.
- **179.6 Mb is net** (chrX is −0.6 Mb). GRCh38 centromeres are modeled sequence rather than N gaps, so its 150.6 Mb of N sequence undercounts what T2T resolved.
- **rDNA is modeled.** CHM13's rDNA arrays are model sequences [1]; window counts there scale with the modeled copy number.
- **chrY is from HG002**, not CHM13 [2].
- Reference assemblies only; nothing here is about any individual.

## Provenance

The published numbers were computed on 2026-06-14 by a Claude Code agent running `./crabbox.sh t2t` through the direct islo CLI path (commit `5357a12`): one [islo.dev](https://islo.dev) sandbox downloads and splits `hs1.fa.gz` and is saved as a snapshot; 24 forks of the snapshot each run `compute.py` on one chromosome. The GRCh38 values are the per-chromosome receipts of the earlier [genomics-sandboxes](https://github.com/zozo123/genomics-sandboxes) run of the same kernel ([live page](https://zozo123.github.io/genomics-sandboxes/)).

An [OpenClaw Crabbox](https://github.com/openclaw/crabbox) shard path was added to `crabbox.sh` later (`61d1488`) and has not been re-run.

| phase | time |
|---|---:|
| T2T warm-up | 77.6 s |
| snapshot save | 15.0 s |
| 24 shards, 8 concurrent (MAXP=8) | 60.0 s |
| snapshot size | 1009.5 MB |

## Reproduce

### Locally, one chromosome (no sandbox)

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install numpy==2.4.6
mkdir -p ref
curl -fsSL -o ref/chr21.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz
GX_REF=ref python3 scripts/compute.py chr21 | python3 -c 'import json,sys; print(json.load(sys.stdin)["cpg_islands"])'
# 4484
```

The full JSON output matches the published GRCh38 chr21 receipt except `compute_sec`. The file's md5 is `184df2bd9b812b6e6b6da16c6021369e` (UCSC `md5sum.txt`).

### All 24 chromosomes (islo.dev sandboxes)

```bash
./crabbox.sh t2t      # needs the islo CLI, logged in
diff <(python3 -m json.tool data/compare.json) <(python3 -m json.tool data/compare.new.json)
```

- Writes `data/compare.new.json`, never `data/compare.json`, and exits non-zero if any of the 24 shard outputs is missing.
- Shards run through `crabbox run --provider islo` when `crabbox` is installed and `ISLO_API_KEY` is set (not yet exercised), otherwise through the islo CLI (the path that produced the published numbers).
- GRCh38 receipts are fetched from genomics-sandboxes at a pinned commit unless `./crabbox.sh grch38` has recomputed all 24 locally.
- Pinned: image `python:3.12.13-slim`, `numpy==2.4.6` (what the floating tags resolved to on the run date; the original run recorded no versions), crabbox flags checked against 0.36.0, and the UCSC md5 of every downloaded FASTA.

## Files

| file | purpose |
|---|---|
| `index.html`, `styles.css`, `script.js` | static GitHub Pages site, rendered from `data/compare.json` |
| `data/compare.json` | published per-chromosome results |
| `crabbox.sh` | islo.dev sandbox harness: `t2t` and `grch38` modes |
| `scripts/compute.py` | per-chromosome kernel |
| `scripts/reduce_compare.py` | GRCh38 vs T2T reducer, writes `data/compare.new.json` |
| `scripts/t2t_warmup.sh` | T2T warm-up: download `hs1.fa.gz`, md5-check, split per chromosome |
| `scripts/warmup.sh` | GRCh38 warm-up: download hg38 chromosomes, md5-check |
| `og-card.html`, `og.png` | social card template and its 1200×630 render |

## References

1. Nurk S. et al. The complete sequence of a human genome. *Science* 376:44–53 (2022). [doi:10.1126/science.abj6987](https://doi.org/10.1126/science.abj6987)
2. Rhie A. et al. The complete sequence of a human Y chromosome. *Nature* 621:344–354 (2023). [doi:10.1038/s41586-023-06457-y](https://doi.org/10.1038/s41586-023-06457-y)
3. Gardiner-Garden M., Frommer M. CpG islands in vertebrate genomes. *J Mol Biol* 196:261–282 (1987). [doi:10.1016/0022-2836(87)90689-9](https://doi.org/10.1016/0022-2836(87)90689-9)
4. Takai D., Jones P.A. Comprehensive analysis of CpG islands in human chromosomes 21 and 22. *PNAS* 99:3740–3745 (2002). [doi:10.1073/pnas.052410099](https://doi.org/10.1073/pnas.052410099)

## License

Code: MIT, see [LICENSE](LICENSE). Reference data: UCSC Genome Browser (hg38) and the T2T Consortium (CHM13v2.0); see their terms.
