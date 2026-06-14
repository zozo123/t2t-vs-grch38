# What T2T changed in the human reference

**A genome-wide GRCh38 -> T2T-CHM13 delta, computed by an agent-run sandbox harness on [islo.dev](https://islo.dev).**

Live site: https://zozo123.github.io/t2t-vs-grch38/

This is the follow-up to the reference-broadcast demo:
https://zozo123.github.io/genomics-sandboxes/

## Result

T2T-CHM13 closes many gaps in GRCh38. Running the same CpG-island kernel across both references shows:

| metric | value |
|---|---:|
| Newly usable sequence in T2T | 179.6 Mb |
| GRCh38 gap sequence removed | 150.6 Mb |
| GRCh38 candidate CpG islands | 264,816 |
| T2T candidate CpG islands | 305,308 |
| Extra candidate islands | 40,492 |
| Candidate-island density in added sequence | 225.4/Mb |

The added sequence is not one signal:

- chrY and chr9 add large amounts of sequence with little candidate-island gain, consistent with newly resolved heterochromatin and repeat-rich sequence.
- chr13, chr15, chr21, and chr22 gain many candidate islands, consistent with acrocentric/rDNA-rich sequence and with the known tendency of simple CpG-island rules to over-call GC/CpG-rich repeats.

This is not framed as a novel biological discovery. It is a compact, reproducible measurement of how completing the reference changes what a transparent sequence-composition caller sees.

## Harness

The harness is a Claude Code agent running:

```bash
./crabbox.sh t2t
```

The harness:

1. Downloads UCSC `hs1.fa.gz` (T2T-CHM13v2.0) and splits it by chromosome.
2. Warms an islo sandbox with the reference and kernel.
3. Saves a VM snapshot.
4. Restores/forks one shard per chromosome.
5. Reduces T2T per-chromosome JSON against the GRCh38 genome-wide receipts.

Measured run:

| phase | time |
|---|---:|
| T2T warm-up | 77.6 s |
| snapshot save | 15.0 s |
| 24-way fan-out | 60.0 s |
| snapshot size | 1009.5 MB |

## Files

| file | purpose |
|---|---|
| `index.html`, `styles.css`, `script.js` | static GitHub Pages paper |
| `data/compare.json` | comparison receipts used by the figures |
| `crabbox.sh` | genome-wide harness with `run` and `t2t` modes |
| `scripts/t2t_warmup.sh` | T2T warm-up: download, split, index |
| `scripts/compute.py` | per-chromosome MAP kernel |
| `scripts/reduce_compare.py` | GRCh38 vs T2T reducer |

New islo.dev accounts can use coupon `YOSSI150` for 150 free credits.

## Caveats

The CpG-island calls are candidates from a simple sequence-composition rule. The 1987 Gardiner-Garden and Frommer criterion is transparent and useful for this demonstration, but it is known to over-call in some GC-rich repeat contexts. The result should be read as a reference-composition comparison, not as promoter annotation, methylation measurement, or medical inference.
