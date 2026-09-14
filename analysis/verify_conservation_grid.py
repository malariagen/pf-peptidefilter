"""Sweep the haplotype gene rule to find an on-grid conservation-stringent scenario.

Supports the O5 decision: strain_conservation=0.995 is off the app 0.01 slider grid,
so the stringent arm tightens haplotype_gene_pc instead, at threshold 0.99.
Run with the pinned environment (Python 3.10 / pandas 1.5.1).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import pepfilter_headless as H  # noqa: E402

OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(OUT, "conservation_grid.txt")
_lines = []


def say(s=""):
    print(s)
    _lines.append(s)

genes, peptides = H.load_data()
H.validate(genes, peptides)

ID = "gene_id" if "gene_id" in genes.columns else genes.columns[0]

def gset(cfg):
    res = H.run(cfg, genes, peptides)
    return set(res[2][ID])

base = H.m253()
m253 = gset(base)
say(f"\nM253 baseline: {len(m253)} genes  (haplotype_gene_pc=95, strain_conservation=0.99)\n")

say("OPTION 2 — tighten the gene rule on-grid, threshold held at 0.99")
say(f"{'gene_pc':>8} {'column':>8} {'genes':>7} {'vs M253':>9} {'subset':>7}")
for pc in [95, 96, 97, 98, 99, 100]:
    cfg = dict(base); cfg["haplotype_gene_pc"] = pc
    g = gset(cfg)
    say(f"{pc:>8} {'_p%02d'%(100-pc):>8} {len(g):>7} {len(g)-len(m253):>+9} {str(g<=m253):>7}")

say("\nFor comparison — the current off-grid arm and its neighbours")
for thr in [0.99, 0.995, 0.999, 1.0]:
    cfg = dict(base); cfg["strain_conservation"] = thr
    g = gset(cfg)
    grid = "on-grid" if abs(thr*100 - round(thr*100)) < 1e-9 else "OFF-GRID"
    say(f"  strain_conservation={thr:<6} {len(g):>4} genes  {grid}")

with open(LOG, "w") as fh:
    fh.write("\n".join(_lines) + "\n")
print(f"\nwrote {LOG}")
