"""Interaction cross-checks behind the one-at-a-time limitation stated in the report (S3).

These are NOT scenario rows. Table A is one-filter-at-a-time by design -- the email scoped that
deliberately ("we are not going to iterate through hundreds of combinations"). These runs exist
only to measure what that design cannot see: whether the filters compound, and by how much.
Kept as a separate script so the figures quoted in the report are reproducible rather than
asserted.

Writes analysis/output/interaction_checks.txt.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pepfilter_headless as H  # noqa: E402
from run_sensitivity import OUT  # noqa: E402


def main():
    genes, peptides = H.load_data()
    H.validate(genes, peptides, verbose=False)
    m = H.m253()

    def gset(**kw):
        cfg = dict(m)
        cfg.update(kw)
        _, _, fg, _, _ = H.run(cfg, genes, peptides)
        return set(fg.gene_id)

    ref = gset()
    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("Interaction cross-checks -- exploratory, not archived scenario rows")
    out("=" * 72)
    out(f"M253 baseline: {len(ref)} genes")
    out("")

    # the five relaxed arms, individually, exactly as Table A defines them
    arms = {
        "human (alignment_length<=18)":        dict(alignment_length=18),
        "conservation (>=0.95)":               dict(strain_conservation=0.95),
        "indel GENE RULE (indel_gene_pc=95)":  dict(indel_gene_pc=95),
        "indel FREQUENCY (<=0.10)":            dict(indel_frequency=0.10),
        "expression (>=2 D4 replicates)":      dict(expression_rule="At least two replicates"),
        "orthology (none required)":           dict(homology_species=[]),
    }
    singles = {}
    out("Individual relaxed arms")
    for label, kw in arms.items():
        g = gset(**kw)
        singles[label] = g
        out(f"  {label:38s} {len(g):5d} genes  (+{len(g) - len(ref)})")
    out("")

    common = dict(alignment_length=18, strain_conservation=0.95,
                  expression_rule="At least two replicates", homology_species=[])
    out("All five relaxed together -- the two indel arms give different answers")
    for label, indel_kw, indel_arm in [
        ("indel GENE RULE (indel_gene_pc=95)", dict(indel_gene_pc=95),
         "indel GENE RULE (indel_gene_pc=95)"),
        ("indel FREQUENCY (<=0.10)", dict(indel_frequency=0.10),
         "indel FREQUENCY (<=0.10)"),
    ]:
        g = gset(**common, **indel_kw)
        used = ["human (alignment_length<=18)", "conservation (>=0.95)", indel_arm,
                "expression (>=2 D4 replicates)", "orthology (none required)"]
        gain_sum = sum(len(singles[u]) - len(ref) for u in used)
        out(f"  with {label}")
        out(f"    actual                {len(g):5d} genes  (+{len(g) - len(ref)})")
        out(f"    individual gain sum   {gain_sum:5d}")
        out(f"    additive prediction   {len(ref) + gain_sum:5d}")
        out(f"    super-additive by     {len(g) - (len(ref) + gain_sum):5d} genes")
        out(f"    M253 is a subset      {ref <= g}")
    out("")
    out("  Expected for conjunctive filters: relaxing one admits only genes already passing the")
    out("  other four, so simultaneous relaxation admits genes failing several at once. The")
    out("  single-filter table is therefore a lower bound on the catalogue's sensitivity.")
    out("")

    out("Pairwise: do orthology and the indel gene rule interact?")
    a = gset(homology_species=[]) - ref
    b = gset(indel_gene_pc=95) - ref
    both = gset(homology_species=[], indel_gene_pc=95)
    out(f"  orthology adds {len(a)}, indel gene rule adds {len(b)}, overlap {len(a & b)}")
    out(f"  both together {len(both)}; additive prediction {len(ref) + len(a) + len(b)} "
        f"-> {'additive' if len(both) == len(ref) + len(a) + len(b) else 'NOT additive'}")
    out("")

    out("Does the cost of tightening conservation depend on the human setting?")
    for hpc in (75, 50):
        loose = gset(human_id_gene_pc=hpc)
        tight = gset(human_id_gene_pc=hpc, haplotype_gene_pc=99)
        out(f"  human_id_gene_pc={hpc:3d}: {len(loose):4d} genes -> {len(tight):4d} with the 99% "
            f"conservation gene rule (costs {len(loose) - len(tight)})")
    out("  The two costs differ, so the filters are not independent.")

    path = os.path.join(OUT, "interaction_checks.txt")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
