"""
Every figure quoted in the report that is not already produced by
``run_sensitivity.py``, ``run_antigens.py`` or ``run_interaction_checks.py``.

Those three scripts cover the fifteen frozen scenarios, the antigen benchmark and the
one-at-a-time-versus-joint comparison.  A number of figures quoted in the report's
narrative sections came from ad-hoc runs of the same harness and had no shipped script,
so they could not be re-derived by a reader.  This script is that missing one: it emits
``output/report_figures.txt``, one labelled block per claim, so every narrative figure
can be checked against a re-run.

Gene counts come from the app's own ``filter_datasets`` via ``pepfilter_headless``.
Where a block compares the filter's verdict against a count of peptides -- the
percentile-versus-count block below -- the verdict is still the app's; only the
count it is compared against is computed here, which is the whole point of the
comparison.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pepfilter_headless as H          # noqa: E402
from run_sensitivity import ANTIGENS, OUT, IDENTITY_ONLY, build_scenarios  # noqa: E402

REPORT_FIGURES = os.path.join(OUT, "report_figures.txt")

# the gene-level indel column the filter tests at the default 100% gene rule
INDEL_COL = "max_indel_freq_type_any_fs_any_p100"

_lines = []


def say(s=""):
    print(s)
    _lines.append(s)


def head(n, title):
    say()
    say(f"{'=' * 78}")
    say(f"{n}. {title}")
    say(f"{'=' * 78}")


def n_genes(params, genes, peptides):
    _, _, fg, _, _ = H.run(params, genes, peptides)
    return fg


def counts(params, genes, peptides):
    fg = n_genes(params, genes, peptides)
    return len(fg), int(peptides["gene_id"].isin(fg.gene_id).sum())


# --------------------------------------------------------------------------------------
def block_percentile_vs_count(genes, peptides):
    """The conservation gene rule: interpolated percentile against a count of peptides.

    The filter's verdict is the app's (a conservation-isolated run).  The count it is
    compared against -- the true share of a gene's peptides at >=0.99 -- is computed
    here from the peptide table, because no app setting reports it.
    """
    head(1, "Percentile versus count, conservation filter (report section 04)")

    retained = set(n_genes(H.isolate("conservation"), genes, peptides).gene_id)
    m253 = set(n_genes(H.m253(), genes, peptides).gene_id)

    passing = peptides.groupby("gene_id")["hap_top_hap_freq"].apply(
        lambda s: float((s >= 0.99).sum()) / len(s))
    frac = passing.reindex(genes.gene_id.values).values
    p05 = genes["hap_top_hap_freq_p05"].values
    gid = genes.gene_id.values
    keep = np.array([g in retained for g in gid])

    say(f"  conservation isolated, gene rule 95 / threshold 0.99 : {len(retained)} genes")
    say()
    say("  Of the published 253:")
    in253 = np.array([g in m253 for g in gid])
    say(f"    kept on fewer than 95% of their peptides at >=0.99 : "
        f"{int(((frac < 0.95) & in253).sum())}")
    say()
    say("  Genome-wide, the two readings disagree as follows.")
    say("  Non-strict count reading (a gene 'meets the count' if fraction >= 0.95),")
    say("  which is the reading the report uses:")
    admitted_below = (keep & (frac < 0.95))
    rejected_above = (~keep & (frac >= 0.95))
    nan_p05 = np.isnan(p05)
    say(f"    admitted although fewer than 95% of peptides pass : {int(admitted_below.sum())}")
    say(f"    rejected although at least 95% of peptides pass   : {int(rejected_above.sum())}")
    say(f"      of which the percentile column is present (the interpolation disagrees "
        f"with the count) : {int((rejected_above & ~nan_p05).sum())}")
    for g, f, v in zip(gid[rejected_above & ~nan_p05], frac[rejected_above & ~nan_p05],
                       p05[rejected_above & ~nan_p05]):
        say(f"        {g}  true fraction {f:.4f}  hap_top_hap_freq_p05 {v:.5f}")
    say(f"      of which the percentile column is missing (NaN, so the test is False "
        f"whatever the peptides say) : {int((rejected_above & nan_p05).sum())}")
    for g, f in zip(gid[rejected_above & nan_p05], frac[rejected_above & nan_p05]):
        say(f"        {g}  true fraction {f:.4f}  hap_top_hap_freq_p05 NaN")
    say()
    say("  Strict count reading (fraction > 0.95), for completeness:")
    say(f"    admitted although 95% or fewer pass : {int((keep & (frac <= 0.95)).sum())}")
    say(f"    rejected although more than 95% pass: {int((~keep & (frac > 0.95)).sum())}")
    say()
    say("  Worked example quoted in the report:")
    ex = "PF3D7_0308800"
    pep_ex = peptides[peptides.gene_id == ex]["hap_top_hap_freq"]
    row = genes[genes.gene_id == ex].iloc[0]
    say(f"    {ex}: {len(pep_ex)} peptides, {int((pep_ex >= 0.99).sum())} at >=0.99 "
        f"({100.0 * (pep_ex >= 0.99).mean():.1f}%), hap_top_hap_freq_p05 "
        f"{row['hap_top_hap_freq_p05']:.4f}, in the 253: {ex in m253}")


def block_indel_arms(genes, peptides):
    head(2, "Indel filter: every control, to test for a stringent arm (departure 2)")
    m = H.m253()
    say(f"  baseline (any type, any effect, <=0.05, 100% rule) : "
        f"{counts(m, genes, peptides)[0]} genes")
    say()
    say("  frequency threshold swept down:")
    for t in [0.05, 0.04, 0.03, 0.02, 0.01, 0.0]:
        cfg = dict(m); cfg["indel_frequency"] = t
        say(f"    <={t:<5} : {counts(cfg, genes, peptides)[0]} genes")
    say(f"    lowest non-zero value present in the gene metric: "
        f"{genes[INDEL_COL].replace(0, np.nan).min():.4f}")
    say(f"  gene rule: indel_gene_pc default is {H.app_defaults()['indel_gene_pc']} "
        f"(maximum offered is {max(H.UI_CONFIG['indel_gene_pc']['options'])})")
    say()
    say("  which indels are counted (type x effect):")
    for typ in H.UI_CONFIG["indel_type"]["options"]:
        for fs in H.UI_CONFIG["indel_frameshfits"]["options"]:
            cfg = dict(m); cfg["indel_type"] = typ; cfg["indel_frameshfits"] = fs
            say(f"    type={typ:<16s} effect={fs:<16s} : "
                f"{counts(cfg, genes, peptides)[0]} genes")


def block_r21_rebuttals(genes, peptides):
    head(3, "R2.1: the two stage-independence runs (report section 03)")
    spo = dict(H.m253()); spo["expression_stage"] = ["Sporozoite Stage"]
    g, p = counts(spo, genes, peptides)
    fg = n_genes(spo, genes, peptides)
    say(f"  sporozoite stage substituted for day 4 : {g} genes / {p} peptides, "
        f"antigens retained: {sorted(k for k, v in ANTIGENS.items() if v in set(fg.gene_id)) or 'none'}")

    noexp = dict(H.m253()); noexp["expression_stage"] = []; noexp["homology_species"] = []
    g, p = counts(noexp, genes, peptides)
    fg = n_genes(noexp, genes, peptides)
    say(f"  expression and orthology both removed  : {g} genes / {p} peptides, "
        f"antigens retained: {sorted(k for k, v in ANTIGENS.items() if v in set(fg.gene_id)) or 'none'}")

    hi = dict(H.m253()); hi.update(H.NEUTRAL)
    ref = H.m253()
    for k in ["human_id_gene_pc", "identity_percent", "alignment_length",
              "indel_gene_pc", "indel_frequency"]:
        hi[k] = ref[k]
    fg = n_genes(hi, genes, peptides)
    got = {k: v for k, v in ANTIGENS.items() if v in set(fg.gene_id)}
    say(f"  human similarity and indels only       : {len(fg)} genes, "
        f"antigens retained: {sorted(got) or 'none'}")


def block_corrections(genes, peptides):
    head(4, "The six manuscript corrections (report section 04)")
    d = H.app_defaults()

    hc = dict(d); hc["indel_gene_pc"] = 0; hc["indel_frequency"] = 1.0
    say(f"  correction 1  conservation + human as the app ran it (human rule 75%) : "
        f"{counts(hc, genes, peptides)[0]} genes   [report: 1,274]")
    hc95 = dict(hc); hc95["human_id_gene_pc"] = 95
    say(f"                the same sentence run as written (human rule 95%)      : "
        f"{counts(hc95, genes, peptides)[0]} genes   [report: 158]")

    # The "simultaneously ... and" reading of the same sentence: not three percentile
    # tests ANDed, but one subset of peptides required to clear every criterion at once.
    # No app setting expresses that, so it is counted here over the app's own 1,274 set.
    set1274 = set(n_genes(hc, genes, peptides).gene_id)
    pep = peptides.assign(
        c=peptides["hap_top_hap_freq"] >= 0.99,
        i=peptides["blt_pident"] <= 80,
        l=peptides["blt_length"] <= 15)
    sub = pep[pep.gene_id.isin(set1274)]
    joint3 = sub.groupby("gene_id").apply(lambda t: (t.c & t.i & t.l).mean() >= 0.95)
    joint2 = sub.groupby("gene_id").apply(lambda t: (t.c & t.i).mean() >= 0.95)
    say(f"                joint-subset reading, all three criteria at once        : "
        f"{int(joint3.sum())} of {len(set1274)}   [report: 92]")
    say(f"                joint-subset reading, length criterion dropped         : "
        f"{int(joint2.sum())} of {len(set1274)}   [report: 1,133]")

    co = dict(d); co.update(H.NEUTRAL)
    co["haplotype_gene_pc"] = 95; co["strain_conservation"] = 0.99
    say(f"  correction 2  conservation only, 95% rule, >=0.99  : "
        f"{counts(co, genes, peptides)}   [report: 2,107 / 84,556]")
    co99 = dict(co); co99["haplotype_gene_pc"] = 99
    say(f"                conservation only, 99% rule, >=0.99  : "
        f"{counts(co99, genes, peptides)}   [report: 1,517 / 42,208]")
    co95 = dict(co); co95["strain_conservation"] = 0.95
    say(f"                conservation only, 95% rule, >=0.95  : "
        f"{counts(co95, genes, peptides)[0]} genes   [report: 3,824]")

    ind = H.isolate("indel")
    say(f"  correction 3  indel isolated at <=0.05 (isolated run, NaN leak included) : "
        f"{counts(ind, genes, peptides)[0]} genes")
    say(f"                genes at or below 0.05 on the indel metric alone           : "
        f"{int((genes[INDEL_COL] <= 0.05).sum())}   [report: 3,457]")
    say(f"                genes strictly below 0.05 on the indel metric alone        : "
        f"{int((genes[INDEL_COL] < 0.05).sum())}   "
        f"(the '<5%' and '<=5%' readings agree)")
    r95 = dict(H.m253()); r95["indel_gene_pc"] = 95
    say(f"                M253 with the indel gene rule at 95%                      : "
        f"{counts(r95, genes, peptides)}   [report: 291 / 13,772]")

    both = dict(H.m253()); both["homology_species"] = ["P. vivax", "P. berghei"]
    say(f"  correction 4  M253 with P. vivax only        : "
        f"{counts(H.m253(), genes, peptides)[0]} genes   [report: 253]")
    say(f"                M253 with P. vivax AND P. berghei : "
        f"{counts(both, genes, peptides)[0]} genes   [report: 251]")

    say()
    say("  correction 5  'jointly': how often the identity-passing and length-passing")
    say("                peptide subsets of a gene differ.")
    idp = peptides.assign(i=peptides["blt_pident"] <= 80, l=peptides["blt_length"] <= 15)
    byg = idp.groupby("gene_id")[["i", "l"]].sum()
    same_set = idp.groupby("gene_id").apply(lambda t: bool((t["i"] == t["l"]).all()))
    diff_frac = int((byg["i"] != byg["l"]).sum())
    diff_set = int((~same_set).sum())
    say(f"                genes whose passing counts differ            : {diff_frac}   [report: 4,840]")
    say(f"                genes whose passing subsets differ at all    : {diff_set}   [report: 4,846]")
    say(f"                equal counts but different peptides          : {diff_set - diff_frac}   [report: 6]")

    m253 = set(n_genes(H.m253(), genes, peptides).gene_id)
    both_h = idp[idp.gene_id.isin(m253)].groupby("gene_id").apply(
        lambda t: float((t["i"] & t["l"]).mean()))
    low = both_h[both_h < 0.75].sort_values()
    say(f"                of the 253, under 75% of peptides pass both human")
    say(f"                conditions at once                          : {len(low)}   [report: 10]")
    for g, f in low.items():
        say(f"                  {g}  {100.0 * f:.1f}%")


def block_conservation_boundary(genes):
    head(5, "How close the 253 sit to the 0.99 conservation cut (report section 05, S5)")
    m = set(pd.read_csv(os.path.join(OUT, "configs", "M253_genes.csv")).gene_id) \
        if os.path.exists(os.path.join(OUT, "configs", "M253_genes.csv")) else None
    if m is None:
        m = set(_M253_GENES)
    v = genes.set_index("gene_id").loc[sorted(m), "hap_top_hap_freq_p05"]
    for w, want in [(0.001, 18), (0.002, 33), (0.005, 89)]:
        say(f"  of the 253, within {w} of 0.99 : "
            f"{int(((v >= 0.99) & (v < 0.99 + w)).sum())}   [report: {want}]")
    say(f"  closest to the cut : {v[v >= 0.99].min():.5f}   [report: 0.99007]")
    allv = genes["hap_top_hap_freq_p05"]
    say(f"  genome-wide in [0.985, 0.995) : "
        f"{int(((allv >= 0.985) & (allv < 0.995)).sum())}   [report: 1,218]")


def block_expression(genes):
    head(6, "Day-2/4/5/6 expression, and the library-depth reading of it (S6)")
    say("  Genes at >=1 CPM, by stage. 'Median CPM' is the median over all replicate")
    say("  values pooled, not the median of per-gene medians; both are shown.")
    say(f"  {'stage':<12s} {'all reps':>9s} {'>=1 rep':>9s} {'median (pooled)':>17s} "
        f"{'median (of gene medians)':>26s} {'zero in all reps':>18s}")
    stages = [("Day 2", "exp_d2_cpm_r%d", 3), ("Day 4", "exp_d4_cpm_r%d", 3),
              ("Day 5", "exp_d5_cpm_r%d", 3), ("Day 6", "exp_d6_cpm_r%d", 3),
              ("Sporozoite", "exp_sporozoite_cpm_r%d", 5)]
    for label, tmpl, nrep in stages:
        cols = [tmpl % i for i in range(1, nrep + 1)]
        sub = genes[cols]
        allrep = int((sub >= 1).all(axis=1).sum())
        anyrep = int((sub >= 1).any(axis=1).sum())
        pooled = float(np.nanmedian(sub.values))
        of_med = float(np.nanmedian(sub.median(axis=1).values))
        zero = int((sub == 0).all(axis=1).sum())
        say(f"  {label:<12s} {allrep:>9d} {anyrep:>9d} {pooled:>17.2f} {of_med:>26.2f} "
            f"{zero:>18d}")
    say()
    say("  Summed CPM per replicate (a transcriptome-wide CPM total should be 1e6):")
    for label, tmpl, nrep in stages:
        tot = " / ".join(f"{genes[tmpl % i].sum():,.0f}" for i in range(1, nrep + 1))
        say(f"    {label:<12s} {tot}")


def block_defects(genes, peptides):
    head(7, "Defects found in shipped columns")
    true_n = peptides.groupby("gene_id").size().reindex(genes.gene_id.values).values
    delta = genes["num_peptides"].values - true_n
    say(f"  num_peptides: genes where it is exactly one short : "
        f"{int((delta == -1).sum())} of {len(genes)}   [report: 4,423 of 4,937]")
    say(f"               distinct differences observed        : {sorted(set(delta.tolist()))}")
    say(f"               column total vs peptide rows         : "
        f"{genes['num_peptides'].sum():,} vs {len(peptides):,}   "
        f"[report: 374,902 vs 379,325]")
    m253 = set(n_genes(H.m253(), genes, peptides).gene_id)
    sub = genes[genes.gene_id.isin(m253)]
    say(f"               over the published catalogue         : "
        f"{sub['num_peptides'].sum():,} vs "
        f"{int(peptides.gene_id.isin(m253).sum()):,}   [report: 9,861 vs 10,088]")
    say()
    pcols = [c for c in genes.columns if c.startswith("hap_top_hap_freq_p")]
    nan_all = genes[genes[pcols].isna().all(axis=1)]
    say(f"  hap_top_hap_freq percentile columns entirely NaN : {len(nan_all)} genes")
    say("  For each, what its own peptides actually say (this column is tested by the")
    say("  conservation filter, so a NaN is a silent rejection):")
    for _, r in nan_all.iterrows():
        pep = peptides[peptides.gene_id == r.gene_id]["hap_top_hap_freq"]
        have = pep.notna().sum()
        say(f"    {r.gene_id:<14s} {len(pep):>3d} peptides, {have:>3d} with a haplotype "
            f"frequency, min {pep.min() if have else float('nan'):.4f}, "
            f"share >=0.99 {100.0 * (pep >= 0.99).mean():.1f}%, in the 253: "
            f"{r.gene_id in m253}")


def block_tightening(genes, peptides):
    head(8, "Direction of change: relaxations and tightenings against the 253")
    m253 = set(n_genes(H.m253(), genes, peptides).gene_id)
    say(f"  {'scenario':<26s} {'genes':>6s} {'direction':>11s} "
        f"{'253 kept':>9s} {'strict subset':>14s} {'superset':>9s}")
    n_tight = n_relax = 0
    for label, filt, _, params in build_scenarios():
        if label in ("APP", "M253"):
            continue
        s = set(n_genes(params, genes, peptides).gene_id)
        direction = "tightening" if len(s) < 253 else ("relaxation" if len(s) > 253 else "equal")
        n_tight += direction == "tightening"
        n_relax += direction == "relaxation"
        say(f"  {label:<26s} {len(s):>6d} {direction:>11s} {len(s & m253):>9d} "
            f"{str(s < m253):>14s} {str(s > m253):>9s}")
    say()
    say(f"  tightenings (fewer genes than 253) : {n_tight}")
    say(f"  relaxations (more genes than 253)  : {n_relax}")


def block_antigens(genes, peptides):
    head(9, "Where the eight antigens sit on each metric (report section 05, S2)")
    tb = pd.read_csv(os.path.join(OUT, "supplementary_table_B_antigens.csv")).set_index("antigen")
    idx = genes.set_index("gene_id")

    best = genes[INDEL_COL].min()
    tie = int((genes[INDEL_COL] == best).sum())
    ceiling_pct = tb["indel_genome_percentile"].max()
    say(f"  best (lowest) value of {INDEL_COL} : {best} ({tie} genes tie)")
    say(f"  the best attainable indel genome percentile is therefore {ceiling_pct}%")
    say()
    say(f"  {'antigen':<9s} {'gene':<14s} {'indel value':>12s} {'indel pct':>10s} "
        f"{'at ceiling':>11s} {'cons pct':>9s} {'ident pct':>10s} {'pv':>6s} {'pb':>6s} {'pk':>6s}")
    at_ceiling = []
    for name, gid in ANTIGENS.items():
        r = idx.loc[gid]
        b = tb.loc[name]
        v = r[INDEL_COL]
        hit = bool(v == best)
        if hit:
            at_ceiling.append(name)
        say(f"  {name:<9s} {gid:<14s} {v:>12.4f} "
            f"{b['indel_genome_percentile']:>9.2f}% {str(hit):>11s} "
            f"{b['conservation_genome_percentile']:>8.2f}% "
            f"{b['human_identity_genome_percentile']:>9.2f}% "
            f"{str(bool(r['hom_homolog_in_pv'])):>6s} "
            f"{str(bool(r['hom_homolog_in_pb'])):>6s} "
            f"{str(bool(r['hom_homolog_in_pk'])):>6s}")
    say()
    say(f"  antigens exactly at the best indel value : {len(at_ceiling)} "
        f"({', '.join(at_ceiling)})")
    say(f"  antigens below it : "
        f"{', '.join(n for n in ANTIGENS if n not in at_ceiling)}")
    say()
    cons = tb["conservation_genome_percentile"]
    say(f"  conservation percentile: in the bottom 1% : "
        f"{int((cons < 1).sum())} ({', '.join(cons[cons < 1].index)})")
    say(f"                           in the bottom 2% : "
        f"{int((cons < 2).sum())} ({', '.join(cons[cons < 2].index)})")
    say(f"                           in the bottom 3% : {int((cons < 3).sum())}")
    say(f"                           sixth lowest is  : "
        f"{cons.sort_values().index[5]} at {cons.sort_values().iloc[5]:.2f}%")
    say()
    say("  Recorded orthologues. The filter tests a boolean column per species; a False")
    say("  is the absence of a recorded one-to-one orthologue in the shipped table, not")
    say("  evidence that no counterpart exists in that species.")
    for sp, col in [("P. vivax", "hom_homolog_in_pv"), ("P. berghei", "hom_homolog_in_pb"),
                    ("P. knowlesi", "hom_homolog_in_pk")]:
        missing = [n for n, g in ANTIGENS.items() if not bool(idx.loc[g, col])]
        say(f"    no recorded {sp:<12s} orthologue : {', '.join(missing) or 'none'}")


_M253_GENES = []


def main():
    genes, peptides = H.load_data()
    say("Report figures: every narrative number not produced by the other scripts.")
    say(f"pandas {pd.__version__}, numpy {np.__version__}")
    say()
    say("Harness validation (the app's own figures, from pepfilter_headless.validate):")
    H.validate(genes, peptides, verbose=False)
    say("  all six published figures reproduce")

    block_percentile_vs_count(genes, peptides)
    block_indel_arms(genes, peptides)
    block_r21_rebuttals(genes, peptides)
    block_corrections(genes, peptides)
    global _M253_GENES
    _M253_GENES = list(n_genes(H.m253(), genes, peptides).gene_id)
    block_conservation_boundary(genes)
    block_expression(genes)
    block_defects(genes, peptides)
    block_tightening(genes, peptides)
    block_antigens(genes, peptides)

    with open(REPORT_FIGURES, "w") as fh:
        fh.write("\n".join(_lines) + "\n")
    print(f"\nwrote {REPORT_FIGURES}")


if __name__ == "__main__":
    main()
