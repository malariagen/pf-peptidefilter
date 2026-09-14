"""
Vaccine-antigen benchmark (R2.1): how the eight requested antigens behave under each
filter of Pf-PeptideFilter.

Nothing here reimplements a filter condition. Each PASS/FAIL is a run of the app with
every other filter held at a permissive threshold, and each "% of peptides passing" is
read out by sweeping the app's own *gene filtering rule* slider to find the highest
percentage-of-peptides requirement the gene still satisfies.

Two things that reading invites and the numbers do not support:

  * "neutralised" is not "switched off". The permissive filters are still applied, so a
    gene with a NaN metric fails them anyway -- see the caveat on
    ``pepfilter_headless.isolate``. Immaterial to all eight antigens (none has a NaN
    metric), but an isolated gene count is not the filter in isolation.
  * the swept value is a *slider setting*, not the share of peptides passing. It is the
    correct floor on the slider's own 39-value grid for the app's percentile column --
    that was checked for all 8 x 3 cells -- but the gene rule tests a linearly
    interpolated percentile, so it can sit either side of the grid floor of the true
    peptide fraction. It does for 4 of those 24 cells: TRAP indel reads 93 where the raw
    fraction is 92.98% (floor 92), MSP1 indel 92 where it is 91.86% (floor 91), EBA175
    human 85 where only 83.33% of peptides pass both conditions (floor 82), and CelTOS
    conservation 40 where 50.00% pass (floor 50). No PASS/FAIL verdict changes.

Interpretation: this measures how established antigens behave under the specific design
criteria implemented by Pf-PeptideFilter. It is not an assessment of whether they are
good or bad vaccine targets.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pepfilter_headless as H  # noqa: E402
from run_sensitivity import ANTIGENS, OUT  # noqa: E402

# the values the app's "% of gene peptides that need to pass" sliders actually offer
PC_OPTIONS = H.UI_CONFIG["haplotype_gene_pc"]["options"]

# gene-rule slider belonging to each filter
RULE_KEY = {"human": "human_id_gene_pc",
            "conservation": "haplotype_gene_pc",
            "indel": "indel_gene_pc"}

REQUIREMENT = {"human": 75, "conservation": 95, "indel": 100}   # the M253 settings


def sweep_gene_rule(filt, genes, peptides):
    """Highest '% of peptides must pass' setting at which each gene still survives.

    Equivalent to the proportion of the gene's peptides that meet the filter, but
    obtained purely by driving the app.
    """
    key = RULE_KEY[filt]
    best = {}
    for pc in PC_OPTIONS:
        cfg = H.isolate(filt)
        cfg[key] = pc
        _, _, fg, _, _ = H.run(cfg, genes, peptides)
        for gid in set(fg.gene_id):
            if best.get(gid, -1) < pc:
                best[gid] = pc
    return best


def sweep_threshold(genes, peptides, key, values, filt):
    """Most stringent threshold value at which each gene still survives."""
    best = {}
    for val in values:
        cfg = H.isolate(filt)
        cfg[key] = val
        _, _, fg, _, _ = H.run(cfg, genes, peptides)
        for gid in set(fg.gene_id):
            if best.get(gid, -1) < val:
                best[gid] = val
    return best


def _fast_peptides(peptides, genes):
    """A small peptide frame for the sweeps.

    Gene-level results are computed from the gene table alone -- the peptide table is
    filtered separately and never feeds back -- so shrinking it does not change any gene
    set. Asserted below rather than assumed.
    """
    small = peptides.head(1000).copy()
    # Checked across the configurations actually swept, not just M253: the gene rule, the
    # threshold and the isolation settings all take different code paths in filter_datasets.
    probes = [("m253", H.m253())]
    for filt in ("human", "conservation", "indel"):
        probes.append((f"isolate-{filt}", H.isolate(filt)))
    hard = H.isolate("conservation"); hard["haplotype_gene_pc"] = 99
    probes.append(("isolate-conservation-pc99", hard))
    for label, cfg in probes:
        _, _, full, _, _ = H.run(cfg, genes, peptides)
        _, _, cheap, _, _ = H.run(cfg, genes, small)
        assert set(full.gene_id) == set(cheap.gene_id), (
            f"gene results depend on the peptide table under {label}; "
            "the sweep shortcut is not safe")
    return small


def true_peptide_fractions(peptides, gid):
    """The actual proportion of a gene's peptides that satisfy each peptide-level condition.

    This is what plan Step 14 asks for -- "the relevant percentage of peptides passing" -- and it
    is NOT the same as the ``*_gene_rule_pct_tolerated`` columns, which are settings on the app's
    39-value gene-rule slider. The gene rule tests a linearly interpolated percentile, so the
    setting can sit either side of this fraction; the two disagree in 4 of the 24 antigen x filter
    cells. Both are reported, because they answer different questions.

    Read straight off the peptide table at the M253 peptide-level thresholds. This is a
    descriptive read of the input data, not a reimplementation of a filter: no gene is retained or
    dropped on the strength of it, and every PASS/FAIL in this table still comes from
    filter_datasets().
    """
    sp = peptides[peptides.gene_id == gid]
    if not len(sp):
        return {}
    return {
        "human": float(((sp.blt_pident <= 80) & (sp.blt_length <= 15)).mean() * 100),
        "conservation": float((sp.hap_top_hap_freq >= 0.99).mean() * 100),
        "indel": float((sp.max_indel_freq_type_any_fs_any <= 0.05).mean() * 100),
    }


def genome_percentile_ranks(genes):
    """Each gene's rank among all 4,937 genes on each gene-level metric, as a percentile.

    Answers a question pass/fail against one threshold cannot: not just "does this antigen fail
    the conservation filter" but "how unusual is it?". Without this the benchmark is close to
    circular -- the tool selects for conserved genes, the eight antigens are known for antigenic
    variation, so their exclusion is nearly guaranteed by construction. A percentile rank turns a
    uniform block of "excluded" into a graded result.

    Directionality is handled per metric: higher haplotype frequency is better, lower identity,
    alignment length and indel frequency are better. Every column reads the same way -- "this gene
    is better than X% of all 4,937 genes on this metric" -- so a low number is a bad gene.

    Descriptive statistics over the gene table, not a filter -- no retention decision uses it.
    """
    g = genes.set_index("gene_id")
    spec = {
        "conservation": ("hap_top_hap_freq_p05", False),
        "human_identity": ("blt_pident_p75", True),
        "human_length": ("blt_length_p75", True),
        "indel": ("max_indel_freq_type_any_fs_any_p100", True),
    }
    out = {}
    for name, (col, lower_is_better) in spec.items():
        v = g[col].dropna()
        # "better than X% of all genes", stated in that direction for every metric. Counting
        # strictly-worse genes rather than using rank(pct=) keeps ties out of the numerator and
        # makes the direction explicit: an inverted comparison here would read 99.7 for a gene
        # that is in fact worse than 99.7% of the genome.
        worse = (lambda x: (v < x).mean()) if not lower_is_better else (lambda x: (v > x).mean())
        out[name] = v.map(lambda x: round(float(worse(x)) * 100, 1))
    return out


def main():
    genes, peptides = H.load_data()
    H.validate(genes, peptides)
    print()
    sweep_peptides = _fast_peptides(peptides, genes)
    print("  sweep shortcut verified: gene sets identical with a reduced peptide table\n")

    # --- pass/fail for each filter in isolation, at its M253 setting ---------
    passes = {}
    for filt in ["human", "conservation", "indel", "expression", "orthology"]:
        _, _, fg, _, _ = H.run(H.isolate(filt), genes, peptides)
        passes[filt] = set(fg.gene_id)
        print(f"  isolated {filt:13s} {len(passes[filt]):5d} genes pass", flush=True)

    # --- how far each gene is from the gene-rule requirement -----------------
    rule_headroom = {}
    for filt in ["human", "conservation", "indel"]:
        rule_headroom[filt] = sweep_gene_rule(filt, genes, sweep_peptides)
        print(f"  swept gene rule for {filt}", flush=True)

    # --- most stringent conservation threshold each gene tolerates ----------
    # the app's strain-conservation slider is min 0.0, max 1.0, step 0.01, so the grid
    # below is exactly what a user can select -- 0.995 used to be appended here, which is
    # off that grid and so not a value the app can produce
    cons_step = H.UI_CONFIG["strain_conservation"]["step"]
    cons_grid = [round(i * cons_step, 3)
                 for i in range(int(round(H.UI_CONFIG["strain_conservation"]["max_value"]
                                         / cons_step)) + 1)]
    cons_thresh = sweep_threshold(
        genes, sweep_peptides, "strain_conservation", cons_grid, "conservation")
    print("  swept conservation threshold", flush=True)

    # --- expression evidence and orthology, via the app's own options -------
    expr_level = {}
    for rule, name in [("At least one replicate", 1), ("At least two replicates", 2),
                       ("All replicates", 3)]:
        cfg = H.isolate("expression")
        cfg["expression_rule"] = rule
        _, _, fg, _, _ = H.run(cfg, genes, sweep_peptides)
        for gid in set(fg.gene_id):
            expr_level[gid] = max(expr_level.get(gid, 0), name)
    # --- final status under the two reference configurations (plan Step 14/16) ---
    _, _, fg_app, _, _ = H.run(H.app_defaults(), genes, peptides)
    _, _, fg_m253, _, _ = H.run(H.m253(), genes, peptides)
    app_set, m253_set = set(fg_app.gene_id), set(fg_m253.gene_id)

    ortho = {}
    for sp in ["P. vivax", "P. berghei", "P. knowlesi"]:
        cfg = H.isolate("orthology")
        cfg["homology_species"] = [sp]
        _, _, fg, _, _ = H.run(cfg, genes, sweep_peptides)
        ortho[sp] = set(fg.gene_id)
    print("  swept expression and orthology", flush=True)

    # --- assemble ------------------------------------------------------------
    gene_names = genes.set_index("gene_id")["gene_name"].to_dict()
    npep = peptides.groupby("gene_id").size().to_dict()

    pctile = genome_percentile_ranks(genes)

    rows = []
    for name, gid in ANTIGENS.items():
        true_pct = true_peptide_fractions(peptides, gid)
        failed = [f for f in ["human", "conservation", "indel", "expression", "orthology"]
                  if gid not in passes[f]]
        rows.append({
            "antigen": name,
            "gene_id": gid,
            "gene_name_in_app": gene_names.get(gid),
            "n_peptides": npep.get(gid, 0),

            "human_gene_rule_pct_tolerated": rule_headroom["human"].get(gid, 0),
            "human_true_pct_peptides_passing": round(true_pct.get("human", float("nan")), 1),
            "human_requirement_pct": REQUIREMENT["human"],
            "human": "PASS" if gid in passes["human"] else "FAIL",

            "conservation_gene_rule_pct_tolerated": rule_headroom["conservation"].get(gid, 0),
            "conservation_true_pct_peptides_passing": round(true_pct.get("conservation", float("nan")), 1),
            "conservation_requirement_pct": REQUIREMENT["conservation"],
            "conservation_max_threshold_tolerated": cons_thresh.get(gid, float("nan")),
            "conservation_requirement_threshold": 0.99,
            "conservation": "PASS" if gid in passes["conservation"] else "FAIL",

            "indel_gene_rule_pct_tolerated": rule_headroom["indel"].get(gid, 0),
            "indel_true_pct_peptides_passing": round(true_pct.get("indel", float("nan")), 1),
            "indel_requirement_pct": REQUIREMENT["indel"],
            "indel": "PASS" if gid in passes["indel"] else "FAIL",

            "d4_replicates_expressed": expr_level.get(gid, 0),
            "d4_replicates_required": 3,
            "expression_D4": "PASS" if gid in passes["expression"] else "FAIL",

            "ortholog_P_vivax": gid in ortho["P. vivax"],
            "ortholog_P_berghei": gid in ortho["P. berghei"],
            "ortholog_P_knowlesi": gid in ortho["P. knowlesi"],
            "orthology_Pvivax": "PASS" if gid in passes["orthology"] else "FAIL",

            "conservation_genome_percentile": pctile["conservation"].get(gid),
            "human_identity_genome_percentile": pctile["human_identity"].get(gid),
            "indel_genome_percentile": pctile["indel"].get(gid),

            "APP_status": "retained" if gid in app_set else "excluded",
            "M253_status": "retained" if gid in m253_set else "excluded",
            "n_filters_failed": len(failed),
            "M253_exclusion_reason": ", ".join(failed) if failed else "-",
        })

    table_b = pd.DataFrame(rows)
    table_b.to_csv(os.path.join(OUT, "supplementary_table_B_antigens.csv"), index=False)

    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 60)
    print("\n" + table_b[[
        "antigen", "gene_id", "n_peptides",
        "human_gene_rule_pct_tolerated", "human_true_pct_peptides_passing", "human",
        "conservation_gene_rule_pct_tolerated", "conservation_true_pct_peptides_passing",
        "conservation_max_threshold_tolerated", "conservation",
        "indel_gene_rule_pct_tolerated", "indel_true_pct_peptides_passing", "indel",
        "conservation_genome_percentile",
        "d4_replicates_expressed", "expression_D4",
        "ortholog_P_vivax", "APP_status", "M253_status",
        "M253_exclusion_reason"]].to_string(index=False))
    return table_b


if __name__ == "__main__":
    main()
