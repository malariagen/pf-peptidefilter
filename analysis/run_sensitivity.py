"""
Threshold-sensitivity analysis (R2.3) and vaccine-antigen benchmark (R2.1).

Every figure in the outputs comes from a run of the app's own ``filter_datasets``.
Scenario definitions are frozen in SCENARIOS / REJECTED below.

Outputs (written to analysis/output/):
    supplementary_table_A_sensitivity.csv   one row per scenario
    supplementary_table_B_antigens.csv      per-antigen behaviour under each filter
    antigen_retention_matrix.csv            antigens x scenarios, retained/excluded
    rejected_scenarios.csv                  scenarios tested and discarded, with reasons
    filter_step_summaries.txt               the app's own per-filter gene counts
"""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pepfilter_headless as H  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUT, exist_ok=True)

IDENTITY_ONLY = "Use only identity filter (deactivate length filter)"

ANTIGENS = {                      # PlasmoDB identifiers, verified against gene_name
    "CSP":    "PF3D7_0304600",
    "AMA1":   "PF3D7_1133400",
    "RH5":    "PF3D7_0424100",
    "TRAP":   "PF3D7_1335900",
    "CelTOS": "PF3D7_1216600",
    "EBA175": "PF3D7_0731500",
    "MSP1":   "PF3D7_0930300",
    "LSA3":   "PF3D7_0220000",
}


def build_scenarios():
    app = H.app_defaults()
    m = H.m253(app)

    def v(**kw):
        cfg = dict(m)
        cfg.update(kw)
        return cfg

    return [
        # (label, filter varied, exact setting, params)
        ("APP",                     "-",            "application defaults",              dict(app)),
        ("M253",                    "-",            "published baseline",                dict(m)),

        ("HUMAN-RELAXED",           "human",        "alignment length <=18 aa",          v(alignment_length=18)),
        ("HUMAN-STRINGENT",         "human",        "alignment length <=12 aa",          v(alignment_length=12)),
        ("HUMAN-IDONLY-RELAXED",    "human",        "identity <=90%, length filter off", v(human_id_joint_rule=IDENTITY_ONLY, identity_percent=90)),
        ("HUMAN-IDONLY-BASELINE",   "human",        "identity <=80%, length filter off", v(human_id_joint_rule=IDENTITY_ONLY, identity_percent=80)),
        ("HUMAN-IDONLY-STRINGENT",  "human",        "identity <=70%, length filter off", v(human_id_joint_rule=IDENTITY_ONLY, identity_percent=70)),

        ("CONSERVATION-RELAXED",    "conservation", "haplotype frequency >=0.95",        v(strain_conservation=0.95)),
        ("CONSERVATION-STRINGENT",  "conservation", "99% of peptides at >=0.99 (gene rule)"      , v(haplotype_gene_pc=99)),

        ("INDEL-RELAXED-FREQ",      "indel",        "indel frequency <=0.10",            v(indel_frequency=0.10)),
        ("INDEL-RELAXED-RULE",      "indel",        "95% of peptides must pass",         v(indel_gene_pc=95)),

        ("EXPRESSION-RELAXED",      "expression",   ">=1 CPM in >=2 replicates (D4)",    v(expression_rule="At least two replicates")),
        ("EXPRESSION-STRINGENT",    "expression",   ">=5 CPM in all replicates (D4)",    v(cpm_cutoff=5)),

        ("ORTHOLOGY-RELAXED",       "orthology",    "no orthology requirement",          v(homology_species=[])),
        ("ORTHOLOGY-STRINGENT",     "orthology",    "P. vivax AND P. berghei",           v(homology_species=["P. vivax", "P. berghei"])),
    ]


def build_rejected():
    """Scenarios from the original plan that were tested and discarded."""
    m = H.m253()

    def v(**kw):
        cfg = dict(m)
        cfg.update(kw)
        return cfg

    return [
        ("HUMAN-RELAXED (identity 90%)",   v(identity_percent=90),
         "no effect: under the joint rule the <=15 aa length criterion is binding at every identity value"),
        ("HUMAN-STRINGENT (identity 70%)", v(identity_percent=70),
         "no effect: same reason"),
        ("INDEL-STRINGENT (<=0.01)",       v(indel_frequency=0.01),
         "no effect: the gene-level indel metric is a step function, identical for any cutoff 0.00-0.05"),
        ("CONSERVATION-STRINGENT (=1.00)", v(strain_conservation=1.00),
         "degenerate: empties the catalogue (0 genes)"),
        ("CONSERVATION-STRINGENT (=0.995)", v(strain_conservation=0.995),
         "off-grid: the app's haplotype-frequency slider has a 0.01 step, so 0.995 is not "
         "selectable; superseded by tightening the gene rule to 99%, which is on-grid"),
    ]


def jaccard(a, b):
    union = a | b
    return len(a & b) / len(union) if union else float("nan")


def main():
    genes, peptides = H.load_data()

    print("Validating the harness against the preprint:")
    H.validate(genes, peptides)
    print()

    scenarios = build_scenarios()
    results, step_summaries = {}, {}

    cfg_dir = os.path.join(OUT, "configs")
    run_dir = os.path.join(OUT, "run_outputs")
    os.makedirs(cfg_dir, exist_ok=True)
    os.makedirs(run_dir, exist_ok=True)

    for label, filt, setting, params in scenarios:
        gsum, _, fgenes, fpeptides, fmap = H.run(params, genes, peptides)
        gset = set(fgenes.gene_id)
        pset = set(peptides.loc[peptides.gene_id.isin(gset), "peptide_id"])
        results[label] = dict(filter=filt, setting=setting, genes=gset,
                              peptides=pset, filter_map=fmap)
        step_summaries[label] = gsum

        # frozen settings for this run, exactly as handed to the app
        with open(os.path.join(cfg_dir, f"config_{label}.json"), "w") as fh:
            json.dump({"scenario": label, "filter_varied": filt,
                       "exact_setting": setting, "settings": params}, fh, indent=2)

        # gene- and peptide-level catalogues
        fgenes[["gene_id", "gene_name", "chromosome", "start", "end", "num_peptides"]] \
            .to_csv(os.path.join(run_dir, f"{label}_genes.csv.gz"), index=False)
        peptides.loc[peptides.gene_id.isin(gset),
                     ["peptide_id", "gene_id", "gene_name", "start", "end", "ref_3D7_peptide"]] \
            .to_csv(os.path.join(run_dir, f"{label}_peptides.csv.gz"), index=False)
        # the app's own independently filtered peptide set, for completeness
        fpeptides[["peptide_id", "gene_id"]] \
            .to_csv(os.path.join(run_dir, f"{label}_peptides_filtered_independently.csv.gz"),
                    index=False)

        print(f"  {label:24s} {len(gset):5d} genes  {len(pset):6d} peptides", flush=True)

    ref_g, ref_p = results["M253"]["genes"], results["M253"]["peptides"]
    app_g, app_p = results["APP"]["genes"], results["APP"]["peptides"]

    rows = []
    for label, _, _, _ in scenarios:
        r = results[label]
        G, P = r["genes"], r["peptides"]
        retained = [a for a, gid in ANTIGENS.items() if gid in G]
        rows.append({
            "scenario": label,
            "filter_varied": r["filter"],
            "exact_setting": r["setting"],
            "genes": len(G),
            "peptides": len(P),
            "M253_genes_shared": len(G & ref_g),
            "M253_genes_pct_retained": round(100 * len(G & ref_g) / len(ref_g), 1) if ref_g else float("nan"),
            "M253_gene_jaccard": round(jaccard(G, ref_g), 3),
            "M253_peptides_shared": len(P & ref_p),
            "M253_peptides_pct_retained": round(100 * len(P & ref_p) / len(ref_p), 1) if ref_p else float("nan"),
            "M253_peptide_jaccard": round(jaccard(P, ref_p), 3),
            "APP_genes_shared": len(G & app_g),
            "APP_genes_pct_retained": round(100 * len(G & app_g) / len(app_g), 1) if app_g else float("nan"),
            "APP_gene_jaccard": round(jaccard(G, app_g), 3),
            "APP_peptides_shared": len(P & app_p),
            "APP_peptides_pct_retained": round(100 * len(P & app_p) / len(app_p), 1) if app_p else float("nan"),
            "APP_peptide_jaccard": round(jaccard(P, app_p), 3),
            "genes_gained_vs_M253": len(G - ref_g),
            "genes_lost_vs_M253": len(ref_g - G),
            "peptides_gained_vs_M253": len(P - ref_p),
            "peptides_lost_vs_M253": len(ref_p - P),
            # NOT the string "None": that is in pandas' default na_values list, so the column
            # round-trips to NaN on read and the headline R2.1 result -- no benchmark antigen
            # retained under any scenario -- silently becomes missing data for anyone loading the
            # supplementary table in pandas.
            "benchmark_antigens_retained": ", ".join(retained) if retained else "none retained",
        })
    table_a = pd.DataFrame(rows)
    table_a.to_csv(os.path.join(OUT, "supplementary_table_A_sensitivity.csv"), index=False)

    matrix = pd.DataFrame({
        label: {a: ("retained" if gid in results[label]["genes"] else "excluded")
                for a, gid in ANTIGENS.items()}
        for label, _, _, _ in scenarios
    })
    matrix.to_csv(os.path.join(OUT, "antigen_retention_matrix.csv"))

    rej_rows = []
    for label, params, reason in build_rejected():
        _, _, fg, _, _ = H.run(params, genes, peptides)
        rej_rows.append({"scenario_as_planned": label, "genes": len(fg),
                         "genes_vs_M253": len(fg) - len(ref_g), "reason_rejected": reason})
        print(f"  [rejected] {label:34s} {len(fg):5d} genes", flush=True)
    pd.DataFrame(rej_rows).to_csv(os.path.join(OUT, "rejected_scenarios.csv"), index=False)

    with open(os.path.join(OUT, "filter_step_summaries.txt"), "w") as fh:
        for label, _, setting, _ in scenarios:
            fh.write(f"### {label}  ({setting})\n")
            fh.write(step_summaries[label].to_string(index=False) + "\n\n")

    pd.set_option("display.width", 250)
    print("\n" + table_a[["scenario", "exact_setting", "genes", "peptides",
                          "M253_gene_jaccard", "genes_gained_vs_M253",
                          "genes_lost_vs_M253", "benchmark_antigens_retained"]].to_string(index=False))
    return table_a


if __name__ == "__main__":
    main()
