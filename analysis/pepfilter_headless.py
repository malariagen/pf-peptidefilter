"""
Headless driver for the Pf-PeptideFilter app.

The app collects widget values into ``st.session_state`` and passes them to
``pep.filter.filter_datasets``.  That function touches streamlit in exactly two
places (both plain reads of ``st.session_state``), so substituting a dict for the
streamlit module leaves the filtering logic untouched.

Every number produced downstream therefore comes from the app's own code, its own
``config/filters.json`` and its own ``data/*.csv.gz`` -- not from a reimplementation.

Validated against four figures published in the preprint; see ``validate()``.
"""
import copy
import json
import os
import sys
import types

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# --- stand in for streamlit ------------------------------------------------
# pep.filter uses st.session_state as a plain mapping and nothing else.
_st = types.ModuleType("streamlit")
_st.session_state = {}
_st.cache_data = lambda *a, **k: (lambda f: f)
sys.modules.setdefault("streamlit", _st)

sys.path.insert(0, REPO)
from pep.filter import filter_datasets  # noqa: E402
from pep.ui import UI_CONFIG            # noqa: E402

with open(os.path.join(REPO, "config", "filters.json")) as fh:
    FILTER_CONFIG = json.load(fh)

COMPONENT_KEYS = list(UI_CONFIG.keys())


def app_defaults():
    """The settings a user sees on first opening the app, straight from UI_CONFIG."""
    out = {}
    for key, cfg in UI_CONFIG.items():
        if "index" in cfg:                                    # radio
            out[key] = cfg["options"][cfg["index"]]
        elif isinstance(cfg.get("default"), list):            # multiselect
            out[key] = list(cfg["default"])
        else:                                                 # slider / select_slider
            out[key] = cfg["default"]
    return out


def load_data():
    genes = pd.read_csv(os.path.join(REPO, "data", "gene-metrics-filtering.csv.gz"),
                        low_memory=False)
    peptides = pd.read_csv(os.path.join(REPO, "data", "peptide-metrics-filtering.csv.gz"),
                           low_memory=False)
    return genes, peptides


def run(params, genes, peptides):
    """Call the app's filter_datasets with `params` in place of the widget state."""
    _st.session_state.clear()
    _st.session_state.update(copy.deepcopy(params))
    return filter_datasets(genes, peptides, FILTER_CONFIG, COMPONENT_KEYS)


# --- named configurations --------------------------------------------------

def m253(base=None):
    """The published catalogue: app defaults + Day-4 expression + P. vivax orthology."""
    cfg = dict(base or app_defaults())
    cfg["expression_stage"] = ["Liver Stage Day 4"]
    cfg["homology_species"] = ["P. vivax"]
    return cfg


# Thresholds that let every gene through, so a single filter can be isolated
# using nothing but settings the app itself offers.
NEUTRAL = {
    "human_id_gene_pc": 0, "identity_percent": 100, "alignment_length": 20,
    "haplotype_gene_pc": 0, "strain_conservation": 0.0,
    "indel_gene_pc": 0, "indel_frequency": 1.0,
    "expression_stage": [], "homology_species": [],
}


def isolate(filter_name, base=None):
    """Settings with every filter neutralised except `filter_name`, held at its M253 value.

    Caveat: "neutralised" means held at a threshold every gene clears, not switched off --
    the app offers no way to remove a filter, so the permissive filters are still applied.
    A gene whose metric is NaN therefore still fails them, because `NaN <= x` and
    `NaN >= x` are both False in pandas. Seven of the 4,937 genes have a NaN in at least
    one gene-level metric, so an isolated run can report up to six fewer genes than the
    filter under test alone would admit (indel: 3,451 here vs 3,457 by the metric alone;
    human: 3,693 vs 3,695; conservation: unaffected). None of the seven is in the
    published 253 and none is a benchmark antigen, so no reported verdict changes -- but
    an isolated count is a count under that filter *plus* NaN exclusion, not the filter
    in isolation.
    """
    ref = m253(base)
    cfg = dict(ref)
    cfg.update(NEUTRAL)
    keep = {
        "human":        ["human_id_gene_pc", "identity_percent", "alignment_length"],
        "conservation": ["haplotype_gene_pc", "strain_conservation"],
        "indel":        ["indel_gene_pc", "indel_frequency"],
        "expression":   ["expression_stage"],
        "orthology":    ["homology_species"],
    }[filter_name]
    for k in keep:
        cfg[k] = ref[k]
    return cfg


def validate(genes, peptides, verbose=True):
    """Reproduce four numbers published in the preprint. Raises if any disagree."""
    d = app_defaults()

    cons_only = dict(d); cons_only.update(NEUTRAL)
    cons_only["haplotype_gene_pc"] = 95; cons_only["strain_conservation"] = 0.99

    cons_99 = dict(cons_only); cons_99["haplotype_gene_pc"] = 99

    # human + conservation, with the indel filter neutralised (the preprint's 1,274
    # is the count before indels are applied)
    hum_cons = dict(d)
    hum_cons["indel_gene_pc"] = 0
    hum_cons["indel_frequency"] = 1.0

    # Peptide counts are the "belonging to the retained genes" definition, which is what
    # the paper uses throughout: the independently filtered peptide set is 326,701 for both
    # of the conservation-only rows, because peptide filtering strips the percentile suffix
    # and so cannot see the gene rule at all.
    checks = [
        ("conservation only, 95% of peptides", cons_only, 2107, 84556),
        ("conservation only, 99% of peptides", cons_99,   1517, 42208),
        ("conservation + human identity",      hum_cons,  1274, None),
        ("all five filters (M253)",            m253(d),    253, 10088),
    ]
    ok = True
    for label, params, want_genes, want_peps in checks:
        _, _, fg, _, _ = run(params, genes, peptides)
        got_genes = len(fg)
        got_peps = int(peptides["gene_id"].isin(fg.gene_id).sum())
        good = got_genes == want_genes and (want_peps is None or got_peps == want_peps)
        ok &= good
        if verbose:
            exp = f"{want_genes}" + (f" / {want_peps}" if want_peps else "")
            got = f"{got_genes}" + (f" / {got_peps}" if want_peps else "")
            print(f"  [{'ok' if good else 'FAIL'}] {label:38s} expected {exp:>14s}  got {got:>14s}")
    if not ok:
        raise AssertionError("harness does not reproduce the published numbers")
    return ok
