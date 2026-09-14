# Audit prompt — Pf-PeptideFilter R2 revision package

Paste everything below this line into Codex (or a fresh Claude session) with the repository
checked out. **Read-only review**: report findings, do not edit the package.

---

You are auditing a finished analysis package that answers two reviewer comments on a preprint.
Treat every claim in it as a hypothesis to be tested, not as context to be accepted. Your job is
to find things that are *wrong*, not to summarise.

## What exists

Repository: `malariagen/pf-peptidefilter`, branch `r2-sensitivity-analysis`.
The package is entirely under `analysis/`. The app itself — `pep/`, `config/`, `data/`,
`requirements.txt` — is unmodified and must stay that way.

| path | what it is |
|---|---|
| `analysis/report.html` → `analysis/output/R2-sensitivity-report.pdf` | the deliverable, 8 sections |
| `analysis/pepfilter_headless.py` | drives the app's own `filter_datasets()` from a dict instead of Streamlit widgets; `validate()` asserts six published figures |
| `analysis/run_sensitivity.py` | the 15 scenarios → Table A, retention matrix, rejected scenarios |
| `analysis/run_antigens.py` | per-antigen behaviour → Table B, percentile ranks, threshold sweeps |
| `analysis/run_interaction_checks.py` | multi-filter cross-checks (not archived scenario rows) |
| `analysis/make_figure.py`, `analysis/make_report_pdf.py` | figure, PDF render |
| `analysis/verify_provenance.sh`, `analysis/verify_conservation_grid.py` | provenance and grid checks |
| `analysis/MANUSCRIPT-CORRECTIONS.md` | six copy-ready wording corrections |
| `analysis/DRAFT-EMAIL-TO-JACOB.md` | covering note |
| `analysis/output/` | both supplementary tables, 15 frozen configs, 45 per-run catalogues, provenance record |

`analysis/private/` is excluded from git and is not part of the audit.

## Environment — do not skip this

The host Python may ship a pandas that the analysis was not written against. The pinned stack is
**Python 3.10, pandas 1.5.1, numpy 1.26.4**. Build an isolated environment from the repository's
own `requirements.txt` and run everything there. A previous auditor found the host had pandas 2.3.3
and correctly refused to trust results from it.

WeasyPrint needs a *separate* venv (`analysis/requirements-report.txt`, pinned to 70.0) because it
conflicts with the pinned pandas. Pagination moves between versions and the report is still being
edited, so **establish the current page count yourself with `pdfinfo`** rather than trusting any
figure quoted here.

## Rules

- Every number you report must come from the app's own `filter_datasets()` in `pep/filter.py`,
  driven through `analysis/pepfilter_headless.py`. **Never reimplement the filtering logic in
  standalone pandas** and then claim the package is wrong — that is how two earlier false alarms
  were produced. Recomputing a *gene metric* directly from `data/*.csv.gz` is fine and encouraged;
  recomputing the *filter chain* is not.
- Do not modify anything under `pep/`, `config/`, `data/` or `requirements.txt`.
- Do not commit, push, or create an artifact.
- Where you disagree with a number, give the figure you got, the exact configuration that produced
  it, and the command to reproduce it.

## The mechanism everything depends on

Get this right before auditing anything else, because most of the package's findings rest on it.

Each filter is measured per peptide, then summarised to the gene by selecting a **precomputed
linear-interpolation percentile column** — never by counting how many peptides passed:

| filter | column tested | test | at the published config |
|---|---|---|---|
| human identity | `blt_pident_p{human_id_gene_pc}` | `<= identity_percent` | `_p75 <= 80` |
| human length | `blt_length_p{human_id_gene_pc}` | `<= alignment_length` | `_p75 <= 15` |
| conservation | `hap_top_hap_freq_p{100 - haplotype_gene_pc}` | `>= strain_conservation` | `_p05 >= 0.99` |
| indels | `max_indel_freq_type_any_fs_any_p{indel_gene_pc}` | `<= indel_frequency` | `_p100 <= 0.05` |

Note the conservation inversion (`100 - pc`), which the other filters do not have. Note also that
the two human tests are **separate marginal tests ANDed together** — the app never asks whether the
*same* 75% of peptides passes both.

## Traps that have already caught people

Concrete, not hypothetical. Each of these produced a wrong result in this work at some point.

1. **`isolate()` leaks NaN.** "Neutralised" means held at a permissive threshold, not switched off.
   `NaN <= x` and `NaN >= x` are both False in pandas, so genes with a NaN in an *unrelated* metric
   still drop. An isolated indel run gives **3,451**; the indel metric alone gives **3,457**. Both
   are correct for different questions. Check which one a given claim is about before calling it an
   error.
2. **`"None"` is in pandas' default `na_values`.** The antigen column deliberately reads
   `none retained`, not `None`, so the headline R2.1 result does not round-trip to NaN.
3. **Percentile direction is easy to invert.** A "genome percentile" that reads 99.7 for CSP looks
   excellent but means "worse than 99.7% of the genome". The shipped columns count strictly-worse
   genes; CSP should read **0.3**.
4. **Ties cap the percentile scales.** Conservation reaches 99.6%, but human identity tops out at
   56.7% (2,138 genes tie at the best value) and indel frequency at 30.0% (3,457 tie). Any claim
   about an antigen being "mid-range" must account for what the scale can reach.
5. **Percentile denominators differ.** 4,930 genes carry a conservation value, 4,933 an identity
   value, all 4,937 an indel value. A footnote saying "against all 4,937" was wrong for two of the
   three.
6. **WeasyPrint does not fragment flex/grid containers across page breaks** — it silently drops the
   overflow. Print CSS forces block layout for paged media. Verify no section is missing from the
   PDF rather than assuming.
7. **`pdftotext` prefixes a form feed to headings that start a page**, so `grep '^Heading'` misses
   them. Two "missing section" alarms were artefacts of this.
8. **`.eyebrow` spans are letter-spaced by CSS**, so extracted text reads `A D D I T I O N A L` and
   defeats naive greps.
9. **Slider grids constrain what is selectable.** `strain_conservation` is min 0.0 / max 1.0 /
   step 0.01, so 0.995 and 0.999 cannot be chosen in the app. Gene-rule sliders offer 39 discrete
   values. A scenario that cannot be reproduced by clicking the published app is a defect.
10. **Anchored string replacement in `report.html` is dangerous.** Em dashes appear both as
    `&mdash;` and as literal `—`, and quotes both as entities and as literals. Assert before
    writing.

## Part 1 — technical

**T1.** Build the pinned environment. Run `analysis/pepfilter_headless.py`'s `validate()`. Confirm
all six published figures reproduce: 2,107 / 84,556 and 1,517 / 42,208 (conservation only at two
gene rules), 1,274 (conservation + human identity), 253 / 10,088 (the published configuration).

**T2.** Regenerate everything from the scripts. Compare the regenerated `output/` against what is
committed, cell by cell, not by eyeballing. Report any drift.

**T3.** Replay all 15 frozen configs in `output/configs/` and confirm each reproduces its own row in
Supplementary Table A, including both Jaccard columns and the direction-of-change counts.

**T4.** Independently verify the percentile mechanism from the raw tables: for a sample of genes,
confirm the app's verdict follows the percentile column and not a peptide count, and quantify the
disagreement genome-wide in both directions. The package claims **99 genes admitted** on fewer than
95% of their peptides and **2 rejected** despite meeting the count, and **zero** disagreement for
the indel rule. A previous auditor reported 4 rather than 2 for the second figure; settle it and say
which boundary convention you used.

**T5.** Verify the four departures from the plan are forced rather than chosen. Specifically:
identity should return 253 genes at every setting from ≤60% to ≤100%; the indel frequency threshold
should return 253 at 0.05, 0.04, 0.03, 0.02, 0.01 and 0.00; the indel type restrictions should
*loosen* the filter (280 / 260 / 294 / 297 against 253); conservation at 1.00 should empty the
catalogue. Departure 3 is a judgement call, not a data fact — say whether you think it is defensible.

**T6.** Verify the six manuscript corrections in §04. In particular the three readings of the
published 1,274 sentence: **1,274** as the app ran it, **158** run as written with 95% attached to
both filters, and **92** for the joint-subset reading. The 92 counts genes where the same ≥95% of
peptides clears conservation ≥0.99 *and* identity ≤80% *and* length ≤15 aa; check that the report
specifies all three, since the two-criterion reading gives 1,133.

**T7.** Run `verify_provenance.sh`. Confirm the app files used by the analysis are byte-identical to
live `main`, and that the differences against `prod` are all outside the filtering path.

**T8.** Confirm no tracked file outside `analysis/` is modified, and that `git add --dry-run .`
stages nothing from `analysis/private/` and nothing named `Comments for the authors.docx`.

**T9.** Run every script in `analysis/`. Any that fails, or any claim of the form "X was checked"
with no shipped file a reader could re-run, is a finding. One such claim was previously unsupported.

**T10.** Render the PDF in the WeasyPrint venv. Confirm all eight sections survive into it, no table
is truncated, and every one of the 52 internal `§NN` links resolves to a real section id.

## Part 2 — scientific

**S1.** Is the antigen benchmark fair? The panel mixes pre-erythrocytic and blood-stage antigens
against a liver-stage expression filter. The package argues the conclusion survives because the
sporozoite-stage filter gives 695 genes and still retains none, and removing expression and
orthology entirely gives 1,142 and still retains none. Test both.

**S2.** Is the benchmark circular? The tool selects for conservation and these antigens are famous
for variation. The package argues percentile ranks "blunt" the objection without dissolving it.
Assess whether that framing is honest, and whether the percentile scales are comparable across the
three criteria given the tie ceilings in trap 4.

**S3.** One-filter-at-a-time understates total sensitivity. The package reports 2,048 genes for all
five relaxed against an additive prediction of 958. Verify, and say whether the limitation is stated
plainly enough.

**S4.** Conservation comes from 8,492 African isolates of Pf7's 16,203. Is the scope limitation
stated where a reader would need it?

**S5.** The 0.99 threshold cannot be checked against sampling error — no coverage column ships.
The package reports 18 of the 253 within 0.001 of the cut. Verify, and assess whether the caveat is
proportionate.

**S6.** The day-4 expression default. The package argues day 2 is unusable because of library depth,
not parasite biology, citing one replicate totalling 15,682 CPM. Check the reasoning.

**S7.** Are there claims in the report that the computation does not support — interpretation
presented as output? Check the R2.1 "working as intended" passage especially.

**S8.** Anything scientifically wrong that nobody has flagged. This is the most valuable part of the
audit.

## Part 3 — against the brief

The plan and the originating email are private and not in the repository. §01 of the report quotes
them verbatim, item by item, in a numbered 23-row table. Audit the *internal consistency* of that
table: does each "Delivered, and where" actually exist at the path or section named, and does each
"What it showed" match what that section says? Flag any row whose evidence you cannot locate.

Four rows are marked "answered differently". Check that each has its reason stated and its
supporting run shipped.

## Part 4 — the plain-language walkthrough

Separately from the findings, write a walkthrough a non-specialist could follow: what the analysis
did, what it found, what changes to the paper it implies, and what remains uncertain. No jargon
without a gloss. This is used to check the report is comprehensible, so say plainly which parts of
it you found hard to follow and why.

## Report back

1. **Findings**, ranked by severity. For each: the claim, what you measured, the exact
   configuration, the command to reproduce, and whether it changes a result or only wording.
2. **The walkthrough** from Part 4.
3. **What you could not check**, and why.

Do not fix anything. If you believe a change is needed, say what and where.

One note on tone: the report and the covering note go to the colleague who wrote both the app and
the analysis plan. They must not acquire process language — no "audit round", no "this revision",
no explanations of how the app works, and no third-person references to him. If you suggest wording,
match that.
