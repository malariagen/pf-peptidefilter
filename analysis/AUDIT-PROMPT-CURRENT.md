# Audit prompt — Pf-PeptideFilter R2 revision package

Paste everything below the line into a fresh window opened at
`/nfs/users/nfs_a/ab69/pf-peptidefilter`. In Claude Code, start with `/clear`.

Runnable in Codex or Claude. If run in both, do not let them coordinate.

---

You are auditing a bioinformatics analysis prepared for a manuscript revision — Reviewer 2 comments
**R2.1** (benchmark licensed and leading malaria vaccine antigens) and **R2.3** (add threshold
sensitivity analysis) on the Pf-PeptideFilter preprint, doi 10.64898/2025.12.15.694343. The work is
in `analysis/` and `analysis/output/`, untracked, on branch `add-manuscript-and-wiki`.

The package has been checked before and errors were found and fixed each time. **Every claim below
is therefore a claim, not a foundation.** Previous passes have shipped confident falsehoods —
including a percentile quoted as "bottom 2%" for a gene at 2.028%, a load-bearing claim that "only
conservation discriminates" when in fact three filters do, a percentile column briefly shipped
inverted, a results table silently destroyed by a careless string replacement, and a true finding
retracted after someone checked the wrong git branch. **Reproduce before you believe. Assert before
you edit.** If you find yourself agreeing with everything here, you have not re-derived enough.

Three parts, all equally weighted:

1. **Technical** — does the code do what it says; do the numbers reproduce; is the write-up
   internally consistent?
2. **Scientific** — fair comparisons, warranted inferences, disclosed limitations, honest framing?
3. **Against the brief** — does it answer what Jacob and the reviewer actually asked for?

## Environment

Python 3.12 **cannot** install `pandas==1.5.1` (no cp312 wheel); 3.10 and 3.11 can. Verified:
**Python 3.10.12, pandas 1.5.1, numpy 1.26.4**.

```bash
/usr/bin/python3 -m venv <scratch>/pf151          # /usr/bin/python3 is 3.10.12
<scratch>/pf151/bin/pip install "pandas==1.5.1" "numpy<2" matplotlib
/usr/bin/python3 -m venv <scratch>/wp
<scratch>/wp/bin/pip install -r analysis/requirements-report.txt   # weasyprint==70.0
```

WeasyPrint needs its own venv. Pagination is version-dependent, so establish the current page
count yourself with `pdfinfo` rather than trusting any figure quoted in a document — and never
cite a page number without naming the WeasyPrint version.

## Rules

- Every reported number must come from the app's own `filter_datasets()` in `pep/filter.py`, driven
  through `analysis/pepfilter_headless.py`. **Never reimplement filtering logic in standalone pandas
  as the source of a reported number** — an explicit early instruction from the user that still
  stands. Re-deriving directly from `data/*.csv.gz` is not merely allowed but *wanted*, as an
  independent cross-check.
- **Do not modify `pep/`, `config/`, `data/` or `requirements.txt`.** The provenance claim depends on
  those being byte-identical to live `main`, and `analysis/verify_provenance.sh` checks it —
  `requirements.txt` included. Analysis-only dependencies belong in
  `analysis/requirements-report.txt`.
- **Do not commit or push.**
- **`analysis/private/` must never be committed, pushed or copied off this machine.** It holds
  Jacob's email and planning notes, protected by a self-ignoring `.gitignore` (`*`) and by
  `.git/info/exclude`. Quote from it in your findings — that is what it is for — but never reproduce
  the email into a tracked file. If you move it, re-verify with `git check-ignore -v` and
  `git add --dry-run analysis/ | grep private`.
- When editing `report.html`, anchor every replacement on a string you have verified is present and
  assert the match count. A span replacement between two anchors will silently swallow everything
  between them — that is how a results table was destroyed. Note also that `§` in a Python literal
  has caused anchors to fail; prefer ASCII-only anchors.
- **The report and the email are read by the user's supervisor.** They should read as finished work.
  Do not add references to auditing, review rounds, revision history, or things having been "fixed"
  or "restored" — describe what the analysis found, not how the document evolved.
- Republishing the artifact is Claude-only (`Artifact` tool, pass
  `url: https://claude.ai/code/artifact/7ead69c8-37c7-46b6-9963-d54abff1233f`). **If you are Codex,
  skip it and say so.** Never create a second artifact.

## What exists

**Inputs, do not modify:** `pep/filter.py` (the engine; `filter_datasets()` is the single source of
every number), `pep/ui.py` (`UI_CONFIG` — every default, option list and slider range, i.e. what a
user can actually select), `config/filters.json` (18 filter definitions), and the two data files:
4,937 genes × 434 columns, 379,325 peptides × 45 columns. `manuscript/Manuscript.pdf` is the
preprint (`pdftotext -layout`). `analysis/private/jacob-email.md` is the authoritative statement of
what was asked for; `analysis/private/analysis-plan.md` is the attached plan, Steps 1–20, whose
M253 paragraph is the corrected version with the original preserved in a `<details>` block beneath.

**Scripts:** `pepfilter_headless.py` (substitutes a dict for `st.session_state`; provides
`app_defaults()`, `m253()`, `isolate()`, and `validate()`, which asserts published numbers before
any analysis runs), `run_sensitivity.py` (15 scenarios → Table A, retention matrix,
`rejected_scenarios.csv`, 15 frozen configs, 45 per-run catalogues), `run_antigens.py` (8-antigen
benchmark → Table B), `make_figure.py`, `run_interaction_checks.py` (exploratory),
`verify_conservation_grid.py`, `verify_provenance.sh`, `make_report_pdf.py`.

**Outputs:** Table A (22 cols × 15 rows), Table B (32 cols × 8 rows), the 8 × 15 retention matrix,
`rejected_scenarios.csv`, `interaction_checks.txt`, `COLUMN-DEFINITIONS.md`, `provenance.txt`, the
figure, and `report.html` → the PDF. Supervisor-facing: `DRAFT-EMAIL-TO-JACOB.md` and
`MANUSCRIPT-CORRECTIONS.md`.

## The mechanism everything depends on

Every "% of gene peptides that need to pass" control selects a **precomputed linear-interpolation
percentile column**. It does not count peptides (`pep/filter.py:287-315`):

| filter | column | test | at M253 |
|---|---|---|---|
| human identity | `blt_pident_p{human_id_gene_pc}` | `<= 80` | `_p75` |
| human length | `blt_length_p{human_id_gene_pc}` | `<= 15` | `_p75` |
| conservation | `hap_top_hap_freq_p{100 - haplotype_gene_pc}` | `>= 0.99` | `_p05` |
| indels | `max_indel_freq_type_any_fs_any_p{indel_gene_pc}` | `<= 0.05` | `_p100` |

Claimed consequences, all to be verified: the two human filters are separate marginal tests ANDed
together (the app never asks whether the *same* 75% of peptides passes both); a percentile is not a
count; only the indel rule (`_p100`, the maximum) has prose that is literally true; and the
conservation `100 − pc` inversion is correct and correctly used.

---

# Part 1 — technical

**T1. Re-derive independently.** Reimplement the five filters as direct masks over
`data/*.csv.gz`, never calling `filter_datasets()`, and compare **as sets** — not counts — against
the archived `run_outputs/`, for all 15 scenarios at gene *and* peptide level. Re-derive every
number `validate()` asserts the same way. Report any disagreement: the app-driven number is
authoritative, but a disagreement is a bug.

**T2. Replay the archive.** Feed all 15 frozen configs in `output/configs/` back through
`pepfilter_headless.run()`, confirm each matches its Table A row, and recompute all 22 columns
cell-by-cell.

**T3. The percentile mechanism.** Verify the suffix mapping, the direction of each test, and the
`100 − pc` inversion. Confirm an off-by-one would be detectable. Quantify how far the percentile
form diverges from its own label, in both directions, per filter.

**T4. Table B.** Confirm every `*_gene_rule_pct_tolerated` is the correct **floor** on the app's
39-value grid rather than the nearest; that pass/fail is monotone in the setting; that
`*_true_pct_peptides_passing` matches a direct count from the peptide table; that
`*_genome_percentile` runs in the stated direction and uses the right denominator; and that
`M253_exclusion_reason` lists **every** failing filter, not just the first.

**T5. Figure.** Every plotted value must match Table B; each reference line must be the correct
requirement; filled/hollow must agree with PASS/FAIL; sort order consistent across panels; axis
labels must not claim a quantity the values are not.

**T6. Render integrity.** WeasyPrint silently drops flex/grid content at page breaks — it has twice
removed whole sections from this report. After rendering, confirm all **eight** sections (the
plain-language walkthrough, §A, §01–§06) and all **17** settings in §02's completeness table are
present. `pdftotext` prefixes a form feed to any heading that starts a page, so match on heading
**text**, not line start. Before any republish, diff structural counts (`<section`, `<table`, `<tr`,
`.callout`, `.tablebox`, `<figure>`, `<li>`, `<blockquote>`) against the live version.

**T7. Consistency sweep.** Every number in `report.html`, `DRAFT-EMAIL-TO-JACOB.md` and
`MANUSCRIPT-CORRECTIONS.md` must match the current outputs. This is where stale text has repeatedly
survived a regeneration. Check especially: scenario and configuration counts, the
conservation-stringent arm, the two peptide definitions, percentile values and their denominators,
and any figure quoted in more than one document.

**T8. Provenance.** Re-run `verify_provenance.sh`; confirm `main` and `prod` have not moved and
re-verify rather than assume if either has. Confirm `git status --porcelain` shows no modified
tracked file and that `analysis/private/` is not stageable.

**T9. Verify the scripts are runnable and honest.** Every script in `analysis/` should run to
completion in the pinned environment. A crashing or vestigial script in a reproducibility bundle is
a defect. Check also that any claim in the email or report about a check having been performed
corresponds to something actually shipped in `analysis/` — a reader must be able to re-run it.

# Part 2 — scientific

For each, decide: **real / not real / needs a caveat / needs new analysis**, and say whether it
changes a conclusion.

**S1. Is the antigen benchmark a fair test?** The eight span two life stages (CSP, TRAP, CelTOS,
LSA3 pre-erythrocytic; AMA1, MSP1, RH5, EBA175 blood-stage) while the baseline requires liver-stage
day-4 expression. Which failures are stage-explained and which are not? Does the headline survive
substituting the sporozoite-stage filter, or removing expression and orthology altogether? Which
filters discriminate on stage-independent grounds, and how many antigens does each account for?

**S2. Is it circular?** The tool selects conserved genes; these antigens are known for variation.
Does the percentile-rank framing genuinely escape that, or only restate it? Is the whole gene set
the right comparator, or should it be matched for surface exposure or immune selection? Check what
each percentile scale can actually reach — ties at the best value cap the maximum attainable score,
and that changes how "extreme" and "mid-range" should be read.

**S3. Do the inferences outrun the computation?** Go through the report's scientific review, §03 and
§04 sentence by sentence and ask of each: does the computation establish this, or is it
interpretation? The report claims to mark that difference; check it does, consistently.

**S4. Threshold and sampling robustness.** Is 0.99 defensible given what can and cannot be
determined about sampling error from the shipped data? Is the fragility correctly quantified?

**S5. Expression filter validity.** Is CPM here parasite-relative within host-dominated libraries?
What does that do to cross-stage comparison, to the ≥1 CPM rule, and to the stringent arm? Is the
Day-4 justification stated in a way the data supports — and is the same criticism true *within* day
4 across its three replicates?

**S6. Geographic scope.** Conservation is African-only. Disclosed? Correctly scoped? Free of
uncited biological claims?

**S7. Design limitations.** One-at-a-time cannot detect compounding. Jacob explicitly scoped
combinatorics out, so the question is whether the limitation is *stated*, not whether it is fixed.
Is the interaction evidence correct and framed as exploratory?

**S8. Statistics and framing.** No inferential statistics are used — confirm none are needed and
none are implied. Is Jaccard versus containment interpreted correctly for nested sets of very
different sizes? Is every term a reader needs defined on first use? Does anything imply a verdict on
whether these antigens are good vaccine targets, which Jacob flagged as **Important**?

**S9. Anything else.** Every pass so far has found something new. Expect to.

# Part 3 — against the brief

`analysis/private/jacob-email.md` and `analysis/private/analysis-plan.md` are authoritative. Read
both line by line. Quote them; do not paraphrase. **Verify every quotation the report attributes to
them, and every quotation it attributes to the preprint** — one preprint quote was previously
silently normalised, so check them character by character against `pdftotext` output.

Points that carry consequences:

- The email rules out combinatorial search: *"we are not going to iterate through hundreds of
  combinations."* One-at-a-time is a decision. Do not propose a factorial sweep as a remedy.
- The email pre-authorises arbitrary thresholds: *"a bit arbitrary but just need something
  sensible."* So the specific values are not a defect — but **does every scenario have a stated
  rationale, and does it hold for the value actually used** rather than the one it replaced?
- Plan Step 3 forbids varying alignment length and the Joint rule; the analysis does both. Steps 4
  and 5 contain similar prohibitions on changing the gene-peptide proportion; check whether those
  are breached too, and whether each breach is disclosed or only some.
- Step 14 asks for the *percentage of peptides passing*. Confirm that quantity is present, correct,
  and clearly distinguished from the slider setting.
- Steps 2 and 19 require pre-registration; the email asks for empirical vetting. These conflict.
  Does the report state which it did, and avoid implying pre-registration?
- Steps 8, 9 and 15 specify **12 runs**; there are 15. Is the difference explained?
- Step 20 asks for a Results summary with five named components. One is drafted in the report.
  **This is the only text destined for the manuscript — audit it hardest.** Check all five
  components, every number, that it is marked as proposed rather than presented as manuscript text,
  and above all that it does not restate any error the rest of the document exists to correct.
- Also confirm: the email's four named outputs; Step 6's "prespecified sensitivity-analysis value"
  caveat; Step 12's identifier table; Step 17's heatmap and whatever replaced it; Step 18's
  interpretation rules; Step 19's reproducibility list.
- R2.1 asked about *"licensed and leading"* antigens. Only one of the eight underlies a licensed
  vaccine. Does the report imply otherwise?

# Part 4 — the plain-language walkthrough

The report opens with a plain-language section for a reader who did not build this and has not read
`filter.py`. Read it cold and judge it as that reader. It must explain what was asked, what the app
does (including the percentile gene-rule mechanism, with a concrete worked example), what the
baseline is, how the sensitivity analysis and antigen benchmark work, what was found, and what the
analysis cannot tell you. Verify every number in it, define-on-first-use for every term, and check
its counts agree with the rest of the document.

# Report back

1. **A table of every claim checked** — verified / refuted / could not verify, with the number *you*
   computed, never the number quoted at you here.
2. **Every discrepancy**, with the code or data that demonstrates it.
3. **Anything not on this list.**
4. **Anything on this list you judge to be wrong.**
5. **Your plain judgement**: is this sound enough to go into a manuscript revision, and what would
   you change first?

Fix what is clearly broken and within scope — stale numbers, wrong wording, a crashing script. For
anything that changes a **result** or a **scientific claim**, report it and propose the fix rather
than applying it unilaterally. Then re-render the PDF and verify per T6; Claude, republish to the
same artifact URL; Codex, skip that and say so.
