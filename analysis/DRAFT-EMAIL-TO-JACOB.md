# Draft covering note — edit freely before sending

**Subject:** R2.1 / R2.3 analysis done — report attached, six wording corrections needed

---

Hi Jacob,

The combined R2.1 / R2.3 run is finished. Report attached; it answers each item in the plan
in turn, so I'll keep this short and let it speak for the results.

**Headlines.**

- **R2.1 — none of the eight antigens is retained, under any of the fifteen configurations.**
  Strain conservation is the criterion responsible in every case. Six of the eight are in the
  bottom 2% of the genome for conservation. It isn't a life-stage artefact: the sporozoite-stage
  filter still retains none, and removing expression and orthology altogether still retains none.
- **R2.3 — the catalogue ranges 176 to 1,142 genes** around the published 253. Most sensitive to
  expression and conservation, least to indel frequency.
- **Four of the plan's ten scenarios had to be replaced.** Three measured nothing at all and one
  emptied the catalogue — the case you anticipated. Replacements are tested and selectable in the
  app; the rejected ones are kept with reasons in `rejected_scenarios.csv`. §01 of the report has
  the detail, and there are four departures from the plan overall, each with the run that forced it.

**Six wording corrections, and five are already in print.** This is the part I'd flag hardest.
Every *count* in the preprint reproduces exactly — it's the prose describing how the counts were
reached that's wrong, and it traces to one thing: the gene-level controls are labelled "% of
peptides that must pass", but the code tests a precomputed percentile column and never counts
peptides. §04 of the report has all six as a numbered checklist, and
`analysis/MANUSCRIPT-CORRECTIONS.md` has them copy-ready — sentence as published, replacement
beneath.

**Three things I need from you.**

1. The Day-4 citation for the expression default. The report shows the app's own expression table
   supports the choice, so an external reference may not be strictly necessary, but you'll know
   better.
2. Sign-off on the Step 20 Results passage. It's drafted in §04 and marked as proposed — please
   treat it as a starting point rather than finished prose.
3. A decision on where `analysis/` lives. It's committed to `r2-sensitivity-analysis` locally but
   not pushed. If the revision cites it, it needs a citable home — and note the report quotes your
   email verbatim in places, so worth a look before anything goes public.

**One limitation worth a sentence in the Discussion:** conservation is measured in African isolates
only — 8,492 of Pf7's 16,203 genomes — so the catalogue is an African-frequency catalogue rather
than a global one. The sensitivity analysis can't test this, because the app ships a single
African-derived metric.

Everything reproduces from the app's own `filter_datasets()`; nothing under `pep/`, `config/`,
`data/` or `requirements.txt` was touched. An independent reviewer re-derived all fifteen gene
sets, both peptide definitions and every numeric cell of Supplementary Table A with a separate
implementation, and they match.

Happy to walk through any of it.

Andrew

---

## Notes for me, not for sending

- Attach: `analysis/output/R2-sensitivity-report.pdf` (36 pp), plus
  `supplementary_table_A_sensitivity.csv` and `supplementary_table_B_antigens.csv` if he wants the
  raw tables.
- `MANUSCRIPT-CORRECTIONS.md` numbering matches §04's checklist 1–6, so the two can be worked side
  by side.
- Don't send `analysis/private/` or the reviewer `.docx` — both are excluded from git, and the
  `.docx` is a confidential peer review.
- If he asks why the human filter varies alignment length rather than identity: identity returns
  exactly 253 genes at every setting from ≤60% to ≤100%, so the scenario as specified in Step 3
  measures nothing. That's departure 1 and the one most worth discussing.
