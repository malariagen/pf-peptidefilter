# Manuscript corrections — copy-ready

Six wording corrections for the revision. **Every count in the paper reproduces exactly** — all six
published figures are now asserted by `analysis/pepfilter_headless.validate()` before any analysis
runs. What needs fixing is the prose describing *how* those counts were obtained: five sentences
describe a configuration that is not the one that produces the published numbers.

Five of the six are already in print. That is worth saying plainly in the response to R2, because
R2.3 is a reproducibility question and these are precisely the sentences a reviewer re-deriving the
work would hit first.

Every number below was re-derived through the app's own `filter_datasets()` on the repository's own
`config/filters.json` and `data/*.csv.gz` (Python 3.10.12, pandas 1.5.1, numpy 1.26.4), against a
working tree verified byte-identical to live `main` `3a9af72`.

## Why all of these are the same underlying mistake

Every "% of gene peptides that need to pass" control in the app selects a **precomputed percentile
column**; it does not count peptides. `pep/filter.py` appends the suffix and the column is a
linearly interpolated quantile of the peptide-level metric:

| filter | column tested | operator | at the published settings |
|---|---|---|---|
| human identity | `blt_pident_p{human_id_gene_pc}` | `<= 80` | `_p75` |
| human alignment length | `blt_length_p{human_id_gene_pc}` | `<= 15` | `_p75` |
| strain conservation | `hap_top_hap_freq_p{100 − haplotype_gene_pc}` | `>= 0.99` | `_p05` |
| indel frequency | `max_indel_freq_type_any_fs_any_p{indel_gene_pc}` | `<= 0.05` | `_p100` |

So "at least N% of peptides must pass" is a *description* of a percentile test, not the test itself,
and the two come apart:

- The two human filters are separate marginal tests combined with `AND`. The app never asks whether
  the **same** 75% of peptides satisfies both. **10 of the 253 retained genes** have fewer than 75%
  of their peptides satisfying both conditions at once (`PF3D7_1203200`, 70.0%, is the extreme), and
  the identity-passing and length-passing peptide subsets differ for **4,846 of 4,937** genes.
- A percentile is not a count. **19 of the 253** are retained on fewer than 95% of their peptides
  meeting the conservation threshold — `PF3D7_0308800` gets in on 9 of 10 peptides, because its
  interpolated 5th percentile is 0.9914.
- The indel rule is the one exception: `_p100` is the maximum, where percentile and count agree
  exactly, so "all peptides within each gene" is literally true on all 4,937 genes.

---

## 1. Results — the 1,274-gene sentence *(in print; wrong by a factor of eight)*

**As published:**

> "Applying both the strain conservation and human identity filters simultaneously–i.e. by requiring
> 95% of peptides per gene to exceed a 99% conservation threshold and remain below 80% identity to
> any human exon with a maximum alignment length of 15 amino acids–reduced the number of retained
> genes to 1,274"

**Two problems.** "95% of peptides per gene" is attached to *both* filters, but the human gene rule
is **75%**, not 95%. And "simultaneously … and" reads as one joint subset of peptides. The three
readings give wildly different answers:

| reading | genes |
|---|---|
| human `pc=75` + conservation `pc=95` at 0.99 — **what was actually run** | **1,274** |
| human `pc=95` + conservation `pc=95` at 0.99 — what the sentence says | 158 |
| ≥95% of peptides passing conservation *and* human identity jointly | 92 |

Of the app's 1,274 genes, only 92 satisfy the joint reading; 1,182 do not.

**Replacement:**

> Applying the strain conservation and human identity filters together reduced the number of
> retained genes to 1,274. Each filter is evaluated as a percentile of its own peptide-level metric
> across the gene's peptides: genes were retained where the 5th percentile of peptide haplotype
> frequency was ≥99% (the percentile form of the application's "at least 95% of peptides" gene
> rule), the 75th percentile of peptide identity to human exons was ≤80%, and the 75th percentile of
> peptide alignment length was ≤15 amino acids. The three criteria are evaluated independently of
> one another, so they need not be satisfied by the same subset of a gene's peptides.

## 2. Results — the 2,107-gene sentence *(in print; parenthetical names the wrong threshold)*

**As published:**

> "When configuring Pf-PeptideFilter to retain only those genes for which at least 95% of peptides
> exhibited a dominant haplotype (i.e. a haplotype present in 95% of isolates or more), a total of
> 2,107 genes and 84,556 peptides were retained"

**Problem.** The parenthetical should say **99%**, not 95%. It has picked up the gene rule's 95% and
applied it to the frequency threshold as well. Verified:

| configuration | genes |
|---|---|
| gene rule 95% of peptides, threshold **0.99** — what gives the published figure | **2,107** |
| gene rule 95% of peptides, threshold 0.95 — what the parenthetical describes | 3,824 |

The counts in the sentence are correct: 2,107 genes / 84,556 peptides both reproduce exactly, as do
1,517 / 42,208 for the 99%-of-peptides row.

**Replacement:** change the parenthetical only —

> "… (i.e. a haplotype present in 99% of isolates or more) …"

## 3. Results — the default-thresholds parenthesis, `< 5%` vs `≤ 5%` *(in print; no numerical effect)*

**As published:** "indel frequency < 5%".

**Problem.** The code tests `<= 0.05`. Here it changes nothing — every gene at or below 0.05 sits at
exactly 0.0 and the next distinct value in the data is 0.057, so strict and non-strict both admit
3,457 genes. Worth aligning anyway, since a reviewer reproducing the work will read the operator.

**Replacement:** "indel frequency ≤ 5%".

## 4. Results — the homology filter names two species *(in print; the run used one)*

**As published:**

> "adding filters for indel frequency, liver-stage gene expression, and cross-species homology
> (presence of one-to-one orthologs in P. vivax, and P. berghei) further narrows this candidate set."

**Problem.** The published configuration requires ***P. vivax* only**. Requiring both gives **251
genes**, not 253 — and the very next sentence correctly says "homology with *P. vivax*", so the two
sentences contradict each other.

**Replacement:**

> adding filters for indel frequency, liver-stage gene expression, and cross-species homology
> (presence of a one-to-one orthologue in *P. vivax*) further narrows this candidate set.

If the intent was to describe the filter's *capability* rather than the run, say so explicitly: "the
homology filter can require orthologues in *P. berghei*, *P. vivax* or *P. knowlesi*; the
configuration used here required *P. vivax*."

## 5. Methods — "averaged across genes" *(in print; it is not a mean)*

**As published**, in Filtering Criteria:

> "gene expression is evaluated at the gene level, whereas human identity is assessed per peptide and
> averaged across genes."

and in Table 1: "Filtering is performed at the peptide level and summarised per gene."

**Problem.** Neither is a mean. The default gene rule is a **75th-percentile test**. "Averaged" is
wrong; "summarised" is too vague to reproduce from.

**Replacement:**

> gene expression is evaluated at the gene level, whereas human identity is computed per peptide and
> summarised to the gene as a percentile of its peptide values — by default the 75th percentile,
> equivalent to requiring that at least 75% of the gene's peptides pass.

And in Table 1: "Filtering is performed at the peptide level and summarised to the gene as a
user-selected percentile of the gene's peptide values (default: 75th percentile)."

## 6. Methods — the expanded configuration paragraph *(revision draft, not yet in print)*

This is the paragraph added for the revision. The corrected version is in
`sensitive-analysis-plan.md` (the original is preserved beneath it under "Superseded text") and in
§02 of `analysis/report.html`. Three clauses changed:

| clause | was | is |
|---|---|---|
| indel gene rule | "at least 95% of peptides" | "all peptides within each gene (100%)" — the app default. 95% gives 291 genes / 13,772 peptides, not 253 / 10,088 |
| human identity | "at least 75% of peptides … applied jointly" | two separate 75th-percentile tests, stated as such |
| conservation | "at least 95% of peptides … met the ≥99% threshold" | the 5th percentile of peptide haplotype frequency ≥99% |

One item in that paragraph is still open: the `[NEED TO ADD REF]` on the Day-4 justification. The
claim is now supported by the application's own expression table, which can be cited alongside any
external reference — genes at ≥1 CPM in all three replicates: **day 2 → 36, day 4 → 1,583, day 5 →
4,375, day 6 → 4,657** (of 4,937 genes; median CPM 0.00 / 22.09 / 64.77 / 67.45). Day 2 is unusable
at this rule and days 5 and 6 are near saturation (89% and 94% of all genes), so day 4 is the only
one of the four that is both well-powered and selective. An external reference is still preferable.

---

## One thing that is *not* a manuscript correction

The app's "joint vs independent" radio button (`human_id_joint_rule`) has **no effect on any result,
ever** — two column tests combined with `AND` in one step give the same set as the same two tests
applied sequentially. Both options return the identical 253 genes, and 1,142 under application
defaults. It changes only how the app groups the two tests in its own step summary. It is therefore
deliberately *not* mentioned in the corrected Methods paragraph: naming it would imply a constraint
that is not applied.

## Separately — a documentation bug in the deployed repository

Not a manuscript issue, but worth reporting to whoever owns the repo. The README **on `prod`** (the
branch the Streamlit app is deployed from, verified at `5861f79`) says at line 17: *"Please use
`requirements.txt` and Python 3.12"*. But `requirements.txt` pins `pandas==1.5.1`, which publishes no
cp312 wheel — PyPI has cp38, cp39, cp310 and cp311 only — so it cannot be installed on the Python
the README asks for without building from source. Anyone following the README's own instructions
fails at `pip install`. Either the pin or the instruction needs to move; 3.11 works today. Note that
`main`'s README is a four-line stub that mentions no Python version, so the bug is only visible on
`prod`.
