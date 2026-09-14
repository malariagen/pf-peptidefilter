# Column definitions — supplementary tables A and B

Written because three columns in Table B are easy to misread. Everything here was produced by the
published app's own `filter_datasets()` (`pep/filter.py`), driven through
`analysis/pepfilter_headless.py`.

## The one thing to know about Table B

`human_gene_rule_pct_tolerated`, `conservation_gene_rule_pct_tolerated` and
`indel_gene_rule_pct_tolerated` are **settings on the app's "% of gene peptides that need to pass"
control** — the most demanding of its 39 options at which the gene still survives that filter. They
are **not** the percentage of the gene's peptides that pass.

These columns were previously named `*_pct_peptides_passing`, which asserted exactly the wrong
thing and has been corrected.

## The quantity the rename removed, now restored

Renaming those columns was right, but it left the table without the thing plan Step 14 actually
asked for: "the relevant **percentage of peptides passing** rather than only a binary PASS/FAIL
result … an antigen in which 94% of peptides pass a 95% requirement; from one in which only 20%
pass."

So Table B now carries **both**, and they answer different questions:

| column | what it is |
|---|---|
| `*_gene_rule_pct_tolerated` | the most demanding of the app's 39 gene-rule settings at which the gene still survives that filter. A **slider setting**. |
| `*_true_pct_peptides_passing` | the actual proportion of that gene's peptides satisfying the peptide-level condition, computed from the peptide table. A **count**. |

They differ in 4 of the 24 antigen × filter cells, in both directions, because the gene rule tests
an interpolated percentile rather than counting. Neither is wrong; quote the true fraction when the
question is "how far from the requirement is this gene", and the slider setting when the question
is "what could a user select and still keep it".

## Genome percentile ranks

`conservation_genome_percentile`, `human_identity_genome_percentile` and `indel_genome_percentile`
give each gene's rank among all 4,937 genes on that metric, stated in one direction throughout:
**"better than X% of all genes"**, so a *low* number is a *bad* gene.

Directionality is handled per metric — higher haplotype frequency is better; lower identity,
alignment length and indel frequency are better — and the column counts strictly-worse genes rather
than using a rank function, so ties stay out of the numerator and the direction is explicit. An
inverted comparison here would report 99.7 for a gene that is in fact worse than 99.7% of the
genome, which is exactly the error this phrasing is meant to prevent.

These exist because pass/fail against a single threshold cannot say whether an antigen is a
borderline miss or a genomic outlier. Six of the eight benchmark antigens sit in the bottom 2% of
all genes for conservation while ranking mid-range on human identity — a graded result that a
uniform block of "excluded" cannot express.

Both additions are descriptive reads of the input tables, not filter logic: no gene is retained or
dropped on the strength of either, and every PASS/FAIL in Table B still comes from
`filter_datasets()`.

The two differ because the control does not count peptides. It selects a precomputed
**linear-interpolation percentile column** (`pep/filter.py:287-315`), so the setting can fall
either side of the grid floor of the true peptide fraction. It does in 4 of 24 cells:

| cell | setting reported | true peptide fraction | grid floor of the true value |
|---|---|---|---|
| CelTOS conservation | 40 | 50.00% | 50 |
| EBA175 human | 85 | 83.33% | 82 |
| TRAP indel | 93 | 92.98% | 92 |
| MSP1 indel | 92 | 91.86% | 91 |

No PASS/FAIL verdict is affected. But do not quote these values as percentages of peptides.

## Table B

| column | meaning |
|---|---|
| `antigen`, `gene_id` | benchmark antigen and its PlasmoDB identifier |
| `gene_name_in_app` | the `gene_name` value as it appears in the app's gene table |
| `n_peptides` | peptides the gene contributes to the peptide table |
| `*_gene_rule_pct_tolerated` | **slider setting, not a peptide fraction** — see above |
| `*_requirement_pct` | the setting the published configuration requires (75 human, 95 conservation, 100 indel) |
| `human`, `conservation`, `indel`, `expression_D4`, `orthology_Pvivax` | PASS / FAIL for that filter, each from an app run with every other filter neutralised |
| `conservation_max_threshold_tolerated` | most stringent haplotype frequency the gene still meets at the 95%-of-peptides rule, swept at the app's own 0.01 step |
| `conservation_requirement_threshold` | the published threshold, 0.99 |
| `d4_replicates_expressed` / `_required` | liver-stage day-4 replicates at ≥1 CPM, and the 3 required |
| `ortholog_P_vivax` / `_berghei` / `_knowlesi` | orthologue present in that species |
| `APP_status`, `M253_status` | retained / excluded under application defaults and under the published catalogue |
| `n_filters_failed` | how many of the five filters the gene fails |
| `M253_exclusion_reason` | **every** filter the gene fails, not only the first one the sequential pipeline removes it at |

"Neutralised" is not "off": the app cannot switch a filter off, so the others are held at permissive
thresholds. `NaN <= x` is False in pandas, so 7 of 4,937 genes with a missing metric drop out
regardless of threshold. None is a benchmark antigen or a member of the 253.

## Table A

One row per scenario. `genes` / `peptides` are the retained genes and the peptides belonging to
them — the same definition the paper uses for 10,088, not the independently filtered peptide set.

For each of the two reference catalogues — `M253_*` (the published 253-gene catalogue) and `APP_*`
(application defaults) — `_shared`, `_pct_retained` and `_jaccard` are given at both gene and
peptide level, plus `genes_gained/lost_vs_M253` and `peptides_gained/lost_vs_M253`.

`benchmark_antigens_retained` reads `none retained` in all 15 rows: no benchmark antigen is retained under
any scenario. The column is kept because it was explicitly requested, and an empty result is a
result.

`exact_setting` is the single setting that differs from the published baseline, except for the three
`HUMAN-IDONLY-*` rows, which also change the identity/length combination mode.
