# Draft covering note — edit freely before sending

**Subject:** R2.1 / R2.3 analysis done — six wording corrections needed, five of them already in print

---

Hi Jacob,

I've run the combined analysis. Full write-up here: **[LINK]** (PDF also attached). Short version
below, and there are three things I'd like you to confirm.

**The plan works, but four of the ten sensitivity scenarios had to be replaced.** Three of them changed
nothing whatsoever — the identical 253 genes came back — and one emptied the catalogue entirely:

- Human identity at ≤90% and ≤70% both return exactly the published gene set. Under the joint
  rule the ≤15 aa alignment-length criterion is doing all the work; identity is inert at every
  value from 60% to 100%. So I varied alignment length instead (≤12 / ≤15 / ≤18 → 183 / 253 /
  376 genes), and added an arm using the app's "deactivate length filter" option, where identity
  does bite (≤70 / ≤80 / ≤90 → 506 / 558 / 587).
- Indel frequency at ≤0.01 also returns the published set. The gene-level metric is a step
  function: every gene at or below 0.05 sits at exactly 0.00 and the next distinct value in the
  data is 0.057, so **3,457 genes** clear the indel filter at any cutoff from 0.00 to 0.05. (An
  isolated run through the app reports 3,451 rather than 3,457, because the app has no way to
  switch a filter off — the other filters are held at permissive thresholds, and six genes whose
  *unrelated* metrics are NaN fail those anyway. None of the six is in the 253. 3,457 is the
  number to quote for the indel filter itself.) There is no stringent setting available on *any*
  of the filter's four axes, which I checked rather than assumed: the frequency threshold is inert
  from 0.05 down to 0.00, the gene rule is already at its 100% maximum, and restricting which
  indels are counted *relaxes* it — "Only Insertions" 280 genes, "Only Deletions" 260, "Only
  Frameshifts" 294, insertions-and-frameshifts 297, against the baseline 253. Counting fewer
  indels lowers each gene's maximum, so more genes pass. The filter is already at maximum
  stringency, so that arm is relaxation-only and I've said so explicitly.
- Conservation at 1.00 gives **zero genes**. Only 19 of 4,937 reach it and none survive the other
  filters. The stringent arm now tightens the *gene rule* instead — `haplotype_gene_pc = 99`,
  **176 genes** — which keeps the threshold at the published 0.99 and stays reproducible in the
  app. The off-grid `strain_conservation = 0.995` option gives 164 genes, but is not selectable
  on the app's 0.01 slider; it is kept in `rejected_scenarios.csv` with the reason.

Everything else behaved, and the consistency checks in your Step 11 hold throughout — verified by
set containment, not by comparing counts: M253 is a subset of every relaxed catalogue and a superset
of every tightened one, at both gene and peptide level, across all fifteen runs. One label needs
reading carefully: HUMAN-IDONLY-STRINGENT is named for its identity threshold (≤70%) but is a
relaxation overall, because switching the length filter off removes a constraint — it gains 253
genes and loses none.

**The headline for R2.3:** the catalogue is most sensitive to expression and conservation —
relaxing the Day-4 replicate rule more than doubles it (596 genes). Indels and orthology depend on
which axis you move: the indel *frequency* threshold barely matters (+7 genes) and tightening
orthology costs 2, but relaxing the indel *gene rule* to 95% adds 38 (+15%) and dropping the
orthology requirement adds 37 (+15%). I'd avoid calling those two "a few percent".

**Six wording corrections before the revision goes out — five already in print.** Every *count* in
the preprint reproduces exactly; it is several of the printed *descriptions* of how those counts
were obtained that need correcting, along with the
expanded methods paragraph that needs work. The numbers you quoted for the human (75%) and
conservation (95%) gene rules are right, but two of the three clauses describe the rule in a way
the code doesn't support, and one number is wrong:

1. **The indel gene rule was 100%, not 95%** — the app default. At 95% you get 291 genes /
   13,772 peptides rather than 253 / 10,088.
2. **"with the identity and alignment-length criteria applied jointly" is not what happens.** The
   code tests two separate percentile columns and ANDs them: 75th-percentile identity ≤80% *and*
   75th-percentile alignment length ≤15 aa. It never asks whether the same 75% of peptides meets
   both. That isn't pedantry — 10 of our 253 genes have fewer than 75% of their peptides passing
   both conditions at once, so a filter written the way the paragraph reads would drop them.
   (Related: the app's "joint vs independent" radio button makes no difference to any result at
   all — two column tests ANDed in one step or applied in two give the same set. It only changes
   how the app groups them in its own step summary, so I've left it out of the Methods rather
   than imply a constraint that isn't applied.)
3. **The conservation clause has the same shape of problem.** "At least 95% of peptides … met the
   ≥99% threshold" is really `hap_top_hap_freq_p05 ≥ 0.99` — a linearly interpolated 5th
   percentile, not a count. 19 of the 253 are retained on fewer than 95% of peptides passing;
   PF3D7_0308800 gets in on 9 of 10. The indel clause needs no repair, because its rule is the
   100th percentile (the maximum), where percentile and count agree exactly.

There's a redrafted paragraph in §02 of the write-up that states all three as computed. It also
makes explicit that 10,088 counts the peptides *belonging to* the 253 genes — the app's
independently filtered peptide count for the same settings is 84,859, and it's worth the paper
saying which it means.

**And I'm afraid five of these are already in print.** Every *count* in the preprint reproduces
exactly — I've now got all six of them asserted before any analysis runs. But two sentences
describing how the counts were obtained don't match the configuration that produces them:

- **The 1,274 sentence** says we applied both filters "by requiring 95% of peptides per gene to
  exceed a 99% conservation threshold and remain below 80% identity…". The human gene rule is
  **75%**, not 95%. Run as that sentence describes it gives **158 genes**. And "simultaneously …
  and" reads as one joint subset of peptides, which gives **92**. So it's 1,274 / 158 / 92
  depending on how you read it, and only the first is what we did.
- **The 2,107 sentence** says the dominant haplotype had to be "present in 95% of isolates or
  more". It should be **99%** — the parenthetical has picked up the gene rule's 95% and applied
  it to the frequency threshold. At a 0.95 threshold the same rule gives **3,824 genes**.
  The counts in that sentence are fine: 2,107 / 84,556 and 1,517 / 42,208 all reproduce.

Two smaller mismatches in the same passage, neither of which changes a number, are also already in print: the paper says
"indel frequency < 5%" where the code tests ≤ 5% (identical here, since everything at or below
0.05 is exactly 0.00), and it describes the homology filter as "one-to-one orthologs in *P. vivax*,
and *P. berghei*" a sentence before saying "homology with *P. vivax*" — the published run is
*P. vivax* only, and requiring both gives 251 genes rather than 253.

I'd suggest we correct all of these in the revision alongside the methods paragraph, and say so
plainly in the response to R2 rather than hoping it isn't noticed — R2.3 is a reproducibility
question, and these are exactly the sentences someone checking reproducibility would hit first. The
results are unaffected; it's the descriptions that need fixing.

**`analysis/MANUSCRIPT-CORRECTIONS.md` has all six of them copy-ready** — sentence as published,
replacement underneath, with the numbers behind each. That's the file to work from rather than the
write-up. I've also replaced the paragraph in `sensitive-analysis-plan.md` with the corrected
version and kept the original beneath it marked "Superseded text", so you can diff it against
whatever is currently in the manuscript.

**On R2.1 — all eight antigens are excluded, in every single scenario.** One wording point for the
response: the reviewer asked about "licensed and leading" antigens, and only **CSP** underlies a
licensed product — both WHO-recommended and prequalified vaccines, RTS,S/AS01 and R21/Matrix-M,
are CSP-based. The other seven are candidates at various stages. Worth saying "one licensed and
seven leading candidates" rather than implying all eight are licensed. So the retention column
reads `none retained` fifteen times over and the retained/excluded matrix (8 × 15, since we ended
up with fifteen scenarios rather than the twelve in the plan) is one flat colour. Neither is worth
publishing as designed.

But I think the result is stronger than it first looks, not weaker. Every one of the eight fails
on **conservation**, and the Discussion already argues exactly this — that CSP, AMA1 and MSP1
underperformed in trials because of antigenic variation. The benchmark puts numbers on a claim
the paper already makes: AMA1 clears the conservation criterion at only the 25% gene-rule setting
against a 95% requirement, and MSP1 at 50%. CSP is the one to state carefully — it misses badly on
the threshold it tolerates (0.34 against 0.99) but **74.4%** of its peptides clear the 0.99
threshold, which is second only to RH5's 80.8%. (The 70 in the table is the gene-rule slider
setting, not a peptide fraction — both columns are in Table B and they are not the same thing.) The two measures rank the eight differently, so I
wouldn't quote either alone. One further caveat on all of these: they are settings of the app's own
gene-rule slider, not exact peptide fractions — the rule tests an interpolated percentile, so a
value can sit a few points either side of the true share (CelTOS reads 40 where 50.0% of its
peptides actually pass). So the answer to "does your tool find the known
antigens?" is: no, and it rejects them for precisely the reason the paper identifies as having
limited those vaccines. That reads as the conservation-first design working, and I'd frame it
that way rather than defensively.

I've replaced the heatmap with a distance-to-threshold figure, which shows how *far* each antigen
is from each criterion — that separates AMA1 (hopeless, 25%) from RH5 (marginal, 80%), which a
pass/fail grid throws away.

**Three things to confirm:**

1. **Human filter** — I've included both arms (alignment length, and identity-with-length-off).
   *I'd keep both*: five rows rather than two, but together they make the point that identity is
   non-limiting under the default rule, which is worth a sentence in the Results. Happy to drop
   the identity-only arm if you'd rather keep the table tight.
2. **Indels** — *I'd keep it relaxation-only* and state plainly that the baseline already sits at
   maximum stringency. Inventing a "stringent" setting that silently equals the baseline would be
   worse than admitting the limit.
3. **Conservation stringent** — now `haplotype_gene_pc = 99`, **176 genes / 5,007 peptides**,
   a strict subset of the 253 (0 gained, 77 lost, Jaccard 0.696). This replaces
   `strain_conservation = 0.995` (164 genes), which was real but not selectable on the app's 0.01
   slider — and there was no on-grid threshold alternative, since the next step above 0.99 is 1.00
   and that empties the catalogue. Tightening the gene rule instead has a distinct advantage: the
   paper *already* reports a 99% conservation gene rule, because that is how the published
   1,517-gene figure is produced, so the Methods can cite a value the analysis already uses rather
   than introducing a new one. Table A, the retention matrix, the frozen configs and the per-run
   outputs have all been regenerated.

**A scientific pass as well as a technical one.** Beyond checking that the numbers reproduce,
I looked at whether the comparisons are fair and the inferences warranted. Five things came out of
that, all written up in §A of the report:

- **Conservation is the only filter all eight antigens fail — though it is not the only one that discriminates them.** Drop conservation and keep human similarity and indels: the app-style run returns 2,397 genes, and only three of the eight re-enter (AMA1, RH5, CelTOS); CSP is excluded by human similarity and four more by indels, both stage-independent. Six of the eight fail the
  liver-stage day-4 expression filter, but seven of the eight are expressed in sporozoites, and
  substituting the sporozoite criterion retains none of them — nor does removing expression and
  orthology altogether. RH5 having no *P. vivax* orthologue is expected biology for a
  *P. falciparum*-restricted gene, not a discriminating test. So "all eight fail conservation" is
  the defensible headline; the other per-filter failures are largely stage-explained and we
  shouldn't present them as five independent failures.
- **The benchmark was close to circular, and percentile ranks fix it.** As pass/fail it only
  confirms the filter does what it is defined to do. Ranked against all 4,937 genes it becomes a
  real result: five of the eight sit in the **bottom 1% of the genome for conservation** (AMA1
  0.02%, CelTOS 0.04%, TRAP 0.14%, CSP 0.30%, MSP1 0.83%), with EBA175 the sixth just outside at
  2.03%, while all of them rank mid-range on human identity. That is the number I'd quote in the Results. New columns in Table B.
- **Conservation is measured in African isolates only** — 8,492 of Pf7's 16,203 genomes. The paper
  states the scope but never flags it as a limitation, and it is one for a vaccine-target
  shortlist. The app ships only the African metric, so no setting of any filter can test
  geographic robustness; that's worth a sentence in the Discussion.
- **"Day 2 is unusable" is a library-depth artefact, not parasite biology.** One day-2 replicate
  totals 15,682 CPM against ~1,000,000 at days 5–6, so there is almost no parasite material in it,
  and 3,367 of 4,937 genes are zero across all three day-2 replicates. Day 4 is still the right
  choice, but for that reason. It also means ≥1 CPM is not a constant stringency across stages.
- **The 0.99 threshold can't be checked against sampling error** — neither table carries a coverage
  or callable-sample count, so there is no denominator. What can be said is that 18 of the 253 sit
  within 0.001 of the cut and 89 within 0.005, so the catalogue is genuinely sensitive to sub-1%
  movements in that threshold.

One more of the same kind, in the published Methods rather than the Results. They say human
identity is "assessed per peptide and averaged across genes", and Table 1 says "filtering is
performed at the peptide level and summarised per gene". Neither is a mean — it's a
75th-percentile test. "Averaged" is wrong and "summarised" is too vague to reproduce from; both
should name the percentile.

And one documentation bug in the repo itself, which someone should pass on to whoever owns it. The
README **on `prod`** — the branch the Streamlit app is deployed from — says "Please use
`requirements.txt` and Python 3.12", but `requirements.txt` pins `pandas==1.5.1`, which publishes no
cp312 wheel. So anyone following the deployed README's own instructions fails at `pip install`.
Either the pin or the instruction needs to move; 3.11 works fine. (Easy one to miss: `main`'s README
is a four-line stub with no Python version in it, so the bug only shows on `prod`.)

Everything is re-runnable: scripts in `analysis/`, outputs in `analysis/output/`, including a
frozen config file and the gene- and peptide-level catalogue for each of the fifteen runs, plus the
source of this write-up (`analysis/report.html`) so the PDF can be regenerated. One caveat to be
straight about: `analysis/` is untracked, so it isn't in the public repo — someone cloning
malariagen/pf-peptidefilter gets the app and the data but not these scripts. If the revision cites
the analysis, it needs to land somewhere citable. All of
it is driven through the app's own `filter_datasets()` against the live repo's data — the four
numbers already in the preprint (2,107 / 1,517 / 1,274 / 253+10,088) reproduce exactly, and that
check runs as an assertion before every analysis. I've also re-derived 253 / 10,088 and all three
of the other published figures straight from `data/*.csv.gz` without calling the app at all, as an
independent cross-check: the gene sets are identical, not just the counts. Environment was Python
3.10.12 with pandas 1.5.1 and numpy 1.26.4, and the provenance check was re-run against live
`main` 3a9af72 / `prod` 5861f79 (unmoved). The PDF renderer is pinned separately in
`analysis/requirements-report.txt` (WeasyPrint 70.0) rather than in the repo's `requirements.txt`,
which is the app's dependency set and is byte-compared against live `main` — worth knowing because
pagination moves between WeasyPrint versions, so don't cite page numbers of the PDF.

Andrew
