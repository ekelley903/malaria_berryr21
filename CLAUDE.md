# Malaria PepSeq Project Notes

## Batch effect investigation (IgA, Plate 2 vs Plate 4)

Source: Alexis's `~/Downloads/Analysis/20250219_explore_IgA_batch.html` (2025-02-19) and follow-up `20250224_top_z_peptides_follow_up.html` (2025-02-24).

- Raw IgA counts cluster into two groups that align almost perfectly with sequencing plate (Plate 2 vs Plate 4), not with any biological covariate tested (sex, location, timepoint, collection date). A random forest predicting plate from raw counts got 0.92% OOB error (~99% accuracy) — a strong batch effect.
- Correction methods tested, evaluated by whether a random forest could still predict plate afterward (lower error = batch signal still present):
  - Raw (baseline): 0.92% error
  - Z-scores: 3.67%
  - Whole-plate scaling: 1.83%
  - Quantile normalization: 2.75%
  - **ComBat: 0% error** (perfect confusion matrix — plate still fully separable)
- **Conclusion: ComBat did not resolve the batch effect.** By this diagnostic it was the worst of all methods tried — post-ComBat, plate was *more* separable than in the raw data or any simpler correction.
- For reference, a permutation-null test (shuffled plate labels) gives ~50% error, confirming that even the best correction attempts (z-score filtering: 7–23% error depending on threshold) still leave a detectable residual batch signal, just less extreme than ComBat's 0%.
- After this result, the analysis pivoted to filtering peptides by mean z-score threshold (>6 or >10) as an alternative to ComBat, which reduced but did not eliminate plate-predictability.

**Open question / to revisit:** no correction method fully eliminated the batch effect. If batch-corrected data is used downstream, treat residual plate signal as a known confound.

**Hypothesis for why ComBat did worst: it was run on raw counts, not normalized data.** Checked the code in `20250219_explore_IgA_batch.html` — the only `ComBat()` call in the whole `Analysis` folder is:
```r
data = as.data.frame(t(IgA_raw_sample_info_w_plates[,9:59491]))
combat_out = ComBat(data, batch = IgA_raw_sample_info_w_plates$plate)
```
This is straight off raw, untransformed counts (no log transform, no library-size/quantile normalization first). The z-score/whole-plate-scaling/quantile-normalization approaches in the same doc were tried as separate alternatives, never combined with ComBat.

Classic `ComBat` (sva package) fits a parametric empirical-Bayes location/scale model assuming roughly-Gaussian features — built for microarray intensities or already-normalized/log-transformed data. Raw PepSeq counts are non-negative, heavy-tailed, and have very different per-plate library-size distributions (shown in her own histograms), so feeding raw counts into a Gaussian L/S model confounds the batch correction with an uncorrected scale effect — plausibly explaining the 0% OOB error (correction distorting rather than removing the batch signal).

**Untested next steps:**
1. Log-transform or quantile-normalize the counts first, *then* run ComBat on the normalized values.
2. Try `ComBat-seq` (negative-binomial-based, designed for raw RNA-seq/count data) directly on raw counts instead of classic `ComBat`.
3. Pass a `mod=` model matrix with biological covariates (sex, location, timepoint) to `ComBat()`. The call above used no `mod` argument (defaults to `NULL`), so no biological variation was protected from removal alongside the batch effect — separate issue from the raw-counts problem, but another gap in the current run worth closing.

**ComBat vs ComBat-seq signatures** (both in the `sva` package):
```r
# classic ComBat — expects continuous/normalized data
ComBat(dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
       mean.only = FALSE, ref.batch = NULL)
# what she ran:
combat_out = ComBat(data, batch = IgA_raw_sample_info_w_plates$plate)

# ComBat-seq — expects raw integer counts directly, no transform needed
ComBat_seq(counts, batch, group = NULL, covar_mod = NULL,
           full_mod = TRUE, shrink = FALSE, shrink.disp = FALSE,
           gene.subset.n = NULL)
# untried alternative:
combat_out = ComBat_seq(counts = as.matrix(data),
                         batch = IgA_raw_sample_info_w_plates$plate,
                         group = IgA_raw_sample_info_w_plates$sex)  # or other biological var to protect
```
Key differences: `ComBat` takes one combined `mod` matrix and returns continuous (possibly negative) adjusted values; `ComBat_seq` splits biological-covariate protection into `group` (main variable of interest) + `covar_mod` (other covariates), takes raw non-negative counts, and returns adjusted integer counts. Using `ComBat_seq` would resolve both open gaps at once: raw-count input (no Gaussian assumption) and covariate protection.

## Follow-up (2026-07-15): ComBat-seq tried — root cause found, not fixable by correction alone

Ran `ComBat_seq` in `2-exploratory/UMAPs/2-20260714_umap.qmd`, closing both gaps flagged above: raw integer counts as input (not Z-scores or norm counts), and `group = cohort` passed to protect Mfera vs Liwonde biological variation from being removed alongside the batch effect. Restricted to the same `>=3`-samples-Z>6-filtered peptide set used throughout that script, for both IgG and IgA.

**Result: ComBat-seq did not remove the IgA plate separation visible in UMAP space**, despite correctly identifying the batches (2 for IgA, 3 for IgG — plate counts are now accurate after fixing a separate well-column/plate mislabeling bug, see `0-sample_tidy` history) and applying covariate protection.

**Root cause identified: plate assignment is confounded with subject ID, and that confound is not explained by any covariate currently in the clinical metadata.** Within each cohort, the two IgA plates cover completely non-overlapping, contiguous ID ranges:

| Cohort | Plate 2 | Plate 4 |
|---|---|---|
| Mfera (`mfera_id`) | 20–71 (n=32) | 75–102 (n=29) |
| Liwonde (`subject_id`) | 30–202 (n=22) | 211–528 (n=26) |

`cohort` and `sex` are both reasonably balanced across the two plates (so protecting `cohort` in ComBat-seq wasn't wrong, it just wasn't the actual confound). Checked whether the ID split tracks a date-based confound — **the answer differs by cohort**:
- **Liwonde**: `subject_id` correlates with collection date (Pearson r=0.71, p=1.6e-8). Plate 2 spans 2017-12-20–2018-11-14; Plate 4 spans 2018-02-06–2019-04-24 — overlapping, but Plate 4 skews markedly later (5+ months beyond Plate 2's range). So part of Liwonde's plate confound is a real, if imperfect, calendar-time effect.
- **Mfera**: `mfera_id` does *not* correlate with `infectiondate` (r=0.12, p=0.37, not significant) — the two plates' infection-date ranges are nearly identical (Plate 2: 2017-03-27–2018-06-18; Plate 4: 2017-04-26–2018-07-10). The same ID-contiguous-block-per-plate pattern shows up here too, but it isn't explained by date. Mechanism unverified — could be wet-lab logistics (tubes pulled from storage/racks in ID order) unrelated to any biological variable, or something not captured in this dataset. Do not assume it's the same "enrollment order" story as Liwonde.

**Why no correction method can fully fix this**: ComBat/ComBat-seq/Harmony/RUV-style corrections can only protect biological signal that's given to them as an explicit covariate (`group`/`mod`). Here, the actual confound isn't captured by any single covariate on hand — not fully explained by cohort, sex, or (for Mfera) even collection date — so there's no `group=`/`mod=` argument that could have protected against it. This is a design/data-availability gap, not something fixable by picking a different correction method. IgA is additionally the most exposed to this problem: unlike IgG (which spans 2 sequencing runs, IM0163 + IM0174, and 3 plates), all 109 real IgA samples come from a single run (IM0163) split across exactly 2 plates — there's no second run to provide independent replication or leverage, and no bridge/reference samples repeated across plates to empirically check whether any correction preserved biology.

**Practical implications flagged for downstream modeling** (raised by the ML collaborator): any model trained on IgA data risks partly learning whatever subject-ID-correlated signal plate carries (confirmed real for Liwonde via date; unverified but present for Mfera), rather than purely the intended biological signal. Standard random k-fold CV would not catch this (train/test folds would share the same plate-label correlation) and would report misleadingly good performance that might not generalize to a new run/plate.

**Recommendations (not yet implemented):**
- Use plate/run-blocked cross-validation for any IgA modeling (never split within a plate) so leakage becomes visible instead of hidden.
- Include `plate_id` as an explicit covariate/random effect in statistical models rather than trying to fully remove it upstream — keeps the confound structure honest.
- Prefer within-subject, within-plate comparisons (e.g. paired serum1/serum2) over between-subject, cross-plate comparisons where possible.
- Treat IgA findings as provisional until validated on independently (ideally randomly) plated samples.
- For future sequencing runs: randomize sample-to-plate assignment with respect to subject ID, and include bridge/reference samples replicated across every plate — without them there's no ground truth to confirm a correction method preserved biology rather than just reducing plate-predictability.

**Good news for the central pre/post pair analysis**: checked whether each subject's serum1/serum2 (Liwonde) or serum1e/serum2e (Mfera) time-pair lands on the same plate or is split across plates — this matters because the paired pre/post comparison is the primary analysis, not the between-subject/cross-plate comparisons the confound above threatens most. Result: **100% same-plate, zero exceptions**, across every pair in both isotypes and both cohorts (IgA: 21/21 Liwonde, 27/27 Mfera; IgG: 23/23 Liwonde, 29/29 Mfera — 100 pairs total). Every subject's two timepoints were plated together.

This means a per-subject, shared plate shift mostly cancels out in any within-subject paired comparison (e.g. serum2 − serum1, or a paired fold-change) — the batch confound is a much smaller threat to the core pre/post analysis than the raw UMAP clustering makes it look. It does **not** protect against: plate-specific technical noise that isn't shared between a subject's two timepoints (e.g. well-position or degradation effects that don't cancel within a pair), or any analysis that compares *across* subjects/plates (e.g. cohort-level fold-change comparisons, or any classifier trained across subjects) — the leakage/confounding risks described above still fully apply there.
