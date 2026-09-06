# Brain Network Group Analysis & Visualization

MATLAB code for the statistics and figure-generation stage of the Stochastic Actor-Oriented Models (SAOM) of brain networks analysis pipeline. Here comparing network metrics across clinical/demographic groups, modeling their relationship with age, and
producing publication-style figures at both the individual-metric and
whole-battery scale.

## What this does

Given results of the SAOM for functional (fMRI) and
structural (DTI) time windows with particular factors and factor weights:

- **Runs weighted statistical comparisons** between groups (diagnosis,
  sex, their interaction, and multi-level severity), using inverse-variance
  weighting so subjects with noisier per-metric estimates contribute less,
  via weighted t-tests, one-way and two-way ANOVAs.
- **Fits age-relationship models** (linear and optional quadratic, with a
  mixed-effects option) per metric per group, and extracts group-specific
  age slopes and their significance.
- **Generates figures at three scales**:
  - *Individual*: one violin/box plot or age-trend scatter per
    metric x comparison, with significance brackets and a summary
    effects box.
  - *Grid*: every metric for one comparison stacked into a single compact
    figure, so a whole battery of tests can be scanned at once.
  - *Aggregate*: Functional vs. Structural side-by-side, and
    all-subnetworks summaries, built from the same per-metric results.
  - An **age-effect heatmap** (metric x term, colored and starred by
    effect size/significance) at the whole-brain, combined, and
    per-subnetwork-aggregate levels.
- **Exports everything twice**: a normal figure, and (for tests with a
  summary significance box) a companion "SigBoxOnly" version with the
  data hidden and just the significance annotation left — sized to overlay
  onto the normal version, for cases where the box's default placement
  needs manual adjustment. Every grid/heatmap/aggregate figure also gets a
  companion `.txt` file listing the exact p-value (and effect size, for the
  heatmaps) behind every bracket/cell in that image.
- **Logs results to CSV** (`MuSigmaResults.csv`, `EffectsResults.csv`,
  `SampleSizeSummary.csv`) alongside the figures, so downstream reporting
  doesn't require re-deriving anything from the plots.

Layout, label placement, and figure sizing are computed from measured text
extents and data ranges rather than hardcoded, so the same code adapts to
however many groups, metrics, or subnetworks a given comparison involves.

## Data availability

**This repository does not include any participant data.** The dataset
this code was written for involves human-subjects neuroimaging and
clinical data and is not publicly shareable. The script expects a handful
of preprocessed structures (`AW`, `As`, `Es`, `Efields` — see the comment
at the top of `ThresholdGraphs1.m`) to already exist in the MATLAB
workspace; the data-preparation stage that builds those from raw
connectome data is a separate, unincluded step. This repository is meant
to demonstrate the statistical-analysis and figure-generation approach,
not to be run end-to-end without access to that data.

`Graphing.mat` *is* included — it's just color palettes, marker symbols,
and group-code label arrays used for consistent plot styling, with no
subject-level information in it.

## Dependencies

- MATLAB (developed against a recent release; uses `string` arrays,
  `arrayfun`/`cellfun`, `annotation`, mixed-effects modeling via
  `fitlme`/`fitglme`, and other modern-ish MATLAB features).
- [`violin.m`](https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot-using-matlab-default-kernel-density-estimation)
  by Holger Hoffmann (2015), included in this repo. Not written by me —
  see the citation notice at the top of that file.

## Structure

`ThresholdGraphs1.m` is a single script file organized into three phases
(data cleaning → statistical tests → graphing), followed by a local-functions
section containing every plotting/statistics helper. It's long because it
covers a lot of different comparisons and plot types with a consistent,
reusable set of helpers rather than duplicating logic per comparison.

## About this project

This code was written over the course of my own graduate research,
iterating on both the statistics and the figure design as the analysis
needs evolved. I'm not a software engineer by training — this repository
is here to show how I approach a real, messy, iterative data-analysis
problem, not to present idealized textbook code.

## License

MIT — see [LICENSE](LICENSE). Note that `violin.m` is third-party code
(see its own header for citation); the license here applies to my own code
in `ThresholdGraphs1.m`.
