# Partition intra- and inter-observer uncertainty (replicate-based)

`observer_error_partitioning()` quantifies how much variance in trackway
parameters is due to observer and replicate effects
(inter-/intra-observer) versus genuine between-track (biological)
differences, using replicated digitizations and mixed-effects models (no
anatomical simulation).

## Usage

``` r
observer_error_partitioning(
  data,
  metadata,
  variables = c("TurnAng", "sdTurnAng", "Distance", "Length", "StLength", "sdStLength",
    "Sinuosity", "Straightness", "TrackWidth", "PaceAng")
)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- metadata:

  A data frame with one row per track (same order as
  `track_param(data)`). Required columns:

  - `replica` — replicate index within (`track`, `observer`).

  - `track` — biological track ID.

  - `observer` — observer ID.

- variables:

  A character vector of parameters to analyze. Default is
  `c("TurnAng","sdTurnAng","Distance","Length","StLength","sdStLength", "Sinuosity","Straightness","TrackWidth","PaceAng")`.

## Value

An `"error_partitioning"` R object consisting of a list containing the
following elements:

- summary:

  Data frame of variance components per variable. Columns: `variable`,
  `component`, `variance`, `percent`, and \\R^2_c\\. Components reflect
  what was estimable (e.g., `track`, `observer`, `observer:track`,
  `Residual`).

- snr:

  Data frame with the signal-to-noise ratio per variable. Columns:
  `variable`, `bio_component` (always `"track"`), `bio_var`,
  `error_var`, and `SNR`. The SNR is \$\$SNR = Var(track) /
  \[Var(observer) + Var(observer:track) + Var(Residual)\].\$\$

- qc:

  Compact quality-control table per variable with `SNR`, a qualitative
  `SNR_rating` (`"weak"`, `"moderate"`, `"strong"`), the `top_component`
  by percent, its `top_percent`, \\R^2_c\\, and whether the observer
  term was estimable (`observer_estimable`).

- analysis_table:

  Data frame used in the analysis: selected variables joined with
  `metadata` (columns typically include `replica`, `track`, and
  `observer`).

- models:

  Fitted mixed-effects models used to extract variance components. For
  linear variables, each entry stores the `lmer` fit and its variance
  components. For angular variables, each entry is a list with separate
  `sin` and `cos` fits.

- formulae:

  Model formula(e) per variable. For angular variables, a list
  containing the `sin` and `cos` formulas.

- track_names:

  Character vector of internal track labels used in the analysis.

## Details

This function partitions the variability of each selected metric into
components attributable to **biology** (consistent differences among
tracks) and to the **human sampling process**: differences between
observers (**inter-observer sampling error**) and differences between
repeated digitizations by the same observer (**intra-observer sampling
error**). Use it to audit reproducibility across people and sessions,
identify metrics that remain stable despite who digitizes them, and
prioritize where protocol changes (clearer landmark definitions,
training, calibration) will most reduce measurement noise. Metrics whose
variability is dominated by biological signal are more interpretable
downstream (testing, clustering, classification), whereas those
dominated by observer/replica effects warrant caution or methodological
refinement.

For each metric the model fits random intercepts for `track`,
`observer`, and their `observer:track` interaction.
Replicate-to-replicate scatter within the same (`observer`, `track`)
cell contributes to the residual (this is the intra-observer sampling
error). Random terms that appear with a single level are dropped
automatically.

\$\$y \sim 1 + (1\|track) + (1\|observer) + (1\|observer:track).\$\$

Angular variables (`TurnAng`, `PaceAng`) are handled in sine–cosine
space to respect circularity and avoid 0°/360° discontinuities: the same
random-effects structure is fit separately to \\\sin(\theta)\\ and
\\\cos(\theta)\\, and the two variance decompositions are averaged
(Fisher, 1995).

Robustness is summarized by the signal-to-noise ratio \$\$SNR =
\frac{Var(track)}{Var(observer) + Var(observer:track) +
Var(Residual)}.\$\$ As a rule of thumb:

- **SNR \< 1** — observer/replicate error exceeds biology; the metric is
  *weak*;

- **SNR ~ 1–2** — biology and error are comparable; *moderate*
  robustness;

- **SNR \> 2** — biology dominates; *strongly robust*.

A compact QC table reports SNR, the qualitative rating, the top variance
component by percent, and the conditional \\R^2_c\\; it also flags
whether the inter-observer term was estimable (\>1 level).

Practical designs:

- *Full inter + intra*: provide multiple observers *and* multiple
  replicates per observer. Both `observer` and `observer:track` are
  estimable; replicate scatter is residual.

- *Intra-only*: set `metadata$observer` to the same value for all rows
  and provide `replica` with \>1 levels. Only `track` remains; replicate
  scatter is residual.

- *Inter-only*: set `metadata$replica` to a single level and supply
  multiple observers. `observer` is estimable; the residual still
  represents replicate-level noise (if present).

Note that estimates can be sensitive when designs are highly unbalanced
or sample sizes are small; singular fits (random variances at/near zero)
may occur (REML) and should be interpreted cautiously. This function
does not simulate anatomical landmark jitter; to quantify sensitivity to
landmark placement itself (e.g., preservation, landmark ambiguity), use
[`anatomical_error_partitioning()`](https://macrofunuv.github.io/QuAnTeTrack/reference/anatomical_error_partitioning.md).

All **replicated digitizations** (across observers and/or replicate
attempts) must be bundled together in this single `data` object (i.e.,
in the same file), with each replicate represented as its own track
entry. Do not split replicates across multiple files/objects, as
`track_param(data)` and `metadata` must align one row per replicate in
the same order.

## Logo

![](figures/Logo.png)

## References

Fisher, N. I. (1995). Statistical Analysis of Circular Data. Cambridge
University Press

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`track_param`](https://macrofunuv.github.io/QuAnTeTrack/reference/track_param.md)

## Author

Humberto G. Ferrón

humberto.ferron@uv.es

Macroevolution and Functional Morphology Research Group
(www.macrofun.es)

Cavanilles Institute of Biodiversity and Evolutionary Biology

Calle Catedrático José Beltrán Martínez, nº 2

46980 Paterna - Valencia - Spain

Phone: +34 (9635) 44477

## Examples

``` r
# Example 1: Full partition (inter + intra via residual)
# Model: (1|track) + (1|observer) + (1|observer:track).
# Var(track)=biology; Var(observer)=inter-observer;
# Var(observer:track)=observer×track; Residual=intra-observer.
tps <- system.file("extdata", "PaluxyRiverObsRep.tps",
                   package = "QuAnTeTrack")
RL <- rep(c("R","L"), length.out = 34)
trks <- tps_to_track(file = tps, scale = 0.004341493,
                     R.L.side = RL, missing = FALSE, NAs = NULL)

n  <- length(trks$Trajectories)
md <- data.frame(
  track    = rep(c("T01","T02"), length.out = n),
  observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
  stringsAsFactors = FALSE
)
md$replica <- ave(seq_len(n),
                  interaction(md$observer, md$track, drop = TRUE),
                  FUN = seq_along)
md <- md[, c("replica","track","observer")]

vars <- c("Distance","Straightness","TurnAng","PaceAng")
res_full <- observer_error_partitioning(data = trks,
                                        metadata = md,
                                        variables = vars)
#> Error in lme4::lFormula(formula = fm, data = dat, REML = TRUE, control = structure(list(    optimizer = "bobyqa", restart_edge = TRUE, boundary.tol = 1e-05,     calc.derivs = NULL, use.last.params = FALSE, checkControl = list(        autoscale = NULL, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "stop",         check.nlev.gtreq.5 = "ignore", check.nlev.gtr.1 = "stop",         check.nobs.vs.nRE = "stop", check.rankX = "message+drop.cols",         check.scaleX = "warning", check.formula.LHS = "stop"),     checkConv = list(check.conv.nobsmax = 10000, check.conv.nparmax = 10,         check.conv.grad = list(action = "warning", tol = 0.002,             relTol = NULL), check.conv.singular = list(action = "ignore",             tol = 1e-04), check.conv.hess = list(action = "warning",             tol = 1e-06)), optCtrl = list(maxfun = 1e+05)), class = c("lmerControl", "merControl"))): 0 (non-NA) cases
res_full$qc
#> Error: object 'res_full' not found

# Example 2: Intra-only, single observer; residual = intra-observer
md_intra <- md
md_intra$observer <- "obs1"
res_intra <- observer_error_partitioning(data = trks,
                                         metadata = md_intra,
                                         variables = vars)
#> Error in lme4::lFormula(formula = fm, data = dat, REML = TRUE, control = structure(list(    optimizer = "bobyqa", restart_edge = TRUE, boundary.tol = 1e-05,     calc.derivs = NULL, use.last.params = FALSE, checkControl = list(        autoscale = NULL, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "stop",         check.nlev.gtreq.5 = "ignore", check.nlev.gtr.1 = "stop",         check.nobs.vs.nRE = "stop", check.rankX = "message+drop.cols",         check.scaleX = "warning", check.formula.LHS = "stop"),     checkConv = list(check.conv.nobsmax = 10000, check.conv.nparmax = 10,         check.conv.grad = list(action = "warning", tol = 0.002,             relTol = NULL), check.conv.singular = list(action = "ignore",             tol = 1e-04), check.conv.hess = list(action = "warning",             tol = 1e-06)), optCtrl = list(maxfun = 1e+05)), class = c("lmerControl", "merControl"))): 0 (non-NA) cases
res_intra$qc
#> Error: object 'res_intra' not found

# Example 3: Inter-only, single replicate; estimates inter-observer
md_inter <- md
md_inter$replica <- 1L
res_inter <- observer_error_partitioning(data = trks,
                                         metadata = md_inter,
                                         variables = vars)
#> Error in lme4::lFormula(formula = fm, data = dat, REML = TRUE, control = structure(list(    optimizer = "bobyqa", restart_edge = TRUE, boundary.tol = 1e-05,     calc.derivs = NULL, use.last.params = FALSE, checkControl = list(        autoscale = NULL, check.nobs.vs.rankZ = "ignore", check.nobs.vs.nlev = "stop",         check.nlev.gtreq.5 = "ignore", check.nlev.gtr.1 = "stop",         check.nobs.vs.nRE = "stop", check.rankX = "message+drop.cols",         check.scaleX = "warning", check.formula.LHS = "stop"),     checkConv = list(check.conv.nobsmax = 10000, check.conv.nparmax = 10,         check.conv.grad = list(action = "warning", tol = 0.002,             relTol = NULL), check.conv.singular = list(action = "ignore",             tol = 1e-04), check.conv.hess = list(action = "warning",             tol = 1e-06)), optCtrl = list(maxfun = 1e+05)), class = c("lmerControl", "merControl"))): 0 (non-NA) cases
res_inter$qc
#> Error: object 'res_inter' not found

# Example 4: Circular metrics only (angles via sine–cosine embedding)
res_ang <- observer_error_partitioning(data = trks,
                                       metadata = md,
                                       variables = c("TurnAng","PaceAng"))
res_ang$summary
#>           variable      component     variance      percent      R2_c
#> TurnAng.1  TurnAng observer:track 9.047652e-11 1.974434e-05 0.9997673
#> TurnAng.2  TurnAng       observer 0.000000e+00 0.000000e+00 0.9997673
#> TurnAng.3  TurnAng          track 4.581145e-04 9.997256e+01 0.9997673
#> TurnAng.4  TurnAng       Residual 1.256725e-07 2.742503e-02 0.9997673
#> PaceAng.1  PaceAng observer:track 2.082826e-50 1.531558e-47 0.9990678
#> PaceAng.2  PaceAng       observer 1.524969e-17 1.121350e-14 0.9990678
#> PaceAng.3  PaceAng          track 1.355629e-01 9.968306e+01 0.9990678
#> PaceAng.4  PaceAng       Residual 4.310243e-04 3.169437e-01 0.9990678
```
