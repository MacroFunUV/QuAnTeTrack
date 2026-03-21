# Partition intra- and inter-observer uncertainty (replicate-based)

`observer_error_partitioning()` quantifies how much variance in trackway
parameters is due to observer and replicate effects
(inter-/intra-observer) versus genuine between-trackway (biological)
differences, using replicated digitizations and mixed-effects models.

## Usage

``` r
observer_error_partitioning(
  data,
  metadata,
  veltrack = NULL,
  variables,
  gauge_size = NA
)
```

## Arguments

- data:

  A `trackway` R object. See
  [`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md).

- metadata:

  A data frame with one row per entry. Required columns:

  - `replica` — replicate index within (`trackway`, `observer`).

  - `trackway` — trackway ID.

  - `observer` — observer ID.

- veltrack:

  Optional. A `trackway velocity` R object (list of lists) as returned
  by
  [`velocity_track()`](https://macrofunuv.github.io/QuAnTeTrack/reference/velocity_track.md),
  required only if velocity variables are requested.

- variables:

  A character vector specifying the movement parameters to be used.
  Valid parameter names include: `"TurnAng"`, `"sdTurnAng"`,
  `"PathLen"`, `"BeelineLen"`, `"StLength"`, `"sdStLength"`,
  `"StrideLen"`, `"PaceLen"`, `"Sinuosity"`, `"Straightness"`,
  `"TrackWidth"`, `"Gauge"`, `"PaceAng"`, `"StepAng"`, `"Velocity"`,
  `"sdVelocity"`, `"MaxVelocity"`, `"MinVelocity"`. See
  [`track_param`](https://macrofunuv.github.io/QuAnTeTrack/reference/track_param.md)
  and
  [`cluster_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/cluster_track.md)
  for details.

- gauge_size:

  Numeric. Pes/manus length (or width) used to compute Gauge as
  `Trackway_width / gauge_size`. See
  [`track_param`](https://macrofunuv.github.io/QuAnTeTrack/reference/track_param.md)
  for details.

## Value

An `"error_partitioning"` R object consisting of a list containing the
following elements:

- summary:

  Data frame of variance components per variable. Columns: `variable`,
  `component`, `variance`, `percent`, and \\R^2_c\\. Components reflect
  what was estimable (e.g., `trackway`, `observer`, `observer:trackway`,
  `Residual`).

- snr:

  Data frame with the signal-to-noise ratio per variable. Columns:
  `variable`, `bio_component` (always `"trackway"`), `bio_var`,
  `error_var`, and `SNR`. The SNR is \$\$SNR = Var(trackway) /
  \[Var(observer) + Var(observer:trackway) + Var(Residual)\].\$\$

- qc:

  Compact quality-control table per variable with `SNR`, a qualitative
  `SNR_rating` (`"weak"`, `"moderate"`, `"strong"`), the `top_component`
  by percent, its `top_percent`, \\R^2_c\\, and whether the observer
  term was estimable (`observer_estimable`).

- analysis_table:

  Data frame used in the analysis: selected variables joined with
  `metadata` (columns typically include `replica`, `trackway`, and
  `observer`).

- models:

  Fitted mixed-effects models used to extract variance components. For
  linear variables, each entry stores the `lmer` fit and its variance
  components. For angular variables, each entry is a list with separate
  `sin` and `cos` fits.

- formulae:

  Model formula(e) per variable. For angular variables, a list
  containing the `sin` and `cos` formulas.

- trackway_names:

  Character vector of internal trackway labels used in the analysis.

## Details

This function partitions the variability of each selected metric into
components attributable to **biology** (consistent differences among
trackways) and to the **human sampling process**: differences between
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

For each metric the model fits random intercepts for `trackway`,
`observer`, and their `observer:trackway` interaction.
Replicate-to-replicate scatter within the same (`observer`, `trackway`)
cell contributes to the residual (this is the intra-observer sampling
error). Random terms that appear with a single level are dropped
automatically.

\$\$y \sim 1 + (1\|trackway) + (1\|observer) +
(1\|observer:trackway).\$\$

Angular variables (`TurnAng`, `PaceAng`) are handled in sine–cosine
space to respect circularity and avoid 0°/360° discontinuities: the same
random-effects structure is fit separately to \\\sin(\theta)\\ and
\\\cos(\theta)\\, and the two variance decompositions are averaged
(Fisher, 1995).

Robustness is summarized by the signal-to-noise ratio \$\$SNR =
\frac{Var(trackway)}{Var(observer) + Var(observer:trackway) +
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
  replicates per observer. Both `observer` and `observer:trackway` are
  estimable; replicate scatter is residual.

- *Intra-only*: set `metadata$observer` to the same value for all rows
  and provide `replica` with \>1 levels. Only `trackway` remains;
  replicate scatter is residual.

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
in the same file), with each replicate represented as its own trackway
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
# Model: (1|trackway) + (1|observer) + (1|observer:trackway).
# Var(trackway)=biology; Var(observer)=inter-observer;
# Var(observer:trackway)=observer×trackway; Residual=intra-observer.
tps <- system.file("extdata", "PaluxyRiverObsRep.tps",
                   package = "QuAnTeTrack")
RL <- rep(c("R","L"), length.out = 34)
trks <- tps_to_track(file = tps, scale = 0.004341493,
                     R.L.side = RL, missing = FALSE, NAs = NULL)

n  <- length(trks$Trajectories)
md <- data.frame(
  trackway    = rep(c("T01","T02"), length.out = n),
  observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
  stringsAsFactors = FALSE
)
md$replica <- ave(seq_len(n),
                  interaction(md$observer, md$trackway, drop = TRUE),
                  FUN = seq_along)
md <- md[, c("replica","trackway","observer")]

vars <- c("BeelineLen","Straightness","TurnAng","PaceAng")
res_full <- observer_error_partitioning(data = trks,
                                        metadata = md,
                                        variables = vars)
res_full$qc
#>       variable        SNR SNR_rating top_component top_percent      R2_c
#> 1   BeelineLen 14362.3060     strong      trackway    99.99304 0.9999304
#> 2 Straightness   158.6878     strong      trackway    99.37378 0.9937378
#> 3      TurnAng  3642.6814     strong      trackway    99.97256 0.9997673
#> 4      PaceAng   314.5134     strong      trackway    99.68306 0.9990678
#>   observer_estimable
#> 1                yes
#> 2                yes
#> 3                yes
#> 4                yes

# Example 2: Intra-only, single observer; residual = intra-observer
md_intra <- md
md_intra$observer <- "obs1"
res_intra <- observer_error_partitioning(data = trks,
                                         metadata = md_intra,
                                         variables = vars)
res_intra$qc
#>       variable        SNR SNR_rating top_component top_percent      R2_c
#> 1   BeelineLen 14362.5544     strong      trackway    99.99304 0.9999304
#> 2 Straightness   158.6867     strong      trackway    99.37377 0.9937377
#> 3      TurnAng  3643.9188     strong      trackway    99.97256 0.9997582
#> 4      PaceAng   314.5132     strong      trackway    99.68306 0.9990678
#>   observer_estimable
#> 1       no (1 level)
#> 2       no (1 level)
#> 3       no (1 level)
#> 4       no (1 level)

# Example 3: Inter-only, single replicate; estimates inter-observer
md_inter <- md
md_inter$replica <- 1L
res_inter <- observer_error_partitioning(data = trks,
                                         metadata = md_inter,
                                         variables = vars)
res_inter$qc
#>       variable        SNR SNR_rating top_component top_percent      R2_c
#> 1   BeelineLen 14362.3060     strong      trackway    99.99304 0.9999304
#> 2 Straightness   158.6878     strong      trackway    99.37378 0.9937378
#> 3      TurnAng  3642.6814     strong      trackway    99.97256 0.9997673
#> 4      PaceAng   314.5134     strong      trackway    99.68306 0.9990678
#>   observer_estimable
#> 1                yes
#> 2                yes
#> 3                yes
#> 4                yes

# Example 4: Circular metrics only (angles via sine–cosine embedding)
res_ang <- observer_error_partitioning(data = trks,
                                       metadata = md,
                                       variables = c("TurnAng","PaceAng"))
res_ang$summary
#>           variable         component     variance      percent      R2_c
#> TurnAng.1  TurnAng observer:trackway 9.047652e-11 1.974434e-05 0.9997673
#> TurnAng.2  TurnAng          observer 0.000000e+00 0.000000e+00 0.9997673
#> TurnAng.3  TurnAng          trackway 4.581145e-04 9.997256e+01 0.9997673
#> TurnAng.4  TurnAng          Residual 1.256725e-07 2.742503e-02 0.9997673
#> PaceAng.1  PaceAng observer:trackway 2.082826e-50 1.531558e-47 0.9990678
#> PaceAng.2  PaceAng          observer 1.524969e-17 1.121350e-14 0.9990678
#> PaceAng.3  PaceAng          trackway 1.355629e-01 9.968306e+01 0.9990678
#> PaceAng.4  PaceAng          Residual 4.310243e-04 3.169437e-01 0.9990678
```
