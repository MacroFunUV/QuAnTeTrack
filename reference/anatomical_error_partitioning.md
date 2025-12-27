# Partition anatomical-fidelity uncertainty (simulation-based)

`anatomical_error_partitioning()` quantifies how much variance in
trackway parameters is due to track landmark placement uncertainty
(degree of anatomical fidelity) versus genuine between-track
(biological) differences, using simulation only (no observer effects).

## Usage

``` r
anatomical_error_partitioning(
  data,
  error_radius,
  variables = c("TurnAng", "sdTurnAng", "Distance", "Length", "StLength", "sdStLength",
    "Sinuosity", "Straightness", "TrackWidth", "PaceAng"),
  n_sim = 200,
  distribution = c("uniform", "gaussian"),
  seed = NULL
)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive tracks.

  - **`tracks`**: A list of data frames containing track coordinates,
    metadata (e.g., image reference, ID), and a marker indicating
    whether the track is actual or inferred.

- error_radius:

  Numeric; single radius (m) applied to all tracks or a numeric vector
  of length `length(data$tracks)` with per-track tolerances.

- variables:

  A character vector of parameters to analyze. Default is
  `c("TurnAng","sdTurnAng","Distance","Length","StLength","sdStLength", "Sinuosity","Straightness","TrackWidth","PaceAng")`.

- n_sim:

  Integer; number of Monte Carlo jitter simulations per track. Default
  is `200`.

- distribution:

  A character string indicating the jitter model. Options are
  `"uniform"` (points uniformly within a disk of radius r; default) or
  `"gaussian"` (SD = r/2 on X and Y, truncated at 3 SD).

- seed:

  Optional integer for reproducibility. Default is `NULL`.

## Value

An `"error_partitioning"` R object consisting of a list containing the
following elements:

- summary:

  Data frame of variance components per variable. Columns: `variable`,
  `component`, `variance`, `percent`, and \\R^2_c\\. Components are
  `track` (biological), `anatomical` (within-track simulation variance),
  and `Residual` (from the model on simulation means).

- snr:

  Data frame with the signal-to-noise ratio per variable. Columns:
  `variable`, `bio_component` (always `"track"`), `bio_var`,
  `error_var`, and `SNR`. The ratio is \$\$SNR = Var(track) /
  \[Var(anatomical) + Var(Residual)\].\$\$

- qc:

  Compact quality-control table per variable with `SNR`, a qualitative
  `SNR_rating` (`"weak"`, `"moderate"`, `"strong"`), the `top_component`
  by percent, its `top_percent`, \\R^2_c\\, and an `observer_note`
  (always `"anatomical-only"` for this function).

- models:

  `NULL`. Placeholder kept for consistency with the observer-inclusive
  workflow. No observer models are fitted here.

- analysis_table:

  Long data frame with all simulated values for every track across all
  Monte Carlo runs (variables joined with track IDs and simulation
  labels).

- formulae:

  A short note describing the simulation-only decomposition
  (within-track anatomical variance from simulations, plus a
  random-intercept model on simulation means: \\y \sim 1 + (1\|track)\\;
  angles handled via sine–cosine).

- track_names:

  Character vector of internal track labels used in the analysis.

## Details

The function estimates how much of the variability in trackway metrics
can be attributed to positional uncertainty in track landmarks, rather
than to true biological differences among trackways. It asks, in
essence: "If landmarks were slightly misplaced by up to `error_radius`,
how much would that affect our computed parameters?" For each trackway,
landmarks are repeatedly perturbed within the specified spatial
tolerance, the medial trajectory is reconstructed, and parameters are
recalculated across Monte Carlo simulations. The resulting variance in
each metric is then decomposed into a between-track component
(biological signal) and a within-track component (anatomical noise).

The `error_radius` represents expected imprecision in reference-point
digitization (e.g., preservation quality, erosion, deformation, or
subjective landmark interpretation). By specifying a realistic
tolerance, users can evaluate how sensitive each metric is to anatomical
or taphonomic uncertainty and identify parameters that remain stable
despite imperfect preservation.

For every track and simulation, landmarks are randomly displaced within
a circular area of radius `r`, using either a uniform model
(`"uniform"`) or a truncated Gaussian on `X` and `Y` with `SD = r/2`
(`"gaussian"`, truncated at ±3 SD). Each perturbed landmark set is
converted to a medial trajectory (midpoints between consecutive
footprints, as in
[`tps_to_track()`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md)),
and parameters are recalculated via
[`track_param()`](https://macrofunuv.github.io/QuAnTeTrack/reference/track_param.md).
For each metric \\Y\\, total variance is partitioned via the law of
total variance into:

- **track (biological)** — variance of simulation means among tracks;

- **anatomical** — mean of within-track variances across simulations;

- **Residual** — residual from the random-intercept model on simulation
  means (no observer terms in this simulation-only setup).

Angular variables (`TurnAng`, `PaceAng`) are embedded in sine–cosine
space to handle circularity and avoid 0°/360° discontinuities (Fisher,
1995).

A signal-to-noise ratio quantifies robustness: \$\$\mathrm{SNR} =
\frac{\mathrm{Var}(\mathrm{track})} {\mathrm{Var}(\mathrm{anatomical}) +
\mathrm{Var}(\mathrm{Residual})}.\$\$ As a rule of thumb:

- **SNR \< 1** — anatomical error exceeds biology (weak);

- **SNR ~ 1–2** — biology and anatomical error are comparable
  (moderate);

- **SNR \> 2** — biology dominates (strong).

A compact QC table reports SNR, the qualitative rating, the top
component by percent, and the conditional \\R^2_c\\.

## Logo

![](figures/Logo.png)

## References

Fisher, N. I. (1995). Statistical Analysis of Circular Data. Cambridge
University Press.

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`track_param`](https://macrofunuv.github.io/QuAnTeTrack/reference/track_param.md),
[`subset_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/subset_track.md)

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
# Example 1: PaluxyRiver, small jitter (2 cm), uniform noise
# Expect relatively stable metrics and higher SNR.
set.seed(1)
ep_small <- anatomical_error_partitioning(
  data         = PaluxyRiver,
  error_radius = 0.02,
  variables    = c("Distance", "Straightness", "TurnAng"),
  n_sim        = 10
)
ep_small$qc
#>       variable        SNR SNR_rating top_component top_percent      R2_c
#> 1     Distance 2734.80404     strong         track    99.96345 0.9996345
#> 2 Straightness   22.61413     strong         track    95.76525 0.9576525
#> 3      TurnAng 1471.30665     strong         track    99.93208 0.9993208
#>                observer_estimable
#> 1 no (anatomical-only, sim-based)
#> 2 no (anatomical-only, sim-based)
#> 3 no (anatomical-only, sim-based)

# Example 2: PaluxyRiver, larger jitter (8 cm) to see SNR drop
# Inflate anatomical uncertainty; SNR should typically decrease.
set.seed(1)
ep_large <- anatomical_error_partitioning(
  data         = PaluxyRiver,
  error_radius = 0.08,
  variables    = c("Distance", "Straightness", "TurnAng"),
  n_sim        = 10
)
ep_large$qc
#>       variable        SNR SNR_rating top_component top_percent      R2_c
#> 1     Distance 165.040390     strong         track    99.39774 0.9939774
#> 2 Straightness   2.187696     strong         track    68.62938 0.6862938
#> 3      TurnAng  90.802197     strong         track    98.91070 0.9891070
#>                observer_estimable
#> 1 no (anatomical-only, sim-based)
#> 2 no (anatomical-only, sim-based)
#> 3 no (anatomical-only, sim-based)

# Example 3: MountTom subset + Gaussian jitter (3 cm, truncated at ±3 SD)
# Demonstrates alternative noise model and a reduced dataset.
sbMountTom <- subset_track(
  MountTom,
  tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18)
)
set.seed(2)
ep_gauss <- anatomical_error_partitioning(
  data         = sbMountTom,
  error_radius = 0.03,
  variables    = c("StLength", "sdStLength", "PaceAng"),
  n_sim        = 10,
  distribution = "gaussian"
)
ep_gauss$snr
#>              variable bio_component     bio_var    error_var        SNR
#> StLength     StLength         track 0.058217558 1.031312e-05 5644.99846
#> sdStLength sdStLength         track 0.001165249 2.973018e-05   39.19415
#> PaceAng       PaceAng         track 0.005343677 2.543369e-04   21.01023

# Example 4: PaluxyRiver with heterogeneous per-track tolerances
# Alternate 2 cm / 5 cm across tracks; highlights mixed preservation quality.
set.seed(3)
rads <- rep(c(0.02, 0.05),
            length.out = length(PaluxyRiver$Footprints))
ep_het <- anatomical_error_partitioning(
  data         = PaluxyRiver,
  error_radius = rads,
  variables    = c("TrackWidth", "Sinuosity", "PaceAng"),
  n_sim        = 10
)
ep_het$summary
#>                variable  component     variance    percent      R2_c
#> TrackWidth.1 TrackWidth      track 3.740549e-02 99.7499311 0.9974993
#> TrackWidth.2 TrackWidth anatomical 9.377402e-05  0.2500689 0.9974993
#> TrackWidth.3 TrackWidth   Residual 0.000000e+00  0.0000000 0.9974993
#> Sinuosity.1   Sinuosity      track 3.086845e-03 98.6831363 0.9868314
#> Sinuosity.2   Sinuosity anatomical 4.119198e-05  1.3168637 0.9868314
#> Sinuosity.3   Sinuosity   Residual 0.000000e+00  0.0000000 0.9868314
#> PaceAng.1       PaceAng      track 8.498442e-02 99.7797710 0.9977977
#> PaceAng.2       PaceAng anatomical 1.875734e-04  0.2202290 0.9977977
#> PaceAng.3       PaceAng   Residual 0.000000e+00  0.0000000 0.9977977

# Example 5: Angles only (circular handling via sine–cosine embedding)
# Focus on orientation metrics; output includes averaged sin/cos components.
set.seed(4)
ep_ang <- anatomical_error_partitioning(
  data         = PaluxyRiver,
  error_radius = 0.04,
  variables    = c("TurnAng", "PaceAng"),
  n_sim        = 10
)
ep_ang$summary
#>           variable  component     variance    percent      R2_c
#> TurnAng.1  TurnAng      track 4.283363e-04 99.8192789 0.9981928
#> TurnAng.2  TurnAng anatomical 7.754955e-07  0.1807211 0.9981928
#> TurnAng.3  TurnAng   Residual 0.000000e+00  0.0000000 0.9981928
#> PaceAng.1  PaceAng      track 8.475987e-02 99.7051248 0.9970512
#> PaceAng.2  PaceAng anatomical 2.506750e-04  0.2948752 0.9970512
#> PaceAng.3  PaceAng   Residual 0.000000e+00  0.0000000 0.9970512
```
