# Calculate combined probabilities of similarity or intersection metrics of tracks

`combined_prob()` calculates the combined probabilities of similarity
and intersection metrics derived from different models. The function
uses simulation data to extract *p*-values, providing insight into the
significance of combined metrics across various similarity assessments.

## Usage

``` r
combined_prob(data, metrics = NULL, H1 = c("Higher", "Lower"))
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- metrics:

  A list of `track similarity` and/or `track intersection` R objects
  derived from different tests. All tests must be based on the same
  number of simulations.

- H1:

  Alternative hypothesis for intersection metrics. One of `"Higher"` or
  `"Lower"`.

## Value

A list containing:

- P_values:

  A matrix of raw *p*-values for the combined metrics across all
  trajectories (pairwise).

- P_values_BH:

  A matrix of BH-adjusted *p*-values for the combined metrics
  (pairwise).

- P_values_global:

  A single numeric value: the overall probability of observing the
  combined metrics across all pairs of trajectories.

## Details

The `combined_prob()` function combines evidence from multiple
trajectory-based metrics, including distance-based similarity measures
(e.g., DTW and Fréchet) and intersection metrics, using simulation-based
hypothesis testing.

For distance-based similarity metrics (such as DTW and Fréchet), smaller
values indicate greater similarity between trajectories. Accordingly,
pairwise *p*-values are computed using a lower-tail Monte Carlo test
with the (+1) correction (Phipson & Smyth, 2010): \$\$p = (1 +
\\\\D\_{sim} \le D\_{obs}\\) / (nsim + 1)\$\$ This tests whether the
observed trajectories are at least as similar as expected under the null
model of random movement.

For intersection metrics, the direction of the test depends on the
alternative hypothesis specified by `H1`. When `H1 = "Lower"`, the test
evaluates whether the observed number of intersections is significantly
lower than expected under random simulations (e.g., coordinated or
gregarious movement). When `H1 = "Higher"`, the test evaluates whether
the observed number of intersections is significantly higher than
expected (e.g., chasing or predatory interactions).

Combined pairwise *p*-values are obtained by jointly evaluating all
metrics across simulations, identifying those simulation replicates in
which all metrics are simultaneously as extreme as the observed values,
given their respective directions. A Monte Carlo test with the (+1)
correction is used, and the resulting pairwise *p*-values are adjusted
for multiple comparisons using the Benjamini–Hochberg (BH) procedure.

In addition, a global combined *p*-value is computed based on a single
global statistic summarizing all pairwise comparisons, providing an
overall assessment of whether the observed set of trajectories departs
from random expectations across all metrics considered.

## Logo

![](figures/Logo.png)

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`simulate_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/simulate_track.md),
[`track_intersection`](https://macrofunuv.github.io/QuAnTeTrack/reference/track_intersection.md),
[`simil_DTW_metric`](https://macrofunuv.github.io/QuAnTeTrack/reference/simil_DTW_metric.md),
[`simil_Frechet_metric`](https://macrofunuv.github.io/QuAnTeTrack/reference/simil_Frechet_metric.md)

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
# Example 1: "Directed" model and similarity metrics.
s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
DTW1 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#> 2026-01-21 09:03:31.501205 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 7.767709
#> Track_2 7.767709       NA
#> ------------------------------------
#> 2026-01-21 09:03:31.508037 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 10.91272
#> Track_2 10.91272       NA
#> ------------------------------------
#> 2026-01-21 09:03:31.514805 Iteration 3
#>  
#> DTW metric
#>         Track_1 Track_2
#> Track_1      NA 9.99559
#> Track_2 9.99559      NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet1 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#> 2026-01-21 09:03:31.888176 Iteration 1
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.7781548
#> Track_2 0.7781548        NA
#> ------------------------------------
#> 2026-01-21 09:03:32.066994 Iteration 2
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.7738161
#> Track_2 0.7738161        NA
#> ------------------------------------
#> 2026-01-21 09:03:32.24394 Iteration 3
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.7933628
#> Track_2 0.7933628        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int1 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1,
  origin.permutation = "None")
#> 2026-01-21 09:03:32.268628 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-01-21 09:03:32.279761 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-01-21 09:03:32.290744 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
combined_prob(PaluxyRiver, metrics = list(DTW1, Frechet1, int1), H1 = "Lower")
#> $P_values
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2     0.5      NA
#> 
#> $P_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2     0.5      NA
#> 
#> $P_values_global
#> [1] 0.5
#> 

# Example 2: "Constrained" model and similarity metrics.
s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
DTW2 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2,
  superposition = "None")
#> 2026-01-21 09:03:32.314099 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 67.17671
#> Track_2 67.17671       NA
#> ------------------------------------
#> 2026-01-21 09:03:32.325573 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 34.32521
#> Track_2 34.32521       NA
#> ------------------------------------
#> 2026-01-21 09:03:32.33231 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 11.66344
#> Track_2 11.66344       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet2 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2,
  superposition = "None")
#> 2026-01-21 09:03:32.66673 Iteration 1
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 3.765697
#> Track_2 3.765697       NA
#> ------------------------------------
#> 2026-01-21 09:03:32.829095 Iteration 2
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 3.752331
#> Track_2 3.752331       NA
#> ------------------------------------
#> 2026-01-21 09:03:33.006393 Iteration 3
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.9300576
#> Track_2 0.9300576        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int2 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2,
  origin.permutation = "Min.Box")
#> 2026-01-21 09:03:33.57337 Permutation 1
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2026-01-21 09:03:33.708011 Permutation 2
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2026-01-21 09:03:33.759503 Permutation 3
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2026-01-21 09:03:33.771085 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-01-21 09:03:33.781963 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-01-21 09:03:33.792753 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
combined_prob(PaluxyRiver, metrics = list(DTW2, Frechet2, int2), H1 = "Lower")
#> $P_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $P_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $P_values_global
#> [1] 0.25
#> 

# Example 3: "Unconstrained" model and similarity metrics.
s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
DTW3 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3,
  superposition = "None")
#> 2026-01-21 09:03:33.815015 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 27.66064
#> Track_2 27.66064       NA
#> ------------------------------------
#> 2026-01-21 09:03:33.821781 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 300.0043
#> Track_2 300.0043       NA
#> ------------------------------------
#> 2026-01-21 09:03:33.828178 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 109.3543
#> Track_2 109.3543       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet3 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3,
  superposition = "None")
#> 2026-01-21 09:03:34.186305 Iteration 1
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 1.540872
#> Track_2 1.540872       NA
#> ------------------------------------
#> 2026-01-21 09:03:34.212247 Iteration 2
#>  
#> Frechet metric
#>         Track_1 Track_2
#> Track_1      NA      -1
#> Track_2      -1      NA
#> ------------------------------------
#> 2026-01-21 09:03:34.34404 Iteration 3
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 9.937678
#> Track_2 9.937678       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int3 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3,
  origin.permutation = "Conv.Hull")
#> 2026-01-21 09:03:34.431226 Permutation 1
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2026-01-21 09:03:34.487242 Permutation 2
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2026-01-21 09:03:34.731728 Permutation 3
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2026-01-21 09:03:34.74358 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-01-21 09:03:34.754545 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-01-21 09:03:34.765285 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
combined_prob(PaluxyRiver, metrics = list(DTW3, Frechet3, int3), H1 = "Lower")
#> $P_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $P_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $P_values_global
#> [1] 0.25
#> 
```
