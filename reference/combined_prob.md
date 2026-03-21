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
#> 2026-03-21 22:02:12.20057 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 9.909364
#> Track_2 9.909364       NA
#> ------------------------------------
#> 2026-03-21 22:02:12.213445 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 12.33171
#> Track_2 12.33171       NA
#> ------------------------------------
#> 2026-03-21 22:02:12.220206 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 13.77713
#> Track_2 13.77713       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet1 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#> 2026-03-21 22:02:12.625372 Iteration 1
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.7791748
#> Track_2 0.7791748        NA
#> ------------------------------------
#> 2026-03-21 22:02:12.80519 Iteration 2
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.8071543
#> Track_2 0.8071543        NA
#> ------------------------------------
#> 2026-03-21 22:02:12.981525 Iteration 3
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.7756086
#> Track_2 0.7756086        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int1 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1,
  origin.permutation = "None")
#> 2026-03-21 22:02:13.01014 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-03-21 22:02:13.021431 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-03-21 22:02:13.032831 Iteration 3
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
#> 2026-03-21 22:02:13.056423 Iteration 1
#>  
#> DTW metric
#>         Track_1 Track_2
#> Track_1      NA 72.6468
#> Track_2 72.6468      NA
#> ------------------------------------
#> 2026-03-21 22:02:13.06316 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 48.04177
#> Track_2 48.04177       NA
#> ------------------------------------
#> 2026-03-21 22:02:13.069801 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 11.36891
#> Track_2 11.36891       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet2 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2,
  superposition = "None")
#> 2026-03-21 22:02:13.414616 Iteration 1
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 4.342065
#> Track_2 4.342065       NA
#> ------------------------------------
#> 2026-03-21 22:02:13.582997 Iteration 2
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 5.041493
#> Track_2 5.041493       NA
#> ------------------------------------
#> 2026-03-21 22:02:13.985888 Iteration 3
#>  
#> Frechet metric
#>           Track_1   Track_2
#> Track_1        NA 0.7737559
#> Track_2 0.7737559        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int2 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2,
  origin.permutation = "Min.Box")
#> 2026-03-21 22:02:14.112799 Permutation 1
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2026-03-21 22:02:14.4418 Permutation 2
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2026-03-21 22:02:14.495591 Permutation 3
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2026-03-21 22:02:14.507159 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-03-21 22:02:14.517945 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-03-21 22:02:14.528648 Iteration 3
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
#> 2026-03-21 22:02:14.550829 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 327.7606
#> Track_2 327.7606       NA
#> ------------------------------------
#> 2026-03-21 22:02:14.557191 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 285.6154
#> Track_2 285.6154       NA
#> ------------------------------------
#> 2026-03-21 22:02:14.563591 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 257.7622
#> Track_2 257.7622       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet3 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3,
  superposition = "None")
#> 2026-03-21 22:02:14.96716 Iteration 1
#>  
#> Frechet metric
#>         Track_1 Track_2
#> Track_1      NA 22.6934
#> Track_2 22.6934      NA
#> ------------------------------------
#> 2026-03-21 22:02:15.009372 Iteration 2
#>  
#> Frechet metric
#>         Track_1 Track_2
#> Track_1      NA 21.7064
#> Track_2 21.7064      NA
#> ------------------------------------
#> 2026-03-21 22:02:15.051733 Iteration 3
#>  
#> Frechet metric
#>          Track_1  Track_2
#> Track_1       NA 20.80205
#> Track_2 20.80205       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int3 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3,
  origin.permutation = "Conv.Hull")
#> 2026-03-21 22:02:15.14455 Permutation 1
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2026-03-21 22:02:15.398383 Permutation 2
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2026-03-21 22:02:15.447344 Permutation 3
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2026-03-21 22:02:15.458887 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-03-21 22:02:15.469596 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2026-03-21 22:02:15.480177 Iteration 3
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
