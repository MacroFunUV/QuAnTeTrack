# Calculate combined probabilities of similarity or intersection metrics of tracks

`combined_prob()` calculates the combined probabilities of similarity
and intersection metrics derived from different models. The function
uses simulation data to extract *p*-values, providing insight into the
significance of combined metrics across various similarity assessments.

## Usage

``` r
combined_prob(data, metrics = NULL)
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

The `combined_prob()` function combines *p*-values derived from multiple
similarity metric tests and intersection tests. It calculates the
combined *p*-values by assessing the probability of observing the
combined metrics across simulated datasets. Pairwise *p*-values use a
Monte Carlo tail test with the (+1) correction and are FDR-adjusted with
Benjamini–Hochberg (BH) across unique trajectory pairs.

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
#> 2025-12-29 22:14:28.810195 Iteration 1
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 12.98785
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:28.813155 Iteration 2
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 18.28923
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:28.815881 Iteration 3
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 17.05159
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet1 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#> 2025-12-29 22:14:29.385976 Iteration 1
#>  
#> Frechet metric
#>         Track_1   Track_2
#> Track_1      NA 0.7781548
#> Track_2      NA        NA
#> ------------------------------------
#> 2025-12-29 22:14:29.558083 Iteration 2
#>  
#> Frechet metric
#>         Track_1   Track_2
#> Track_1      NA 0.7738161
#> Track_2      NA        NA
#> ------------------------------------
#> 2025-12-29 22:14:29.726522 Iteration 3
#>  
#> Frechet metric
#>         Track_1   Track_2
#> Track_1      NA 0.7933628
#> Track_2      NA        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int1 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1,
  origin.permutation = "None")
#> 2025-12-29 22:14:29.755168 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 22:14:29.766293 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 22:14:29.777125 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
combined_prob(PaluxyRiver, metrics = list(DTW1, Frechet1, int1))
#> $P_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $P_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $P_values_global
#> [1] 0.25
#> 

# Example 2: "Constrained" model and similarity metrics.
s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
DTW2 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2,
  superposition = "None")
#> 2025-12-29 22:14:29.792181 Iteration 1
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 121.1052
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:29.794877 Iteration 2
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 62.28218
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:29.797555 Iteration 3
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 20.55951
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet2 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2,
  superposition = "None")
#> 2025-12-29 22:14:30.120713 Iteration 1
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 3.765697
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:30.277174 Iteration 2
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 3.752331
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:30.449051 Iteration 3
#>  
#> Frechet metric
#>         Track_1   Track_2
#> Track_1      NA 0.9300576
#> Track_2      NA        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int2 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2,
  origin.permutation = "Min.Box")
#> 2025-12-29 22:14:31.078809 Permutation 1
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 22:14:31.395402 Permutation 2
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 22:14:31.504161 Permutation 3
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2025-12-29 22:14:31.515815 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 22:14:31.526722 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 22:14:31.537618 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
combined_prob(PaluxyRiver, metrics = list(DTW2, Frechet2, int2))
#> $P_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $P_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $P_values_global
#> [1] 0.25
#> 

# Example 3: "Unconstrained" model and similarity metrics.
s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
DTW3 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3,
  superposition = "None")
#> 2025-12-29 22:14:31.552076 Iteration 1
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 647.6714
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:31.554773 Iteration 2
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 265.3258
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:31.557434 Iteration 3
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 57.93855
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
Frechet3 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3,
  superposition = "None")
#> 2025-12-29 22:14:31.762534 Iteration 1
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 24.93024
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:31.901724 Iteration 2
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 10.06423
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 22:14:32.057201 Iteration 3
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 3.036826
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
int3 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3,
  origin.permutation = "Conv.Hull")
#> 2025-12-29 22:14:32.42854 Permutation 1
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2025-12-29 22:14:32.775352 Permutation 2
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2025-12-29 22:14:32.882187 Permutation 3
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2025-12-29 22:14:32.893986 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 22:14:32.904796 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 22:14:32.915671 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
combined_prob(PaluxyRiver, metrics = list(DTW3, Frechet3, int3))
#> $P_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $P_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $P_values_global
#> [1] 0.25
#> 
```
