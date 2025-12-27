# Similarity metric using Fréchet distance

`simil_Frechet_metric()` computes similarity metrics between two or more
trajectories using the Fréchet distance. It allows for different
superposition methods to align trajectories before calculating the
Fréchet distance metrics. The function also supports testing with
simulations to calculate *p*-values for the Fréchet distance metrics.

## Usage

``` r
simil_Frechet_metric(data, test = FALSE, sim = NULL, superposition = "None")
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- test:

  Logical; if `TRUE`, the function compares the observed Fréchet
  distances against simulated trajectories and calculates *p*-values.
  Default is `FALSE`.

- sim:

  A `track simulation` R object consisting of a list of simulated
  trajectories to use for comparison when `test = TRUE`.

- superposition:

  A character string indicating the method used to align trajectories.
  Options are `"None"`, `"Centroid"`, or `"Origin"`. Default is
  `"None"`.

## Value

A `track similarity` R object consisting of a list containing the
following elements:

- Frechet_distance_metric:

  A numeric matrix of pairwise Fréchet distances between trajectories.
  Each entry represents the Fréchet distance between the corresponding
  pair of trajectories.

- Frechet_distance_metric_p_values:

  (If `test = TRUE`) A numeric matrix of raw pairwise *p*-values,
  computed by Monte Carlo tail tests with the (+1) correction (Phipson &
  Smyth, 2010): \$\$p = (1 + \\\\\text{extreme}\\) / (nsim + 1)\$\$.
  Each entry reflects the probability of observing a Fréchet distance as
  extreme as the observed one, given the null hypothesis of no
  difference.

- Frechet_distance_metric_p_values_BH:

  (If `test = TRUE`) A numeric matrix of Benjamini–Hochberg (BH)
  adjusted *p*-values controlling the false discovery rate (FDR),
  applied across the set of unique pairwise tests (Benjamini & Hochberg,
  1995).

- Frechet_metric_p_values_combined:

  (If `test = TRUE`) A single numeric value representing the combined
  *p*-value across all Fréchet distances (based on the global statistic:
  the sum of pairwise distances). This indicates the overall
  significance of the observed Fréchet distances relative to
  simulations.

- Frechet_distance_metric_simulations:

  (If `test = TRUE`) A list containing matrices of Fréchet distances for
  each simulation iteration, allowing for inspection of the distribution
  of Fréchet distances across randomized scenarios.

## Details

The `simil_Frechet_metric()` function calculates the similarity between
trajectories using the `Frechet()` function from the SimilarityMeasures
package.

The Fréchet distance is a measure of similarity between two curves or
continuous trajectories, which takes into account both the order and
location of points within the trajectories (Besse et al. 2015). The
distance can be described by the analogy of a person walking a dog on an
extendable leash (Aronov et al. 2006). Both the person and the dog move
along their respective trajectories, with each able to adjust their
speed but not retrace their steps. The Fréchet distance is the minimum
leash length required to keep the dog connected to the person throughout
the walk (Cleasby et al., 2019).

Unlike other trajectory comparison techniques, such as Dynamic Time
Warping, the Fréchet distance focuses on the overall shape of the
trajectories rather than matching specific points. As a result, it is
sensitive to noise because all points of the compared trajectories are
considered in its calculation. However, it can still be a powerful tool
for trajectory clustering and comparison, particularly when shape is the
primary concern (Cleasby et al., 2019).

Note that when comparing real trajectories that are very disparate or
those simulated under an unconstrained method, the resulting
trajectories may not be suitable for Fréchet distance calculations. In
such cases, the Fréchet distance is returned as -1 to indicate an
invalid measurement.

The function offers three different superposition methods to align the
trajectories before `Frechet()` is applied:

- `"None"`: No superposition is applied.

- `"Centroid"`: Trajectories are shifted to align based on their
  centroids.

- `"Origin"`: Trajectories are shifted to align based on their starting
  point.

If `test = TRUE`, the function can compute *p*-values by comparing the
observed Fréchet distances with those generated from a set of simulated
trajectories. The *p*-values are calculated for both individual
trajectory pairs and for the entire set of trajectories. Pairwise
*p*-values are computed with a Monte Carlo tail test using the (+1)
correction to avoid zero-values (see Phipson & Smyth, 2010): \$\$p =
(1 + \\\\\text{extreme}\\) / (nsim + 1)\$\$

These raw *p*-values are then adjusted for multiple comparisons across
the set of unique pairwise tests using the Benjamini–Hochberg (BH)
procedure for false discovery rate (FDR) control (Benjamini & Hochberg,
1995). Both the raw and the BH-adjusted *p*-value matrices are returned
in the output object, allowing users to inspect either uncorrected or
corrected results. In addition, a global combined *p*-value is provided,
summarizing the overall deviation from the null across all pairs.

## Logo

![](figures/Logo.png)

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. Journal of
the Royal statistical society: series B (Methodological), 57(1),
289-300.

Cleasby, I. R., Wakefield, E. D., Morrissey, B. J., Bodey, T. W.,
Votier, S. C., Bearhop, S., & Hamer, K. C. (2019). Using time-series
similarity measures to compare animal movement trajectories in ecology.
Behavioral Ecology and Sociobiology, 73, 1-19.

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be
zero: calculating exact P-values when permutations are randomly drawn.
Statistical applications in genetics and molecular biology, 9(1).

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`simulate_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/simulate_track.md),
[`Frechet`](https://rdrr.io/pkg/SimilarityMeasures/man/Frechet.html)

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
# Example 1: Simulating tracks using the "Directed" model and comparing Frechet distance
# in the PaluxyRiver dataset
s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#> 2025-12-27 20:56:00.850077 Iteration 1
#>  
#> Frechet metric
#>         Track_1   Track_2
#> Track_1      NA 0.8806958
#> Track_2      NA        NA
#> ------------------------------------
#> 2025-12-27 20:56:01.040742 Iteration 2
#>  
#> Frechet metric
#>         Track_1   Track_2
#> Track_1      NA 0.7820993
#> Track_2      NA        NA
#> ------------------------------------
#> 2025-12-27 20:56:01.25676 Iteration 3
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 1.137625
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $Frechet_distance_metric
#>           Track_1   Track_2
#> Track_1        NA 0.7821257
#> Track_2 0.7821257        NA
#> 
#> $Frechet_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2      NA      NA
#> 
#> $Frechet_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2      NA      NA
#> 
#> $Frechet_metric_p_values_combined
#> [1] 0.3333333
#> 
#> $Frechet_distance_metric_simulations
#> $Frechet_distance_metric_simulations[[1]]
#>         Track_1   Track_2
#> Track_1      NA 0.8806958
#> Track_2      NA        NA
#> 
#> $Frechet_distance_metric_simulations[[2]]
#>         Track_1   Track_2
#> Track_1      NA 0.7820993
#> Track_2      NA        NA
#> 
#> $Frechet_distance_metric_simulations[[3]]
#>         Track_1  Track_2
#> Track_1      NA 1.137625
#> Track_2      NA       NA
#> 
#> 

# Example 2: Simulating tracks using the "Constrained" model and comparing Frechet distance
# in the PaluxyRiver dataset  using the "Centroid" superposition method
s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2, superposition = "Centroid")
#> 2025-12-27 20:56:01.655815 Iteration 1
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 1.555129
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-27 20:56:01.870507 Iteration 2
#>  
#> Frechet metric
#>         Track_1 Track_2
#> Track_1      NA 2.54824
#> Track_2      NA      NA
#> ------------------------------------
#> 2025-12-27 20:56:02.050792 Iteration 3
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 5.288696
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $Frechet_distance_metric
#>           Track_1   Track_2
#> Track_1        NA 0.6589325
#> Track_2 0.6589325        NA
#> 
#> $Frechet_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $Frechet_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $Frechet_metric_p_values_combined
#> [1] 0
#> 
#> $Frechet_distance_metric_simulations
#> $Frechet_distance_metric_simulations[[1]]
#>         Track_1  Track_2
#> Track_1      NA 1.555129
#> Track_2      NA       NA
#> 
#> $Frechet_distance_metric_simulations[[2]]
#>         Track_1 Track_2
#> Track_1      NA 2.54824
#> Track_2      NA      NA
#> 
#> $Frechet_distance_metric_simulations[[3]]
#>         Track_1  Track_2
#> Track_1      NA 5.288696
#> Track_2      NA       NA
#> 
#> 

# Example 3: Simulating tracks using the "Unconstrained" model and comparing Frechet distance
# in the PaluxyRiver dataset using the "Origin" superposition method
s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3, superposition = "Origin")
#> 2025-12-27 20:56:02.453848 Iteration 1
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 5.015942
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-27 20:56:02.495196 Iteration 2
#>  
#> Frechet metric
#>         Track_1 Track_2
#> Track_1      NA      -1
#> Track_2      NA      NA
#> ------------------------------------
#> 2025-12-27 20:56:02.543024 Iteration 3
#>  
#> Frechet metric
#>         Track_1  Track_2
#> Track_1      NA 17.64431
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $Frechet_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 1.144628
#> Track_2 1.144628       NA
#> 
#> $Frechet_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2      NA      NA
#> 
#> $Frechet_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2      NA      NA
#> 
#> $Frechet_metric_p_values_combined
#> [1] 0.3333333
#> 
#> $Frechet_distance_metric_simulations
#> $Frechet_distance_metric_simulations[[1]]
#>         Track_1  Track_2
#> Track_1      NA 5.015942
#> Track_2      NA       NA
#> 
#> $Frechet_distance_metric_simulations[[2]]
#>         Track_1 Track_2
#> Track_1      NA      -1
#> Track_2      NA      NA
#> 
#> $Frechet_distance_metric_simulations[[3]]
#>         Track_1  Track_2
#> Track_1      NA 17.64431
#> Track_2      NA       NA
#> 
#> 

```
