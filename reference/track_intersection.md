# Calculate intersection metrics in tracks

`track_intersection()` calculates the number of unique intersections
between trajectories. The function also supports testing with
simulations and different permutation procedures for the coordinates of
the simulated trajectories' origins to compute *p*-values. This allows
for a robust assessment of the intersection metrics, enabling users to
evaluate the significance of the observed intersections in relation to
simulated trajectories.

## Usage

``` r
track_intersection(
  data,
  test = NULL,
  H1 = NULL,
  sim = NULL,
  origin.permutation = NULL,
  custom.coord = NULL
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

- test:

  Logical; if `TRUE`, the function compares the observed intersection
  metrics against. Default is `FALSE`.

- H1:

  A character string specifying the alternative hypothesis to be tested.
  Options are `"Lower"` for testing whether the observed intersections
  are significantly lower than the simulated ones (e.g., coordinated or
  gregarious movement), or `"Higher"` for testing whether the observed
  intersections are significantly higher than the simulated ones (e.g.,
  predatory or chasing events).

- sim:

  A `track simulation` R object consisting of a list of simulated
  trajectories to use for comparison when `test = TRUE`.

- origin.permutation:

  A character string specifying the method for permutation of the
  coordinates of the simulated trajectories' origins. Options include
  `"None"`, `"Min.Box"`, `"Conv.Hull"`, or `"Custom"`. Default is
  `"None"`.

- custom.coord:

  A matrix of custom coordinates that define the vertices of an area for
  permutation of the coordinates of the simulated trajectories' origins.

## Value

A `track intersection` R object consisting of a list containing the
following elements:

- Intersection_metric:

  A numeric matrix of unique intersection counts between trajectories.
  Each entry represents the number of unique intersection points between
  the corresponding pair of trajectories.

- Intersection_metric_p_values:

  (If `test = TRUE`) A numeric matrix of raw pairwise *p*-values,
  computed by Monte Carlo tail tests according to `H1`, with the (+1)
  correction (Phipson & Smyth, 2010): \$\$p = (1 + \\\\\text{extreme}\\)
  / (nsim + 1)\$\$. Each entry reflects the probability of observing an
  intersection count as extreme as the observed one, given the null
  hypothesis of no difference.

- Intersection_metric_p_values_BH:

  (If `test = TRUE`) A numeric matrix of Benjamini–Hochberg (BH)
  adjusted *p*-values controlling the false discovery rate (FDR),
  applied across the set of unique pairwise tests (Benjamini & Hochberg,
  1995).

- Intersection_metric_p_values_combined:

  (If `test = TRUE`) A single numeric value representing the combined
  *p*-value across all intersections (based on the global statistic: the
  sum of pairwise intersections). This indicates the overall
  significance of the observed intersections relative to simulations.

- Intersection_metric_simulations:

  (If `test = TRUE`) A list containing matrices of intersection counts
  for each simulation iteration, allowing for inspection of the
  distribution of intersections across randomized scenarios.

## Details

The `track_intersection()` function is designed to calculate the number
of unique intersections between trajectories and to evaluate their
statistical significance through hypothesis testing based on simulated
tracks.

Pairwise *p*-values are computed with a Monte Carlo tail test according
to the selected alternative hypothesis (`H1`), using the (+1) correction
to avoid zero-values (see Phipson & Smyth, 2010): \$\$p = (1 +
\\\\\text{extreme}\\) / (nsim + 1)\$\$

These raw *p*-values are then adjusted for multiple comparisons across
the set of unique pairwise tests using the Benjamini–Hochberg (BH)
procedure for false discovery rate (FDR) control (Benjamini & Hochberg,
1995). Both the raw and the BH-adjusted *p*-value matrices are returned
in the output object, allowing users to inspect either uncorrected or
corrected results. In addition, a global combined *p*-value is provided,
summarizing the overall deviation from the null across all pairs.

This framework provides a robust means of comparing observed
intersections against those expected under random conditions, and is
intended for testing specific behavioral hypotheses related to
trackmaker interactions.

Hypothesis testing is controlled by the `H1` argument, which defines the
**alternative hypothesis**:

- `"Lower"`: tests whether the observed intersections are significantly
  fewer than those generated by simulations. This scenario corresponds
  to hypotheses of **coordinated or gregarious movement**, where
  parallel or group movement would result in fewer intersections than
  expected at random.

- `"Higher"`: tests whether the observed intersections are significantly
  greater than those generated by simulations. This scenario corresponds
  to hypotheses of **predatory or chasing interactions**, where one
  trackmaker follows or crosses another, producing more intersections
  than expected at random.

The choice of `H1` must match the biological hypothesis under
investigation. If it is not explicitly specified when `test = TRUE`, an
error will be raised.

The interpretation of the **combined global *p*-value** returned by the
function is directly influenced by the choice of `H1`, as it determines
whether the test detects a reduction or an increase in intersection
counts relative to the simulated dataset.

In addition to hypothesis testing, the `track_intersection()` function
offers several options for altering the initial positions of simulated
tracks through the `origin.permutation` argument. The available options
include:

- `"None"`: Simulated trajectories are not shifted and are compared
  based on their original starting positions.

- `"Min.Box"`: Trajectories are randomly placed within the **minimum
  bounding box** surrounding the original starting points.

- `"Conv.Hull"`: Trajectories are placed within the **convex hull** that
  encompasses all original starting points, providing a more precise
  representation of the area occupied by the tracks.

- `"Custom"`: Allows users to define a specific region of interest by
  providing a matrix of coordinates (`custom.coord`) that specifies the
  vertices of the desired area. This option is particularly useful when
  certain spatial features or environmental conditions are known to
  constrain movement.

The choice of `origin.permutation` should reflect the nature of the
behavioral hypothesis being tested. For example, using `"None"` is most
appropriate when testing how intersections compare under scenarios where
trackmakers originate from specific locations. In contrast, options like
`"Min.Box"`, `"Conv.Hull"`, or `"Custom"` are suitable when evaluating
how intersections would differ if the tracks originated from a broader
or predefined area.

The `track_intersection()` function also allows for integration with
similarity metrics computed using
[`simil_DTW_metric()`](https://macrofunuv.github.io/QuAnTeTrack/reference/simil_DTW_metric.md)
and
[`simil_Frechet_metric()`](https://macrofunuv.github.io/QuAnTeTrack/reference/simil_Frechet_metric.md).
This combination of intersection counts and similarity metrics can
provide a more comprehensive analysis of how trackmakers interacted,
whether their movements were coordinated or independent, and whether
their interactions were consistent with the hypothesized behavioral
patterns.

Overall, the selection of `H1` and `origin.permutation` should be
carefully considered in light of the specific hypotheses being tested.
By combining intersection metrics with similarity measures, users can
obtain a deeper understanding of the behavioral dynamics underlying the
observed trackways.

## Logo

![](figures/Logo.png)

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. Journal of
the Royal statistical society: series B (Methodological), 57(1),
289-300.

Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be
zero: calculating exact P-values when permutations are randomly drawn.
Statistical applications in genetics and molecular biology, 9(1).

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`simulate_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/simulate_track.md),
[`simil_DTW_metric`](https://macrofunuv.github.io/QuAnTeTrack/reference/simil_DTW_metric.md),
[`simil_Frechet_metric`](https://macrofunuv.github.io/QuAnTeTrack/reference/simil_Frechet_metric.md)

## Examples

``` r
# Example 1: Intersection metrics in the PaluxyRiver dataset.
s1 <- simulate_track(PaluxyRiver, nsim = 5, model = "Directed")
int1 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1,
origin.permutation = "None")
#> 2025-12-29 11:24:24.532472 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:24.544037 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:24.555635 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:24.5676 Iteration 4
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:24.579001 Iteration 5
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
print(int1)
#> $Intersection_metric
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2       1      NA
#> 
#> $Intersection_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA     0.8
#> Track_2      NA      NA
#> 
#> $Intersection_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.8
#> Track_2      NA      NA
#> 
#> $Intersection_metric_p_values_combined
#> [1] 0.2
#> 
#> $Intersection_metric_simulations
#> $Intersection_metric_simulations[[1]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[2]]
#>         Track_1 Track_2
#> Track_1      NA       3
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[3]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[4]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[5]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> 

# Example 2: Using "Min.Box" origin permutation in PaluxyRiver dataset.
s2 <- simulate_track(PaluxyRiver, nsim = 5, model = "Constrained")
int2 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2,
origin.permutation = "Min.Box")
#> 2025-12-29 11:24:24.876117 Permutation 1
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:25.308197 Permutation 2
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:25.477045 Permutation 3
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:25.591635 Permutation 4
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:25.730189 Permutation 5
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2025-12-29 11:24:25.742034 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:25.752998 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:25.763779 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:25.774715 Iteration 4
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:25.785782 Iteration 5
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
print(int2)
#> $Intersection_metric
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2       1      NA
#> 
#> $Intersection_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA     0.8
#> Track_2      NA      NA
#> 
#> $Intersection_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.8
#> Track_2      NA      NA
#> 
#> $Intersection_metric_p_values_combined
#> [1] 0.8
#> 
#> $Intersection_metric_simulations
#> $Intersection_metric_simulations[[1]]
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[2]]
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[3]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[4]]
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[5]]
#>         Track_1 Track_2
#> Track_1      NA       2
#> Track_2      NA      NA
#> 
#> 

# Example 3: Using "Conv.Hull" origin permutation in PaluxyRiver dataset.
s3 <- simulate_track(PaluxyRiver, nsim = 5, model = "Unconstrained")
int3 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3,
origin.permutation = "Conv.Hull")
#> 2025-12-29 11:24:26.179644 Permutation 1
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2025-12-29 11:24:26.327468 Permutation 2
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2025-12-29 11:24:26.436626 Permutation 3
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2025-12-29 11:24:26.545647 Permutation 4
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> 2025-12-29 11:24:26.665055 Permutation 5
#>  
#> Permutation of coordinates at origin using Conv.Hull
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2025-12-29 11:24:26.676775 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.687893 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.698816 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.709784 Iteration 4
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.720803 Iteration 5
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
print(int3)
#> $Intersection_metric
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2       1      NA
#> 
#> $Intersection_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA     0.8
#> Track_2      NA      NA
#> 
#> $Intersection_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.8
#> Track_2      NA      NA
#> 
#> $Intersection_metric_p_values_combined
#> [1] 0.4
#> 
#> $Intersection_metric_simulations
#> $Intersection_metric_simulations[[1]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[2]]
#>         Track_1 Track_2
#> Track_1      NA       2
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[3]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[4]]
#>         Track_1 Track_2
#> Track_1      NA       0
#> Track_2      NA      NA
#> 
#> $Intersection_metric_simulations[[5]]
#>         Track_1 Track_2
#> Track_1      NA       1
#> Track_2      NA      NA
#> 
#> 

# Example 4: Using "Min.Box" origin permutation in MountTom subset.
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
s4 <- simulate_track(sbMountTom, nsim = 5)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.
int4 <- track_intersection(sbMountTom, test = TRUE, H1 = "Higher", sim = s4,
origin.permutation = "Min.Box")
#> 2025-12-29 11:24:26.815009 Permutation 1
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:26.819485 Permutation 2
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:26.823824 Permutation 3
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:26.828039 Permutation 4
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> 2025-12-29 11:24:26.832407 Permutation 5
#>  
#> Permutation of coordinates at origin using Min.Box
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2025-12-29 11:24:26.859666 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.886473 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.914574 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.94186 Iteration 4
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:26.969162 Iteration 5
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
print(int4)
#> $Intersection_metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02        0       NA        0        0        0        0        0
#> Track_03        0        0       NA        0        0        0        0
#> Track_04        0        0        0       NA        0        0        0
#> Track_07        0        0        0        0       NA        0        1
#> Track_08        0        0        0        0        0       NA        0
#> Track_09        0        0        0        0        1        0       NA
#> Track_13        0        0        0        0        0        0        0
#> Track_15        0        0        0        0        0        0        0
#> Track_16        0        0        0        0        0        0        0
#> Track_18        0        0        0        0        0        0        0
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15        0       NA        0        0
#> Track_16        0        0       NA        0
#> Track_18        0        0        0       NA
#> 
#> $Intersection_metric_p_values
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        1        1        1        1        1        1
#> Track_02       NA       NA        1        1        1        1        1
#> Track_03       NA       NA       NA        1        1        1        1
#> Track_04       NA       NA       NA       NA        1        1        1
#> Track_07       NA       NA       NA       NA       NA        1        0
#> Track_08       NA       NA       NA       NA       NA       NA        1
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        1        1        1        1
#> Track_02        1        1        1        1
#> Track_03        1        1        1        1
#> Track_04        1        1        1        1
#> Track_07        1        1        1        1
#> Track_08        1        1        1        1
#> Track_09        1        1        1        1
#> Track_13       NA        1        1        1
#> Track_15       NA       NA        1        1
#> Track_16       NA       NA       NA        1
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_p_values_BH
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        1        1        1        1        1        1
#> Track_02       NA       NA        1        1        1        1        1
#> Track_03       NA       NA       NA        1        1        1        1
#> Track_04       NA       NA       NA       NA        1        1        1
#> Track_07       NA       NA       NA       NA       NA        1        0
#> Track_08       NA       NA       NA       NA       NA       NA        1
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        1        1        1        1
#> Track_02        1        1        1        1
#> Track_03        1        1        1        1
#> Track_04        1        1        1        1
#> Track_07        1        1        1        1
#> Track_08        1        1        1        1
#> Track_09        1        1        1        1
#> Track_13       NA        1        1        1
#> Track_15       NA       NA        1        1
#> Track_16       NA       NA       NA        1
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_p_values_combined
#> [1] 0.6
#> 
#> $Intersection_metric_simulations
#> $Intersection_metric_simulations[[1]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        1        0        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        1
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[2]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[3]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        1        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[4]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        1
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[5]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        0        1        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        1
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        1        0
#> Track_09        0        0        1        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> 

# Example 5: Customized origin permutation in MountTom subset.
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
s5 <- simulate_track(sbMountTom, nsim = 5)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.
area_origin <- matrix(c(50, 5, 10, 5, 10, 20, 50, 20), ncol = 2, byrow = TRUE)
int5 <- track_intersection(sbMountTom, test = TRUE, H1 = "Higher", sim = s5,
origin.permutation = "Custom", custom.coord = area_origin)
#> 2025-12-29 11:24:27.075171 Permutation 1
#>  
#> Permutation of coordinates at origin using Custom
#> ------------------------------------
#> 2025-12-29 11:24:27.078667 Permutation 2
#>  
#> Permutation of coordinates at origin using Custom
#> ------------------------------------
#> 2025-12-29 11:24:27.082028 Permutation 3
#>  
#> Permutation of coordinates at origin using Custom
#> ------------------------------------
#> 2025-12-29 11:24:27.085354 Permutation 4
#>  
#> Permutation of coordinates at origin using Custom
#> ------------------------------------
#> 2025-12-29 11:24:27.088615 Permutation 5
#>  
#> Permutation of coordinates at origin using Custom
#> ------------------------------------
#> PERMUTATION COMPLETED
#> ------------------------------------
#>  
#> 2025-12-29 11:24:27.11556 Iteration 1
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:27.142146 Iteration 2
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:27.169033 Iteration 3
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:27.196093 Iteration 4
#>  
#> Intersect metric
#> ------------------------------------
#> 2025-12-29 11:24:27.223728 Iteration 5
#>  
#> Intersect metric
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
print(int5)
#> $Intersection_metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02        0       NA        0        0        0        0        0
#> Track_03        0        0       NA        0        0        0        0
#> Track_04        0        0        0       NA        0        0        0
#> Track_07        0        0        0        0       NA        0        1
#> Track_08        0        0        0        0        0       NA        0
#> Track_09        0        0        0        0        1        0       NA
#> Track_13        0        0        0        0        0        0        0
#> Track_15        0        0        0        0        0        0        0
#> Track_16        0        0        0        0        0        0        0
#> Track_18        0        0        0        0        0        0        0
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15        0       NA        0        0
#> Track_16        0        0       NA        0
#> Track_18        0        0        0       NA
#> 
#> $Intersection_metric_p_values
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        1        1        1        1        1        1
#> Track_02       NA       NA        1        1        1        1        1
#> Track_03       NA       NA       NA        1        1        1        1
#> Track_04       NA       NA       NA       NA        1        1        1
#> Track_07       NA       NA       NA       NA       NA        1        0
#> Track_08       NA       NA       NA       NA       NA       NA        1
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        1        1        1        1
#> Track_02        1        1        1        1
#> Track_03        1        1        1        1
#> Track_04        1        1        1        1
#> Track_07        1        1        1        1
#> Track_08        1        1        1        1
#> Track_09        1        1        1        1
#> Track_13       NA        1        1        1
#> Track_15       NA       NA        1        1
#> Track_16       NA       NA       NA        1
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_p_values_BH
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        1        1        1        1        1        1
#> Track_02       NA       NA        1        1        1        1        1
#> Track_03       NA       NA       NA        1        1        1        1
#> Track_04       NA       NA       NA       NA        1        1        1
#> Track_07       NA       NA       NA       NA       NA        1        0
#> Track_08       NA       NA       NA       NA       NA       NA        1
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        1        1        1        1
#> Track_02        1        1        1        1
#> Track_03        1        1        1        1
#> Track_04        1        1        1        1
#> Track_07        1        1        1        1
#> Track_08        1        1        1        1
#> Track_09        1        1        1        1
#> Track_13       NA        1        1        1
#> Track_15       NA       NA        1        1
#> Track_16       NA       NA       NA        1
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_p_values_combined
#> [1] 0.4
#> 
#> $Intersection_metric_simulations
#> $Intersection_metric_simulations[[1]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        1        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        1        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[2]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[3]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[4]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        0        0        0        0        0
#> Track_02       NA       NA        0        0        1        0        0
#> Track_03       NA       NA       NA        0        0        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        1        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> $Intersection_metric_simulations[[5]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA        0        1        0        1        0        0
#> Track_02       NA       NA        0        0        0        0        0
#> Track_03       NA       NA       NA        0        1        0        0
#> Track_04       NA       NA       NA       NA        0        0        0
#> Track_07       NA       NA       NA       NA       NA        0        0
#> Track_08       NA       NA       NA       NA       NA       NA        0
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01        0        0        0        0
#> Track_02        0        0        0        0
#> Track_03        0        0        0        0
#> Track_04        0        0        0        0
#> Track_07        0        0        0        0
#> Track_08        0        0        0        0
#> Track_09        0        0        0        0
#> Track_13       NA        0        0        0
#> Track_15       NA       NA        0        0
#> Track_16       NA       NA       NA        0
#> Track_18       NA       NA       NA       NA
#> 
#> 
```
