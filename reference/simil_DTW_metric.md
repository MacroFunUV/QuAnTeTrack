# Similarity metric using Dynamic Time Warping (DTW)

`simil_DTW_metric()` computes similarity metrics between two or more
trajectories using Dynamic Time Warping (DTW). It allows for different
superposition methods to align trajectories before calculating the DTW
metric. The function also supports testing with simulations to calculate
*p*-values for the DTW distance metrics.

## Usage

``` r
simil_DTW_metric(data, test = FALSE, sim = NULL, superposition = "None")
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

  Logical; if `TRUE`, the function compares the observed DTW distances
  against simulated trajectories and calculates *p*-values. Default is
  `FALSE`.

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

- DTW_distance_metric:

  A numeric matrix of pairwise Dynamic Time Warping (DTW) distances
  between trajectories. Each entry represents the DTW distance between
  the corresponding pair of trajectories.

- DTW_distance_metric_p_values:

  (If `test = TRUE`) A numeric matrix of raw pairwise *p*-values,
  computed by Monte Carlo tail tests with the (+1) correction (Phipson &
  Smyth, 2010): \$\$p = (1 + \\\\\text{extreme}\\) / (nsim + 1)\$\$.
  Each entry reflects the probability of observing a DTW distance as
  extreme as the observed one, given the null hypothesis of no
  difference.

- DTW_distance_metric_p_values_BH:

  (If `test = TRUE`) A numeric matrix of Benjamini–Hochberg (BH)
  adjusted *p*-values controlling the false discovery rate (FDR),
  applied across the set of unique pairwise tests (Benjamini & Hochberg,
  1995).

- DTW_metric_p_values_combined:

  (If `test = TRUE`) A single numeric value representing the combined
  *p*-value across all DTW distances (based on the global statistic: the
  sum of pairwise distances). This indicates the overall significance of
  the observed DTW distances relative to simulations.

- DTW_distance_metric_simulations:

  (If `test = TRUE`) A list containing matrices of DTW distances for
  each simulation iteration, allowing for inspection of the distribution
  of DTW distances across randomized scenarios.

## Details

The `simil_DTW_metric()` function calculates the similarity between
trajectories using the Dynamic Time Warping (DTW) algorithm from the dtw
package. The `dtw()` function is used with the `dist.method` argument
set to `"Euclidean"` for computing the local distances between points in
the trajectories.

DTW aligns two time series by minimizing the cumulative distance between
their points, creating an optimal alignment despite variations in length
or temporal distortions. The algorithm constructs a distance matrix
where each element represents the cost of aligning points between the
two series and finds a warping path through this matrix that minimizes
the total distance. The warping path is contiguous and monotonic,
starting from the bottom-left corner and ending at the top-right corner
(Cleasby et al., 2019).

DTW measures are non-negative and unbounded, with larger values
indicating greater dissimilarity between the time series. This method
has been used in various contexts, including ecological studies to
analyze and cluster trajectory data (Cleasby et al., 2019).

Potential limitations and biases of DTW include sensitivity to noise and
outliers, computational complexity, and the need for appropriate
distance metrics. Additionally, DTW may not always account for all
structural differences between trajectories and can be biased by the
chosen alignment constraints. While DTW can handle trajectories of
different lengths due to its elastic nature, having trajectories of
similar lengths can improve the accuracy and interpretability of the
similarity measure. Similar lengths result in a more meaningful
alignment and can make the computation more efficient. When trajectories
differ significantly in length, preprocessing or normalization might be
necessary, and careful analysis is required to understand the alignment
path. The function’s flexibility in handling different lengths allows it
to be applied in various contexts. However, large differences in
trajectory lengths might introduce potential biases that should be
considered when interpreting the results.

The function offers three different superposition methods to align the
trajectories before `DTW()` is applied:

- `"None"`: No superposition is applied.

- `"Centroid"`: Trajectories are shifted to align based on their
  centroids.

- `"Origin"`: Trajectories are shifted to align based on their starting
  point.

If `test = TRUE`, the function can compute *p*-values by comparing the
observed DTW distances with those generated from a set of simulated
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
[`dtw`](https://rdrr.io/pkg/dtw/man/dtw.html)

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
# Example 1: Simulating tracks using the "Directed" model and comparing DTW distance
# in the PaluxyRiver dataset
s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#> 2025-12-29 23:15:38.441995 Iteration 1
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 17.58427
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.444784 Iteration 2
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 16.48403
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.447416 Iteration 3
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 24.70844
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 23.47722
#> Track_2 23.47722       NA
#> 
#> $DTW_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA    0.75
#> Track_2      NA      NA
#> 
#> $DTW_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.75
#> Track_2      NA      NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0.6666667
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>         Track_1  Track_2
#> Track_1      NA 17.58427
#> Track_2      NA       NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>         Track_1  Track_2
#> Track_1      NA 16.48403
#> Track_2      NA       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>         Track_1  Track_2
#> Track_1      NA 24.70844
#> Track_2      NA       NA
#> 
#> 

# Example 2: Simulating tracks using the "Constrained" model and comparing DTW distance
# in the PaluxyRiver dataset
s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2, superposition = "None")
#> 2025-12-29 23:15:38.463222 Iteration 1
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 252.5979
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.466031 Iteration 2
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 95.39256
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.468695 Iteration 3
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 170.2938
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 23.47722
#> Track_2 23.47722       NA
#> 
#> $DTW_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $DTW_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>         Track_1  Track_2
#> Track_1      NA 252.5979
#> Track_2      NA       NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>         Track_1  Track_2
#> Track_1      NA 95.39256
#> Track_2      NA       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>         Track_1  Track_2
#> Track_1      NA 170.2938
#> Track_2      NA       NA
#> 
#> 

# Example 3: Simulating tracks using the "Unconstrained" model and comparing DTW distance
# in the PaluxyRiver dataset
s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3, superposition = "None")
#> 2025-12-29 23:15:38.484208 Iteration 1
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 717.4026
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.48689 Iteration 2
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 738.6735
#> Track_2      NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.489525 Iteration 3
#>  
#> DTW metric
#>         Track_1  Track_2
#> Track_1      NA 766.5934
#> Track_2      NA       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 23.47722
#> Track_2 23.47722       NA
#> 
#> $DTW_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $DTW_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2      NA      NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>         Track_1  Track_2
#> Track_1      NA 717.4026
#> Track_2      NA       NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>         Track_1  Track_2
#> Track_1      NA 738.6735
#> Track_2      NA       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>         Track_1  Track_2
#> Track_1      NA 766.5934
#> Track_2      NA       NA
#> 
#> 

# Example 4: Simulating and comparing DTW distance in the MountTom dataset using the
# "Centroid" superposition method
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
s4 <- simulate_track(sbMountTom, nsim = 3)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.
simil_DTW_metric(sbMountTom, test = TRUE, sim = s4, superposition = "Centroid")
#> 2025-12-29 23:15:38.626829 Iteration 1
#>  
#> DTW metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08 Track_09
#> Track_01       NA 34.77507 40.83665 28.97181  7.98321 39.801883 40.97504
#> Track_02       NA       NA 33.15206 52.11918 20.13191 25.407164 18.04472
#> Track_03       NA       NA       NA 37.95267 23.53295  8.399035 19.73591
#> Track_04       NA       NA       NA       NA 24.44538 40.943513 49.13389
#> Track_07       NA       NA       NA       NA       NA 21.926868 23.58801
#> Track_08       NA       NA       NA       NA       NA        NA 11.16308
#> Track_09       NA       NA       NA       NA       NA        NA       NA
#> Track_13       NA       NA       NA       NA       NA        NA       NA
#> Track_15       NA       NA       NA       NA       NA        NA       NA
#> Track_16       NA       NA       NA       NA       NA        NA       NA
#> Track_18       NA       NA       NA       NA       NA        NA       NA
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 39.12781 26.669164 15.953157 39.613209
#> Track_02 40.79062 36.246694 21.116625 34.865976
#> Track_03 14.93785 20.189044 12.790502  5.126975
#> Track_04 25.34226  9.750902 15.373086 33.302578
#> Track_07 24.07498 16.958717  8.490390 22.822048
#> Track_08 20.46661 23.313536 12.973165 12.013410
#> Track_09 31.60224 30.167942 17.917873 23.002248
#> Track_13       NA 12.313347 12.653304 10.495460
#> Track_15       NA        NA  9.390416 16.660829
#> Track_16       NA        NA        NA 11.736856
#> Track_18       NA        NA        NA        NA
#> ------------------------------------
#> 2025-12-29 23:15:38.673233 Iteration 2
#>  
#> DTW metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08  Track_09
#> Track_01       NA 12.45845 39.22490 66.70123 24.60561  8.035231 24.928465
#> Track_02       NA       NA 35.58791 59.53738 17.91867 12.662517 16.324110
#> Track_03       NA       NA       NA 12.72239 18.52695 20.454961 24.738623
#> Track_04       NA       NA       NA       NA 33.03691 42.737136 39.977566
#> Track_07       NA       NA       NA       NA       NA 14.446761  4.820698
#> Track_08       NA       NA       NA       NA       NA        NA 17.287810
#> Track_09       NA       NA       NA       NA       NA        NA        NA
#> Track_13       NA       NA       NA       NA       NA        NA        NA
#> Track_15       NA       NA       NA       NA       NA        NA        NA
#> Track_16       NA       NA       NA       NA       NA        NA        NA
#> Track_18       NA       NA       NA       NA       NA        NA        NA
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 31.623826 42.359709 28.294056 23.937678
#> Track_02 22.775928 38.565890 25.353980 17.303799
#> Track_03 21.854177  4.311715  5.184952 17.583788
#> Track_04 34.014885 14.593099 19.033953 32.036203
#> Track_07  6.162870 22.039899 12.739054  1.269285
#> Track_08 20.025424 23.749981 12.542857 13.868949
#> Track_09  7.600545 27.684611 17.961896  4.639764
#> Track_13        NA 25.246732 16.600928  6.335406
#> Track_15        NA        NA  7.889243 21.107408
#> Track_16        NA        NA        NA 12.047081
#> Track_18        NA        NA        NA        NA
#> ------------------------------------
#> 2025-12-29 23:15:38.719369 Iteration 3
#>  
#> DTW metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08  Track_09
#> Track_01       NA 19.51285 38.94657 60.13078 40.02200 36.102832 35.892140
#> Track_02       NA       NA 33.11222 45.74144 33.69174 27.286147 23.896833
#> Track_03       NA       NA       NA 17.33850  3.67163  5.953785 15.294390
#> Track_04       NA       NA       NA       NA 13.97407  8.426659 14.435758
#> Track_07       NA       NA       NA       NA       NA  5.192675 14.030102
#> Track_08       NA       NA       NA       NA       NA        NA  8.529183
#> Track_09       NA       NA       NA       NA       NA        NA        NA
#> Track_13       NA       NA       NA       NA       NA        NA        NA
#> Track_15       NA       NA       NA       NA       NA        NA        NA
#> Track_16       NA       NA       NA       NA       NA        NA        NA
#> Track_18       NA       NA       NA       NA       NA        NA        NA
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01  5.477892 11.894416 29.161582 39.077199
#> Track_02 12.051650  6.148277 25.281876 35.845563
#> Track_03 26.606883 22.020042  4.558738  7.450284
#> Track_04 46.078995 37.863686 20.287013 28.371299
#> Track_07 28.137917 23.189137  7.140877 10.253712
#> Track_08 23.745531 18.512185  6.477767 11.670501
#> Track_09 24.559188 18.657639 12.677109 21.094087
#> Track_13        NA  7.365335 18.682095 27.305054
#> Track_15        NA        NA 15.022130 23.319641
#> Track_16        NA        NA        NA  5.848749
#> Track_18        NA        NA        NA        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08 Track_09
#> Track_01        NA  7.155827  6.889059  7.146449  8.550492  9.358030 41.27125
#> Track_02  7.155827        NA  7.367287  7.095034  8.511176  8.404712 44.46323
#> Track_03  6.889059  7.367287        NA  2.752378  5.017897  5.558007 23.61988
#> Track_04  7.146449  7.095034  2.752378        NA  3.748965  5.898481 21.55969
#> Track_07  8.550492  8.511176  5.017897  3.748965        NA  3.362680 27.87452
#> Track_08  9.358030  8.404712  5.558007  5.898481  3.362680        NA 27.11038
#> Track_09 41.271248 44.463229 23.619876 21.559693 27.874517 27.110379       NA
#> Track_13 23.690795 26.769160 11.380548 10.754453 16.245880 16.011649 12.44577
#> Track_15 10.006241  8.652232  4.011654  2.107453  4.076379  5.869619 21.45921
#> Track_16  8.956456  8.004512  6.187040  7.121107  5.378815  3.814577 36.64336
#> Track_18 11.561263  9.939214  4.835886  3.128265  5.173189  6.169448 21.46042
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 23.69080 10.006241  8.956456 11.561263
#> Track_02 26.76916  8.652232  8.004512  9.939214
#> Track_03 11.38055  4.011654  6.187040  4.835886
#> Track_04 10.75445  2.107453  7.121107  3.128265
#> Track_07 16.24588  4.076379  5.378815  5.173189
#> Track_08 16.01165  5.869619  3.814577  6.169448
#> Track_09 12.44577 21.459210 36.643358 21.460416
#> Track_13       NA 11.406049 22.371414 11.843307
#> Track_15 11.40605        NA  7.547234  1.280341
#> Track_16 22.37141  7.547234        NA  8.182026
#> Track_18 11.84331  1.280341  8.182026        NA
#> 
#> $DTW_distance_metric_p_values
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA     0.25     0.25     0.25     0.50     0.50     1.00
#> Track_02       NA       NA     0.25     0.25     0.25     0.25     1.00
#> Track_03       NA       NA       NA     0.25     0.50     0.25     0.75
#> Track_04       NA       NA       NA       NA     0.25     0.25     0.50
#> Track_07       NA       NA       NA       NA       NA     0.25     1.00
#> Track_08       NA       NA       NA       NA       NA       NA     1.00
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01     0.50     0.25     0.25     0.25
#> Track_02     0.75     0.50     0.25     0.25
#> Track_03     0.25     0.25     0.75     0.25
#> Track_04     0.25     0.25     0.25     0.25
#> Track_07     0.50     0.25     0.25     0.50
#> Track_08     0.25     0.25     0.25     0.25
#> Track_09     0.50     0.50     1.00     0.75
#> Track_13       NA     0.50     1.00     0.75
#> Track_15       NA       NA     0.25     0.25
#> Track_16       NA       NA       NA     0.50
#> Track_18       NA       NA       NA       NA
#> 
#> $DTW_distance_metric_p_values_BH
#>          Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01       NA 0.4296875 0.4296875 0.4296875 0.6250000 0.6250000 1.0000000
#> Track_02       NA        NA 0.4296875 0.4296875 0.4296875 0.4296875 1.0000000
#> Track_03       NA        NA        NA 0.4296875 0.6250000 0.4296875 0.8418367
#> Track_04       NA        NA        NA        NA 0.4296875 0.4296875 0.6250000
#> Track_07       NA        NA        NA        NA        NA 0.4296875 1.0000000
#> Track_08       NA        NA        NA        NA        NA        NA 1.0000000
#> Track_09       NA        NA        NA        NA        NA        NA        NA
#> Track_13       NA        NA        NA        NA        NA        NA        NA
#> Track_15       NA        NA        NA        NA        NA        NA        NA
#> Track_16       NA        NA        NA        NA        NA        NA        NA
#> Track_18       NA        NA        NA        NA        NA        NA        NA
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 0.6250000 0.4296875 0.4296875 0.4296875
#> Track_02 0.8418367 0.6250000 0.4296875 0.4296875
#> Track_03 0.4296875 0.4296875 0.8418367 0.4296875
#> Track_04 0.4296875 0.4296875 0.4296875 0.4296875
#> Track_07 0.6250000 0.4296875 0.4296875 0.6250000
#> Track_08 0.4296875 0.4296875 0.4296875 0.4296875
#> Track_09 0.6250000 0.6250000 1.0000000 0.8418367
#> Track_13        NA 0.6250000 1.0000000 0.8418367
#> Track_15        NA        NA 0.4296875 0.4296875
#> Track_16        NA        NA        NA 0.6250000
#> Track_18        NA        NA        NA        NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08 Track_09
#> Track_01       NA 34.77507 40.83665 28.97181  7.98321 39.801883 40.97504
#> Track_02       NA       NA 33.15206 52.11918 20.13191 25.407164 18.04472
#> Track_03       NA       NA       NA 37.95267 23.53295  8.399035 19.73591
#> Track_04       NA       NA       NA       NA 24.44538 40.943513 49.13389
#> Track_07       NA       NA       NA       NA       NA 21.926868 23.58801
#> Track_08       NA       NA       NA       NA       NA        NA 11.16308
#> Track_09       NA       NA       NA       NA       NA        NA       NA
#> Track_13       NA       NA       NA       NA       NA        NA       NA
#> Track_15       NA       NA       NA       NA       NA        NA       NA
#> Track_16       NA       NA       NA       NA       NA        NA       NA
#> Track_18       NA       NA       NA       NA       NA        NA       NA
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 39.12781 26.669164 15.953157 39.613209
#> Track_02 40.79062 36.246694 21.116625 34.865976
#> Track_03 14.93785 20.189044 12.790502  5.126975
#> Track_04 25.34226  9.750902 15.373086 33.302578
#> Track_07 24.07498 16.958717  8.490390 22.822048
#> Track_08 20.46661 23.313536 12.973165 12.013410
#> Track_09 31.60224 30.167942 17.917873 23.002248
#> Track_13       NA 12.313347 12.653304 10.495460
#> Track_15       NA        NA  9.390416 16.660829
#> Track_16       NA        NA        NA 11.736856
#> Track_18       NA        NA        NA        NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08  Track_09
#> Track_01       NA 12.45845 39.22490 66.70123 24.60561  8.035231 24.928465
#> Track_02       NA       NA 35.58791 59.53738 17.91867 12.662517 16.324110
#> Track_03       NA       NA       NA 12.72239 18.52695 20.454961 24.738623
#> Track_04       NA       NA       NA       NA 33.03691 42.737136 39.977566
#> Track_07       NA       NA       NA       NA       NA 14.446761  4.820698
#> Track_08       NA       NA       NA       NA       NA        NA 17.287810
#> Track_09       NA       NA       NA       NA       NA        NA        NA
#> Track_13       NA       NA       NA       NA       NA        NA        NA
#> Track_15       NA       NA       NA       NA       NA        NA        NA
#> Track_16       NA       NA       NA       NA       NA        NA        NA
#> Track_18       NA       NA       NA       NA       NA        NA        NA
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 31.623826 42.359709 28.294056 23.937678
#> Track_02 22.775928 38.565890 25.353980 17.303799
#> Track_03 21.854177  4.311715  5.184952 17.583788
#> Track_04 34.014885 14.593099 19.033953 32.036203
#> Track_07  6.162870 22.039899 12.739054  1.269285
#> Track_08 20.025424 23.749981 12.542857 13.868949
#> Track_09  7.600545 27.684611 17.961896  4.639764
#> Track_13        NA 25.246732 16.600928  6.335406
#> Track_15        NA        NA  7.889243 21.107408
#> Track_16        NA        NA        NA 12.047081
#> Track_18        NA        NA        NA        NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08  Track_09
#> Track_01       NA 19.51285 38.94657 60.13078 40.02200 36.102832 35.892140
#> Track_02       NA       NA 33.11222 45.74144 33.69174 27.286147 23.896833
#> Track_03       NA       NA       NA 17.33850  3.67163  5.953785 15.294390
#> Track_04       NA       NA       NA       NA 13.97407  8.426659 14.435758
#> Track_07       NA       NA       NA       NA       NA  5.192675 14.030102
#> Track_08       NA       NA       NA       NA       NA        NA  8.529183
#> Track_09       NA       NA       NA       NA       NA        NA        NA
#> Track_13       NA       NA       NA       NA       NA        NA        NA
#> Track_15       NA       NA       NA       NA       NA        NA        NA
#> Track_16       NA       NA       NA       NA       NA        NA        NA
#> Track_18       NA       NA       NA       NA       NA        NA        NA
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01  5.477892 11.894416 29.161582 39.077199
#> Track_02 12.051650  6.148277 25.281876 35.845563
#> Track_03 26.606883 22.020042  4.558738  7.450284
#> Track_04 46.078995 37.863686 20.287013 28.371299
#> Track_07 28.137917 23.189137  7.140877 10.253712
#> Track_08 23.745531 18.512185  6.477767 11.670501
#> Track_09 24.559188 18.657639 12.677109 21.094087
#> Track_13        NA  7.365335 18.682095 27.305054
#> Track_15        NA        NA 15.022130 23.319641
#> Track_16        NA        NA        NA  5.848749
#> Track_18        NA        NA        NA        NA
#> 
#> 

# Example 5: Simulating and comparing DTW distance in the MountTom dataset using the
# "Origin" superposition method
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
s5 <- simulate_track(sbMountTom, nsim = 3)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.
simil_DTW_metric(sbMountTom, test = TRUE, sim = s5, superposition = "Origin")
#> 2025-12-29 23:15:38.86796 Iteration 1
#>  
#> DTW metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08  Track_09
#> Track_01       NA 74.58195 78.91397 96.64398 70.79938 42.08277 35.468019
#> Track_02       NA       NA 63.99860 15.82293 23.15382 63.20830 70.717161
#> Track_03       NA       NA       NA 71.73955 34.28973 31.76544 44.275450
#> Track_04       NA       NA       NA       NA 18.85089 78.35770 87.732235
#> Track_07       NA       NA       NA       NA       NA 42.76890 50.107566
#> Track_08       NA       NA       NA       NA       NA       NA  9.465876
#> Track_09       NA       NA       NA       NA       NA       NA        NA
#> Track_13       NA       NA       NA       NA       NA       NA        NA
#> Track_15       NA       NA       NA       NA       NA       NA        NA
#> Track_16       NA       NA       NA       NA       NA       NA        NA
#> Track_18       NA       NA       NA       NA       NA       NA        NA
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 93.33423 37.413138 52.051414 44.623207
#> Track_02 65.57425 64.313729 27.180051 64.332909
#> Track_03 14.28249 36.947737 17.754758 30.247222
#> Track_04 69.93651 80.012922 30.825070 79.119289
#> Track_07 32.87243 46.118442  9.575422 42.182148
#> Track_08 46.37700  5.752637 23.576983  2.286923
#> Track_09 58.16662  5.798243 33.668082 12.004624
#> Track_13       NA 51.806516 18.002230 44.967019
#> Track_15       NA        NA 25.658456  7.668531
#> Track_16       NA        NA        NA 23.322257
#> Track_18       NA        NA        NA        NA
#> ------------------------------------
#> 2025-12-29 23:15:38.913939 Iteration 2
#>  
#> DTW metric
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08 Track_09
#> Track_01       NA 54.23687 45.21874 49.44016 71.50794 66.200748 50.16937
#> Track_02       NA       NA  7.59947 12.07591 40.92355 35.146152 75.64967
#> Track_03       NA       NA       NA 17.30699 24.03585 19.747136 53.78668
#> Track_04       NA       NA       NA       NA 54.70680 48.159534 81.22320
#> Track_07       NA       NA       NA       NA       NA  4.641382 57.88723
#> Track_08       NA       NA       NA       NA       NA        NA 56.10677
#> Track_09       NA       NA       NA       NA       NA        NA       NA
#> Track_13       NA       NA       NA       NA       NA        NA       NA
#> Track_15       NA       NA       NA       NA       NA        NA       NA
#> Track_16       NA       NA       NA       NA       NA        NA       NA
#> Track_18       NA       NA       NA       NA       NA        NA       NA
#>          Track_13  Track_15 Track_16 Track_18
#> Track_01 70.35715 38.246414 53.14282 26.44397
#> Track_02 79.04307 59.444327 37.94523 16.21459
#> Track_03 54.32899 40.387821 20.39138 12.63175
#> Track_04 89.50376 65.038178 48.86099 14.16531
#> Track_07 52.78516 47.239386 11.03906 31.14017
#> Track_08 49.82184 42.374859 11.26962 27.12250
#> Track_09 23.79028  5.888145 35.99928 42.43423
#> Track_13       NA 22.474039 27.90860 47.60924
#> Track_15       NA        NA 25.48473 30.77810
#> Track_16       NA        NA       NA 21.73607
#> Track_18       NA        NA       NA       NA
#> ------------------------------------
#> 2025-12-29 23:15:38.965872 Iteration 3
#>  
#> DTW metric
#>          Track_01 Track_02  Track_03  Track_04 Track_07 Track_08 Track_09
#> Track_01       NA 113.7108 75.852344 122.69581 77.46999 69.19479 43.15262
#> Track_02       NA       NA  8.193052  35.02683 19.17839 44.74527 69.97435
#> Track_03       NA       NA        NA  24.76149 14.27356 28.90253 43.40454
#> Track_04       NA       NA        NA        NA 47.01337 71.48320 65.05010
#> Track_07       NA       NA        NA        NA       NA 19.83662 51.52343
#> Track_08       NA       NA        NA        NA       NA       NA 52.01348
#> Track_09       NA       NA        NA        NA       NA       NA       NA
#> Track_13       NA       NA        NA        NA       NA       NA       NA
#> Track_15       NA       NA        NA        NA       NA       NA       NA
#> Track_16       NA       NA        NA        NA       NA       NA       NA
#> Track_18       NA       NA        NA        NA       NA       NA       NA
#>           Track_13  Track_15 Track_16  Track_18
#> Track_01 89.950716 75.784943 37.31391 75.446867
#> Track_02 26.529842 19.578187 42.36166 23.277897
#> Track_03 20.553526  9.494872 22.86946 16.753630
#> Track_04 58.201258 11.676309 58.04265 50.430771
#> Track_07  5.254308 23.198154 21.30381  5.295366
#> Track_08 21.428358 36.170033 13.73617 15.205568
#> Track_09 62.329274 39.511223 31.33146 50.831737
#> Track_13        NA 31.105605 26.59231  3.924337
#> Track_15        NA        NA 25.54778 24.486652
#> Track_16        NA        NA       NA 18.784094
#> Track_18        NA        NA       NA        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08 Track_09
#> Track_01        NA  8.365556  8.140001 13.288601 17.963246 14.043056 80.66455
#> Track_02  8.365556        NA  9.928547 13.870318 15.997694 11.986057 83.72110
#> Track_03  8.140001  9.928547        NA  4.242441  9.250629  6.131319 46.49768
#> Track_04 13.288601 13.870318  4.242441        NA  5.401950  3.172849 42.61717
#> Track_07 17.963246 15.997694  9.250629  5.401950        NA  5.261381 51.73148
#> Track_08 14.043056 11.986057  6.131319  3.172849  5.261381        NA 54.42241
#> Track_09 80.664547 83.721098 46.497678 42.617174 51.731485 54.422412       NA
#> Track_13 45.007884 47.197964 20.295656 19.901140 29.449500 27.328666 22.41087
#> Track_15 18.730370 17.799527  7.910351  4.165947  2.944420  4.915326 41.82503
#> Track_16 14.932046  9.836672 10.158387  8.969975  6.952380  5.349900 68.97341
#> Track_18 22.179224 20.569492  9.928217  6.396727  5.167036  6.863688 42.96907
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 45.00788 18.730370 14.932046 22.179224
#> Track_02 47.19796 17.799527  9.836672 20.569492
#> Track_03 20.29566  7.910351 10.158387  9.928217
#> Track_04 19.90114  4.165947  8.969975  6.396727
#> Track_07 29.44950  2.944420  6.952380  5.167036
#> Track_08 27.32867  4.915326  5.349900  6.863688
#> Track_09 22.41087 41.825031 68.973412 42.969072
#> Track_13       NA 20.579060 39.239001 21.761749
#> Track_15 20.57906        NA  8.572483  2.784293
#> Track_16 39.23900  8.572483        NA 10.795457
#> Track_18 21.76175  2.784293 10.795457        NA
#> 
#> $DTW_distance_metric_p_values
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA     0.25     0.25     0.25     0.25     0.25     1.00
#> Track_02       NA       NA     0.75     0.50     0.25     0.25     1.00
#> Track_03       NA       NA       NA     0.25     0.25     0.25     0.75
#> Track_04       NA       NA       NA       NA     0.25     0.25     0.25
#> Track_07       NA       NA       NA       NA       NA     0.50     0.75
#> Track_08       NA       NA       NA       NA       NA       NA     0.75
#> Track_09       NA       NA       NA       NA       NA       NA       NA
#> Track_13       NA       NA       NA       NA       NA       NA       NA
#> Track_15       NA       NA       NA       NA       NA       NA       NA
#> Track_16       NA       NA       NA       NA       NA       NA       NA
#> Track_18       NA       NA       NA       NA       NA       NA       NA
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01     0.25     0.25     0.25     0.25
#> Track_02     0.50     0.25     0.25     0.50
#> Track_03     0.50     0.25     0.25     0.25
#> Track_04     0.25     0.25     0.25     0.25
#> Track_07     0.50     0.25     0.25     0.25
#> Track_08     0.50     0.25     0.25     0.50
#> Track_09     0.25     1.00     1.00     0.75
#> Track_13       NA     0.25     1.00     0.50
#> Track_15       NA       NA     0.25     0.25
#> Track_16       NA       NA       NA     0.25
#> Track_18       NA       NA       NA       NA
#> 
#> $DTW_distance_metric_p_values_BH
#>          Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01       NA 0.3819444 0.3819444 0.3819444 0.3819444 0.3819444 1.0000000
#> Track_02       NA        NA 0.8250000 0.6111111 0.3819444 0.3819444 1.0000000
#> Track_03       NA        NA        NA 0.3819444 0.3819444 0.3819444 0.8250000
#> Track_04       NA        NA        NA        NA 0.3819444 0.3819444 0.3819444
#> Track_07       NA        NA        NA        NA        NA 0.6111111 0.8250000
#> Track_08       NA        NA        NA        NA        NA        NA 0.8250000
#> Track_09       NA        NA        NA        NA        NA        NA        NA
#> Track_13       NA        NA        NA        NA        NA        NA        NA
#> Track_15       NA        NA        NA        NA        NA        NA        NA
#> Track_16       NA        NA        NA        NA        NA        NA        NA
#> Track_18       NA        NA        NA        NA        NA        NA        NA
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 0.3819444 0.3819444 0.3819444 0.3819444
#> Track_02 0.6111111 0.3819444 0.3819444 0.6111111
#> Track_03 0.6111111 0.3819444 0.3819444 0.3819444
#> Track_04 0.3819444 0.3819444 0.3819444 0.3819444
#> Track_07 0.6111111 0.3819444 0.3819444 0.3819444
#> Track_08 0.6111111 0.3819444 0.3819444 0.6111111
#> Track_09 0.3819444 1.0000000 1.0000000 0.8250000
#> Track_13        NA 0.3819444 1.0000000 0.6111111
#> Track_15        NA        NA 0.3819444 0.3819444
#> Track_16        NA        NA        NA 0.3819444
#> Track_18        NA        NA        NA        NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08  Track_09
#> Track_01       NA 74.58195 78.91397 96.64398 70.79938 42.08277 35.468019
#> Track_02       NA       NA 63.99860 15.82293 23.15382 63.20830 70.717161
#> Track_03       NA       NA       NA 71.73955 34.28973 31.76544 44.275450
#> Track_04       NA       NA       NA       NA 18.85089 78.35770 87.732235
#> Track_07       NA       NA       NA       NA       NA 42.76890 50.107566
#> Track_08       NA       NA       NA       NA       NA       NA  9.465876
#> Track_09       NA       NA       NA       NA       NA       NA        NA
#> Track_13       NA       NA       NA       NA       NA       NA        NA
#> Track_15       NA       NA       NA       NA       NA       NA        NA
#> Track_16       NA       NA       NA       NA       NA       NA        NA
#> Track_18       NA       NA       NA       NA       NA       NA        NA
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 93.33423 37.413138 52.051414 44.623207
#> Track_02 65.57425 64.313729 27.180051 64.332909
#> Track_03 14.28249 36.947737 17.754758 30.247222
#> Track_04 69.93651 80.012922 30.825070 79.119289
#> Track_07 32.87243 46.118442  9.575422 42.182148
#> Track_08 46.37700  5.752637 23.576983  2.286923
#> Track_09 58.16662  5.798243 33.668082 12.004624
#> Track_13       NA 51.806516 18.002230 44.967019
#> Track_15       NA        NA 25.658456  7.668531
#> Track_16       NA        NA        NA 23.322257
#> Track_18       NA        NA        NA        NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_01 Track_02 Track_03 Track_04 Track_07  Track_08 Track_09
#> Track_01       NA 54.23687 45.21874 49.44016 71.50794 66.200748 50.16937
#> Track_02       NA       NA  7.59947 12.07591 40.92355 35.146152 75.64967
#> Track_03       NA       NA       NA 17.30699 24.03585 19.747136 53.78668
#> Track_04       NA       NA       NA       NA 54.70680 48.159534 81.22320
#> Track_07       NA       NA       NA       NA       NA  4.641382 57.88723
#> Track_08       NA       NA       NA       NA       NA        NA 56.10677
#> Track_09       NA       NA       NA       NA       NA        NA       NA
#> Track_13       NA       NA       NA       NA       NA        NA       NA
#> Track_15       NA       NA       NA       NA       NA        NA       NA
#> Track_16       NA       NA       NA       NA       NA        NA       NA
#> Track_18       NA       NA       NA       NA       NA        NA       NA
#>          Track_13  Track_15 Track_16 Track_18
#> Track_01 70.35715 38.246414 53.14282 26.44397
#> Track_02 79.04307 59.444327 37.94523 16.21459
#> Track_03 54.32899 40.387821 20.39138 12.63175
#> Track_04 89.50376 65.038178 48.86099 14.16531
#> Track_07 52.78516 47.239386 11.03906 31.14017
#> Track_08 49.82184 42.374859 11.26962 27.12250
#> Track_09 23.79028  5.888145 35.99928 42.43423
#> Track_13       NA 22.474039 27.90860 47.60924
#> Track_15       NA        NA 25.48473 30.77810
#> Track_16       NA        NA       NA 21.73607
#> Track_18       NA        NA       NA       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>          Track_01 Track_02  Track_03  Track_04 Track_07 Track_08 Track_09
#> Track_01       NA 113.7108 75.852344 122.69581 77.46999 69.19479 43.15262
#> Track_02       NA       NA  8.193052  35.02683 19.17839 44.74527 69.97435
#> Track_03       NA       NA        NA  24.76149 14.27356 28.90253 43.40454
#> Track_04       NA       NA        NA        NA 47.01337 71.48320 65.05010
#> Track_07       NA       NA        NA        NA       NA 19.83662 51.52343
#> Track_08       NA       NA        NA        NA       NA       NA 52.01348
#> Track_09       NA       NA        NA        NA       NA       NA       NA
#> Track_13       NA       NA        NA        NA       NA       NA       NA
#> Track_15       NA       NA        NA        NA       NA       NA       NA
#> Track_16       NA       NA        NA        NA       NA       NA       NA
#> Track_18       NA       NA        NA        NA       NA       NA       NA
#>           Track_13  Track_15 Track_16  Track_18
#> Track_01 89.950716 75.784943 37.31391 75.446867
#> Track_02 26.529842 19.578187 42.36166 23.277897
#> Track_03 20.553526  9.494872 22.86946 16.753630
#> Track_04 58.201258 11.676309 58.04265 50.430771
#> Track_07  5.254308 23.198154 21.30381  5.295366
#> Track_08 21.428358 36.170033 13.73617 15.205568
#> Track_09 62.329274 39.511223 31.33146 50.831737
#> Track_13        NA 31.105605 26.59231  3.924337
#> Track_15        NA        NA 25.54778 24.486652
#> Track_16        NA        NA       NA 18.784094
#> Track_18        NA        NA       NA        NA
#> 
#> 
```
