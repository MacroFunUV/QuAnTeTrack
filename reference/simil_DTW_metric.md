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
  *p*-value across all DTW distances (based on a global dominance
  criterion, evaluating in how many simulations the observed distances
  are smaller than the simulated ones across all trajectory pairs
  simultaneously). This indicates the overall significance of the
  observed DTW distances relative to simulations.

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
#> 2026-01-21 09:03:50.464242 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 9.979817
#> Track_2 9.979817       NA
#> ------------------------------------
#> 2026-01-21 09:03:50.470958 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 13.50682
#> Track_2 13.50682       NA
#> ------------------------------------
#> 2026-01-21 09:03:50.477431 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 13.54744
#> Track_2 13.54744       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 13.26548
#> Track_2 13.26548       NA
#> 
#> $DTW_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2     0.5      NA
#> 
#> $DTW_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA     0.5
#> Track_2     0.5      NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0.3333333
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>          Track_1  Track_2
#> Track_1       NA 9.979817
#> Track_2 9.979817       NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_1  Track_2
#> Track_1       NA 13.50682
#> Track_2 13.50682       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>          Track_1  Track_2
#> Track_1       NA 13.54744
#> Track_2 13.54744       NA
#> 
#> 

# Example 2: Simulating tracks using the "Constrained" model and comparing DTW distance
# in the PaluxyRiver dataset
s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2, superposition = "None")
#> 2026-01-21 09:03:50.500721 Iteration 1
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 21.97172
#> Track_2 21.97172       NA
#> ------------------------------------
#> 2026-01-21 09:03:50.507319 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 47.63361
#> Track_2 47.63361       NA
#> ------------------------------------
#> 2026-01-21 09:03:50.514346 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 81.83937
#> Track_2 81.83937       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 13.26548
#> Track_2 13.26548       NA
#> 
#> $DTW_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $DTW_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>          Track_1  Track_2
#> Track_1       NA 21.97172
#> Track_2 21.97172       NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_1  Track_2
#> Track_1       NA 47.63361
#> Track_2 47.63361       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>          Track_1  Track_2
#> Track_1       NA 81.83937
#> Track_2 81.83937       NA
#> 
#> 

# Example 3: Simulating tracks using the "Unconstrained" model and comparing DTW distance
# in the PaluxyRiver dataset
s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3, superposition = "None")
#> 2026-01-21 09:03:50.538233 Iteration 1
#>  
#> DTW metric
#>         Track_1 Track_2
#> Track_1      NA  76.501
#> Track_2  76.501      NA
#> ------------------------------------
#> 2026-01-21 09:03:50.544837 Iteration 2
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 223.8046
#> Track_2 223.8046       NA
#> ------------------------------------
#> 2026-01-21 09:03:50.551419 Iteration 3
#>  
#> DTW metric
#>          Track_1  Track_2
#> Track_1       NA 234.1775
#> Track_2 234.1775       NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>          Track_1  Track_2
#> Track_1       NA 13.26548
#> Track_2 13.26548       NA
#> 
#> $DTW_distance_metric_p_values
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $DTW_distance_metric_p_values_BH
#>         Track_1 Track_2
#> Track_1      NA    0.25
#> Track_2    0.25      NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>         Track_1 Track_2
#> Track_1      NA  76.501
#> Track_2  76.501      NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_1  Track_2
#> Track_1       NA 223.8046
#> Track_2 223.8046       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>          Track_1  Track_2
#> Track_1       NA 234.1775
#> Track_2 234.1775       NA
#> 
#> 

# Example 4: Simulating and comparing DTW distance in the MountTom dataset using the
# "Centroid" superposition method
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
s4 <- simulate_track(sbMountTom, nsim = 3)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.
simil_DTW_metric(sbMountTom, test = TRUE, sim = s4, superposition = "Centroid")
#> 2026-01-21 09:03:50.654567 Iteration 1
#>  
#> DTW metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA  5.099072 21.127506 41.176421 19.971324 26.676111 31.489501
#> Track_02  5.099072        NA 18.558782 38.024363 15.308206 22.915798 27.811019
#> Track_03 21.127506 18.558782        NA 23.091582 14.934968 11.961874 14.432559
#> Track_04 41.176421 38.024363 23.091582        NA 26.066571 10.182651  4.660675
#> Track_07 19.971324 15.308206 14.934968 26.066571        NA 11.266963 16.757384
#> Track_08 26.676111 22.915798 11.961874 10.182651 11.266963        NA  4.227394
#> Track_09 31.489501 27.811019 14.432559  4.660675 16.757384  4.227394        NA
#> Track_13 26.558060 24.195164  5.947835 20.186184 19.072238 12.721306 13.150728
#> Track_15 28.172648 23.912173 15.914475 16.904806  9.748940  5.070596  9.696100
#> Track_16 20.892144 17.639465  4.584891 20.550000 11.010142  7.828513 11.257876
#> Track_18 20.763096 16.303730 13.399978 23.749563  2.500395  8.816256 14.301279
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 26.558060 28.172648 20.892144 20.763096
#> Track_02 24.195164 23.912173 17.639465 16.303730
#> Track_03  5.947835 15.914475  4.584891 13.399978
#> Track_04 20.186184 16.904806 20.550000 23.749563
#> Track_07 19.072238  9.748940 11.010142  2.500395
#> Track_08 12.721306  5.070596  7.828513  8.816256
#> Track_09 13.150728  9.696100 11.257876 14.301279
#> Track_13        NA 17.424759  7.180503 17.197824
#> Track_15 17.424759        NA 11.442306  7.863099
#> Track_16  7.180503 11.442306        NA  9.162898
#> Track_18 17.197824  7.863099  9.162898        NA
#> ------------------------------------
#> 2026-01-21 09:03:50.68754 Iteration 2
#>  
#> DTW metric
#>          Track_01  Track_02  Track_03 Track_04  Track_07  Track_08  Track_09
#> Track_01       NA 33.130250 14.553396 40.21132 23.035651 15.969618 22.786523
#> Track_02 33.13025        NA 23.861378 11.06924  8.681039 17.029635 26.555642
#> Track_03 14.55340 23.861378        NA 27.93559 14.017225 10.322758  8.154330
#> Track_04 40.21132 11.069238 27.935586       NA 17.973904 24.226321 27.831762
#> Track_07 23.03565  8.681039 14.017225 17.97390        NA  6.178245 18.522523
#> Track_08 15.96962 17.029635 10.322758 24.22632  6.178245        NA 17.003510
#> Track_09 22.78652 26.555642  8.154330 27.83176 18.522523 17.003510        NA
#> Track_13 23.91321 10.951744 17.198328 20.83676  4.277418  7.752624 23.065134
#> Track_15 18.13111 16.366055 12.364290 24.07451  5.532369  2.196908 18.818904
#> Track_16 20.83465 18.319785  6.490797 21.59820  9.731494  9.227222  7.219053
#> Track_18 25.85540 16.823225 11.122091 15.95729 10.823774 12.679533 10.824884
#>           Track_13  Track_15  Track_16 Track_18
#> Track_01 23.913208 18.131111 20.834654 25.85540
#> Track_02 10.951744 16.366055 18.319785 16.82323
#> Track_03 17.198328 12.364290  6.490797 11.12209
#> Track_04 20.836762 24.074506 21.598200 15.95729
#> Track_07  4.277418  5.532369  9.731494 10.82377
#> Track_08  7.752624  2.196908  9.227222 12.67953
#> Track_09 23.065134 18.818904  7.219053 10.82488
#> Track_13        NA  6.306582 13.761201 15.15307
#> Track_15  6.306582        NA 10.783614 13.84831
#> Track_16 13.761201 10.783614        NA  5.03461
#> Track_18 15.153072 13.848306  5.034610       NA
#> ------------------------------------
#> 2026-01-21 09:03:50.729788 Iteration 3
#>  
#> DTW metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA 26.658775 24.580360 38.751433 19.990236 20.941327 26.880016
#> Track_02 26.658775        NA 11.484458 28.024911 23.047248  4.615657 12.218337
#> Track_03 24.580360 11.484458        NA 17.138166 12.928616  5.384313  2.163627
#> Track_04 38.751433 28.024911 17.138166        NA 21.011350 23.094107 15.638116
#> Track_07 19.990236 23.047248 12.928616 21.011350        NA 13.907852 14.788571
#> Track_08 20.941327  4.615657  5.384313 23.094107 13.907852        NA  7.983450
#> Track_09 26.880016 12.218337  2.163627 15.638116 14.788571  7.983450        NA
#> Track_13 16.079192 11.398434 12.525013 30.198404 16.733943  7.381549 14.848048
#> Track_15  9.203114 16.388248 12.880130 28.646032 12.615872  8.702712 15.423801
#> Track_16 16.267647 12.604499  7.202144 22.634739  9.378127  5.150364  9.784507
#> Track_18 25.090862 17.463661  5.157495  7.461897  9.526669  9.462149  6.339707
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 16.079192  9.203114 16.267647 25.090862
#> Track_02 11.398434 16.388248 12.604499 17.463661
#> Track_03 12.525013 12.880130  7.202144  5.157495
#> Track_04 30.198404 28.646032 22.634739  7.461897
#> Track_07 16.733943 12.615872  9.378127  9.526669
#> Track_08  7.381549  8.702712  5.150364  9.462149
#> Track_09 14.848048 15.423801  9.784507  6.339707
#> Track_13        NA  5.888813  6.867163 15.727535
#> Track_15  5.888813        NA  6.154785 14.516339
#> Track_16  6.867163  6.154785        NA  8.481408
#> Track_18 15.727535 14.516339  8.481408        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA  3.790482  5.551991  6.354408  6.293755  6.545517 27.473859
#> Track_02  3.790482        NA  6.252057  6.520295  6.335469  5.997973 29.219040
#> Track_03  5.551991  6.252057        NA  1.580837  3.291543  3.451254 15.731908
#> Track_04  6.354408  6.520295  1.580837        NA  2.682540  3.785809 14.493610
#> Track_07  6.293755  6.335469  3.291543  2.682540        NA  2.128353 18.268223
#> Track_08  6.545517  5.997973  3.451254  3.785809  2.128353        NA 17.906662
#> Track_09 27.473859 29.219040 15.731908 14.493610 18.268223 17.906662        NA
#> Track_13 19.805609 22.392336  8.618568  7.853303 12.238617 11.680236  9.666931
#> Track_15  8.429343  7.678979  2.405548  1.255909  2.733940  3.683404 14.242579
#> Track_16  5.274394  5.305788  4.664555  5.251493  3.962571  2.726565 23.970570
#> Track_18  9.501039  8.585502  2.914682  1.882578  3.297273  3.841897 14.279791
#>           Track_13   Track_15  Track_16   Track_18
#> Track_01 19.805609  8.4293428  5.274394  9.5010393
#> Track_02 22.392336  7.6789790  5.305788  8.5855020
#> Track_03  8.618568  2.4055480  4.664555  2.9146816
#> Track_04  7.853303  1.2559087  5.251493  1.8825775
#> Track_07 12.238617  2.7339402  3.962571  3.2972731
#> Track_08 11.680236  3.6834037  2.726565  3.8418973
#> Track_09  9.666931 14.2425792 23.970570 14.2797910
#> Track_13        NA  8.4643474 17.494560  8.8889488
#> Track_15  8.464347         NA  5.635359  0.7762141
#> Track_16 17.494560  5.6353588        NA  6.0479343
#> Track_18  8.888949  0.7762141  6.047934         NA
#> 
#> $DTW_distance_metric_p_values
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA     0.25     0.25     0.25     0.25     0.25     0.75
#> Track_02     0.25       NA     0.25     0.25     0.25     0.50     1.00
#> Track_03     0.25     0.25       NA     0.25     0.25     0.25     1.00
#> Track_04     0.25     0.25     0.25       NA     0.25     0.25     0.50
#> Track_07     0.25     0.25     0.25     0.25       NA     0.25     0.75
#> Track_08     0.25     0.50     0.25     0.25     0.25       NA     1.00
#> Track_09     0.75     1.00     1.00     0.50     0.75     1.00       NA
#> Track_13     0.50     0.75     0.50     0.25     0.50     0.75     0.25
#> Track_15     0.25     0.25     0.25     0.25     0.25     0.50     0.50
#> Track_16     0.25     0.25     0.50     0.25     0.25     0.25     1.00
#> Track_18     0.25     0.25     0.25     0.25     0.50     0.25     0.75
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01     0.50     0.25     0.25     0.25
#> Track_02     0.75     0.25     0.25     0.25
#> Track_03     0.50     0.25     0.50     0.25
#> Track_04     0.25     0.25     0.25     0.25
#> Track_07     0.50     0.25     0.25     0.50
#> Track_08     0.75     0.50     0.25     0.25
#> Track_09     0.25     0.50     1.00     0.75
#> Track_13       NA     0.75     1.00     0.25
#> Track_15     0.75       NA     0.25     0.25
#> Track_16     1.00     0.25       NA     0.50
#> Track_18     0.25     0.25     0.50       NA
#> 
#> $DTW_distance_metric_p_values_BH
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA 0.4044118 0.4044118 0.4044118 0.4044118 0.4044118 0.8250000
#> Track_02 0.4044118        NA 0.4044118 0.4044118 0.4044118 0.6250000 1.0000000
#> Track_03 0.4044118 0.4044118        NA 0.4044118 0.4044118 0.4044118 1.0000000
#> Track_04 0.4044118 0.4044118 0.4044118        NA 0.4044118 0.4044118 0.6250000
#> Track_07 0.4044118 0.4044118 0.4044118 0.4044118        NA 0.4044118 0.8250000
#> Track_08 0.4044118 0.6250000 0.4044118 0.4044118 0.4044118        NA 1.0000000
#> Track_09 0.8250000 1.0000000 1.0000000 0.6250000 0.8250000 1.0000000        NA
#> Track_13 0.6250000 0.8250000 0.6250000 0.4044118 0.6250000 0.8250000 0.4044118
#> Track_15 0.4044118 0.4044118 0.4044118 0.4044118 0.4044118 0.6250000 0.6250000
#> Track_16 0.4044118 0.4044118 0.6250000 0.4044118 0.4044118 0.4044118 1.0000000
#> Track_18 0.4044118 0.4044118 0.4044118 0.4044118 0.6250000 0.4044118 0.8250000
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 0.6250000 0.4044118 0.4044118 0.4044118
#> Track_02 0.8250000 0.4044118 0.4044118 0.4044118
#> Track_03 0.6250000 0.4044118 0.6250000 0.4044118
#> Track_04 0.4044118 0.4044118 0.4044118 0.4044118
#> Track_07 0.6250000 0.4044118 0.4044118 0.6250000
#> Track_08 0.8250000 0.6250000 0.4044118 0.4044118
#> Track_09 0.4044118 0.6250000 1.0000000 0.8250000
#> Track_13        NA 0.8250000 1.0000000 0.4044118
#> Track_15 0.8250000        NA 0.4044118 0.4044118
#> Track_16 1.0000000 0.4044118        NA 0.6250000
#> Track_18 0.4044118 0.4044118 0.6250000        NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA  5.099072 21.127506 41.176421 19.971324 26.676111 31.489501
#> Track_02  5.099072        NA 18.558782 38.024363 15.308206 22.915798 27.811019
#> Track_03 21.127506 18.558782        NA 23.091582 14.934968 11.961874 14.432559
#> Track_04 41.176421 38.024363 23.091582        NA 26.066571 10.182651  4.660675
#> Track_07 19.971324 15.308206 14.934968 26.066571        NA 11.266963 16.757384
#> Track_08 26.676111 22.915798 11.961874 10.182651 11.266963        NA  4.227394
#> Track_09 31.489501 27.811019 14.432559  4.660675 16.757384  4.227394        NA
#> Track_13 26.558060 24.195164  5.947835 20.186184 19.072238 12.721306 13.150728
#> Track_15 28.172648 23.912173 15.914475 16.904806  9.748940  5.070596  9.696100
#> Track_16 20.892144 17.639465  4.584891 20.550000 11.010142  7.828513 11.257876
#> Track_18 20.763096 16.303730 13.399978 23.749563  2.500395  8.816256 14.301279
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 26.558060 28.172648 20.892144 20.763096
#> Track_02 24.195164 23.912173 17.639465 16.303730
#> Track_03  5.947835 15.914475  4.584891 13.399978
#> Track_04 20.186184 16.904806 20.550000 23.749563
#> Track_07 19.072238  9.748940 11.010142  2.500395
#> Track_08 12.721306  5.070596  7.828513  8.816256
#> Track_09 13.150728  9.696100 11.257876 14.301279
#> Track_13        NA 17.424759  7.180503 17.197824
#> Track_15 17.424759        NA 11.442306  7.863099
#> Track_16  7.180503 11.442306        NA  9.162898
#> Track_18 17.197824  7.863099  9.162898        NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_01  Track_02  Track_03 Track_04  Track_07  Track_08  Track_09
#> Track_01       NA 33.130250 14.553396 40.21132 23.035651 15.969618 22.786523
#> Track_02 33.13025        NA 23.861378 11.06924  8.681039 17.029635 26.555642
#> Track_03 14.55340 23.861378        NA 27.93559 14.017225 10.322758  8.154330
#> Track_04 40.21132 11.069238 27.935586       NA 17.973904 24.226321 27.831762
#> Track_07 23.03565  8.681039 14.017225 17.97390        NA  6.178245 18.522523
#> Track_08 15.96962 17.029635 10.322758 24.22632  6.178245        NA 17.003510
#> Track_09 22.78652 26.555642  8.154330 27.83176 18.522523 17.003510        NA
#> Track_13 23.91321 10.951744 17.198328 20.83676  4.277418  7.752624 23.065134
#> Track_15 18.13111 16.366055 12.364290 24.07451  5.532369  2.196908 18.818904
#> Track_16 20.83465 18.319785  6.490797 21.59820  9.731494  9.227222  7.219053
#> Track_18 25.85540 16.823225 11.122091 15.95729 10.823774 12.679533 10.824884
#>           Track_13  Track_15  Track_16 Track_18
#> Track_01 23.913208 18.131111 20.834654 25.85540
#> Track_02 10.951744 16.366055 18.319785 16.82323
#> Track_03 17.198328 12.364290  6.490797 11.12209
#> Track_04 20.836762 24.074506 21.598200 15.95729
#> Track_07  4.277418  5.532369  9.731494 10.82377
#> Track_08  7.752624  2.196908  9.227222 12.67953
#> Track_09 23.065134 18.818904  7.219053 10.82488
#> Track_13        NA  6.306582 13.761201 15.15307
#> Track_15  6.306582        NA 10.783614 13.84831
#> Track_16 13.761201 10.783614        NA  5.03461
#> Track_18 15.153072 13.848306  5.034610       NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA 26.658775 24.580360 38.751433 19.990236 20.941327 26.880016
#> Track_02 26.658775        NA 11.484458 28.024911 23.047248  4.615657 12.218337
#> Track_03 24.580360 11.484458        NA 17.138166 12.928616  5.384313  2.163627
#> Track_04 38.751433 28.024911 17.138166        NA 21.011350 23.094107 15.638116
#> Track_07 19.990236 23.047248 12.928616 21.011350        NA 13.907852 14.788571
#> Track_08 20.941327  4.615657  5.384313 23.094107 13.907852        NA  7.983450
#> Track_09 26.880016 12.218337  2.163627 15.638116 14.788571  7.983450        NA
#> Track_13 16.079192 11.398434 12.525013 30.198404 16.733943  7.381549 14.848048
#> Track_15  9.203114 16.388248 12.880130 28.646032 12.615872  8.702712 15.423801
#> Track_16 16.267647 12.604499  7.202144 22.634739  9.378127  5.150364  9.784507
#> Track_18 25.090862 17.463661  5.157495  7.461897  9.526669  9.462149  6.339707
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 16.079192  9.203114 16.267647 25.090862
#> Track_02 11.398434 16.388248 12.604499 17.463661
#> Track_03 12.525013 12.880130  7.202144  5.157495
#> Track_04 30.198404 28.646032 22.634739  7.461897
#> Track_07 16.733943 12.615872  9.378127  9.526669
#> Track_08  7.381549  8.702712  5.150364  9.462149
#> Track_09 14.848048 15.423801  9.784507  6.339707
#> Track_13        NA  5.888813  6.867163 15.727535
#> Track_15  5.888813        NA  6.154785 14.516339
#> Track_16  6.867163  6.154785        NA  8.481408
#> Track_18 15.727535 14.516339  8.481408        NA
#> 
#> 

# Example 5: Simulating and comparing DTW distance in the MountTom dataset using the
# "Origin" superposition method
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
s5 <- simulate_track(sbMountTom, nsim = 3)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.
simil_DTW_metric(sbMountTom, test = TRUE, sim = s5, superposition = "Origin")
#> 2026-01-21 09:03:50.857251 Iteration 1
#>  
#> DTW metric
#>          Track_01 Track_02  Track_03 Track_04 Track_07  Track_08 Track_09
#> Track_01       NA 42.47589 40.311066 72.77273 18.23529 41.408291 31.09833
#> Track_02 42.47589       NA 11.311410 48.94278 35.77647 39.209366 11.42958
#> Track_03 40.31107 11.31141        NA 36.50607 22.45184 21.558099 12.84512
#> Track_04 72.77273 48.94278 36.506072       NA 49.59908 35.673796 48.58947
#> Track_07 18.23529 35.77647 22.451842 49.59908       NA 16.118092 22.63527
#> Track_08 41.40829 39.20937 21.558099 35.67380 16.11809        NA 29.19482
#> Track_09 31.09833 11.42958 12.845122 48.58947 22.63527 29.194819       NA
#> Track_13 50.06331 24.98533 10.696201 21.96203 29.97486 22.923788 23.43347
#> Track_15 39.62142 43.90742 27.076781 43.08614 14.38121  7.169873 32.57983
#> Track_16 34.03540 17.39151  7.417121 38.69365 14.68001 15.434633 10.04317
#> Track_18 43.67521 38.98900 20.980624 31.15723 18.57692  2.899171 29.74348
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 50.06331 39.621425 34.035401 43.675209
#> Track_02 24.98533 43.907421 17.391514 38.989004
#> Track_03 10.69620 27.076781  7.417121 20.980624
#> Track_04 21.96203 43.086141 38.693651 31.157234
#> Track_07 29.97486 14.381213 14.680007 18.576916
#> Track_08 22.92379  7.169873 15.434633  2.899171
#> Track_09 23.43347 32.579832 10.043166 29.743479
#> Track_13       NA 29.787928 15.503752 21.128593
#> Track_15 29.78793        NA 19.882182  9.784580
#> Track_16 15.50375 19.882182        NA 15.744219
#> Track_18 21.12859  9.784580 15.744219        NA
#> ------------------------------------
#> 2026-01-21 09:03:50.890335 Iteration 2
#>  
#> DTW metric
#>          Track_01 Track_02  Track_03 Track_04  Track_07  Track_08  Track_09
#> Track_01       NA 25.61680 34.467612 16.81163 43.303719 33.356983 45.231458
#> Track_02 25.61680       NA 36.071323 39.98856 29.715710 14.085289 25.214765
#> Track_03 34.46761 36.07132        NA 28.32518 24.680292 21.825850 34.196734
#> Track_04 16.81163 39.98856 28.325177       NA 49.142973 41.031904 54.415952
#> Track_07 43.30372 29.71571 24.680292 49.14297        NA 10.397980  8.831106
#> Track_08 33.35698 14.08529 21.825850 41.03190 10.397980        NA 11.115288
#> Track_09 45.23146 25.21476 34.196734 54.41595  8.831106 11.115288        NA
#> Track_13 41.13327 22.74470 29.652505 49.76163  8.453755  7.305674  2.988680
#> Track_15 34.65877 37.05261  1.910981 27.02536 26.503780 23.351369 35.902293
#> Track_16 35.74066 31.80623  7.935338 35.64915 16.394261 14.939977 26.258164
#> Track_18 44.84902 34.83599 20.169524 48.06015  7.875074 14.951753 19.289614
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 41.133269 34.658772 35.740664 44.849020
#> Track_02 22.744698 37.052611 31.806226 34.835989
#> Track_03 29.652505  1.910981  7.935338 20.169524
#> Track_04 49.761629 27.025357 35.649147 48.060148
#> Track_07  8.453755 26.503780 16.394261  7.875074
#> Track_08  7.305674 23.351369 14.939977 14.951753
#> Track_09  2.988680 35.902293 26.258164 19.289614
#> Track_13        NA 31.321012 21.885400 16.647471
#> Track_15 31.321012        NA  9.743024 22.017496
#> Track_16 21.885400  9.743024        NA 12.296557
#> Track_18 16.647471 22.017496 12.296557        NA
#> ------------------------------------
#> 2026-01-21 09:03:50.923916 Iteration 3
#>  
#> DTW metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07 Track_08 Track_09
#> Track_01        NA 50.063875  8.404453 57.131271 48.361205 42.79676 46.04793
#> Track_02 50.063875        NA 33.599038  3.016497 33.927255 12.83598 45.55582
#> Track_03  8.404453 33.599038        NA 43.168906 24.868767 19.61443 25.71570
#> Track_04 57.131271  3.016497 43.168906        NA 43.621175 18.83857 55.58202
#> Track_07 48.361205 33.927255 24.868767 43.621175        NA 14.02200 15.57339
#> Track_08 42.796763 12.835983 19.614427 18.838566 14.022005       NA 25.91430
#> Track_09 46.047934 45.555824 25.715696 55.582017 15.573390 25.91430       NA
#> Track_13 52.770640 29.619186 29.113131 38.485064  7.196171 12.98176 22.28757
#> Track_15 12.969407 38.596760  5.441627 48.295251 25.686237 23.26915 23.66086
#> Track_16 35.808576 34.208985 12.909317 44.319729 12.980594 15.09734 11.27135
#> Track_18  8.690007 33.750546  1.242563 43.406014 23.864216 19.14845 24.51367
#>           Track_13  Track_15 Track_16  Track_18
#> Track_01 52.770640 12.969407 35.80858  8.690007
#> Track_02 29.619186 38.596760 34.20898 33.750546
#> Track_03 29.113131  5.441627 12.90932  1.242563
#> Track_04 38.485064 48.295251 44.31973 43.406014
#> Track_07  7.196171 25.686237 12.98059 23.864216
#> Track_08 12.981760 23.269150 15.09734 19.148449
#> Track_09 22.287568 23.660861 11.27135 24.513669
#> Track_13        NA 31.008401 18.95217 28.274447
#> Track_15 31.008401        NA 12.76701  4.996216
#> Track_16 18.952166 12.767006       NA 11.759344
#> Track_18 28.274447  4.996216 11.75934        NA
#> ------------------------------------
#> ANALYSIS COMPLETED
#> ------------------------------------
#>  
#> $DTW_distance_metric
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08 Track_09
#> Track_01        NA  4.275589  7.090874 10.903521 12.540555  9.621188 49.49228
#> Track_02  4.275589        NA  8.654660 11.618577 11.570116  8.414117 50.23074
#> Track_03  7.090874  8.654660        NA  2.139376  5.526765  3.728201 26.31427
#> Track_04 10.903521 11.618577  2.139376        NA  3.464404  2.340593 24.70882
#> Track_07 12.540555 11.570116  5.526765  3.464404        NA  2.920349 30.42873
#> Track_08  9.621188  8.414117  3.728201  2.340593  2.920349        NA 30.78217
#> Track_09 49.492278 50.230744 26.314273 24.708824 30.428730 30.782172       NA
#> Track_13 36.401978 38.148822 13.819979 13.221657 19.869228 19.725147 16.33798
#> Track_15 14.670720 14.171360  4.001816  2.082973  2.156846  3.205273 23.79589
#> Track_16  8.279528  5.618237  7.192102  6.900596  5.063782  3.566561 40.62266
#> Track_18 16.821003 15.681131  4.982802  3.198363  3.246543  4.229595 24.19802
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 36.40198 14.670720  8.279528 16.821003
#> Track_02 38.14882 14.171360  5.618237 15.681131
#> Track_03 13.81998  4.001816  7.192102  4.982802
#> Track_04 13.22166  2.082973  6.900596  3.198363
#> Track_07 19.86923  2.156846  5.063782  3.246543
#> Track_08 19.72515  3.205273  3.566561  4.229595
#> Track_09 16.33798 23.795890 40.622664 24.198024
#> Track_13       NA 13.675181 29.381051 14.831959
#> Track_15 13.67518        NA  6.782730  1.392146
#> Track_16 29.38105  6.782730        NA  7.991596
#> Track_18 14.83196  1.392146  7.991596        NA
#> 
#> $DTW_distance_metric_p_values
#>          Track_01 Track_02 Track_03 Track_04 Track_07 Track_08 Track_09
#> Track_01       NA     0.25     0.25     0.25     0.25     0.25     1.00
#> Track_02     0.25       NA     0.25     0.50     0.25     0.25     1.00
#> Track_03     0.25     0.25       NA     0.25     0.25     0.25     0.75
#> Track_04     0.25     0.50     0.25       NA     0.25     0.25     0.25
#> Track_07     0.25     0.25     0.25     0.25       NA     0.25     1.00
#> Track_08     0.25     0.25     0.25     0.25     0.25       NA     1.00
#> Track_09     1.00     1.00     0.75     0.25     1.00     1.00       NA
#> Track_13     0.25     1.00     0.50     0.25     0.75     0.75     0.50
#> Track_15     0.50     0.25     0.50     0.25     0.25     0.25     0.50
#> Track_16     0.25     0.25     0.25     0.25     0.25     0.25     1.00
#> Track_18     0.50     0.25     0.50     0.25     0.25     0.50     0.50
#>          Track_13 Track_15 Track_16 Track_18
#> Track_01     0.25     0.50     0.25     0.50
#> Track_02     1.00     0.25     0.25     0.25
#> Track_03     0.50     0.50     0.25     0.50
#> Track_04     0.25     0.25     0.25     0.25
#> Track_07     0.75     0.25     0.25     0.25
#> Track_08     0.75     0.25     0.25     0.50
#> Track_09     0.50     0.50     1.00     0.50
#> Track_13       NA     0.25     1.00     0.25
#> Track_15     0.25       NA     0.25     0.25
#> Track_16     1.00     0.25       NA     0.25
#> Track_18     0.25     0.25     0.25       NA
#> 
#> $DTW_distance_metric_p_values_BH
#>           Track_01  Track_02  Track_03  Track_04  Track_07  Track_08  Track_09
#> Track_01        NA 0.3928571 0.3928571 0.3928571 0.3928571 0.3928571 1.0000000
#> Track_02 0.3928571        NA 0.3928571 0.6111111 0.3928571 0.3928571 1.0000000
#> Track_03 0.3928571 0.3928571        NA 0.3928571 0.3928571 0.3928571 0.8593750
#> Track_04 0.3928571 0.6111111 0.3928571        NA 0.3928571 0.3928571 0.3928571
#> Track_07 0.3928571 0.3928571 0.3928571 0.3928571        NA 0.3928571 1.0000000
#> Track_08 0.3928571 0.3928571 0.3928571 0.3928571 0.3928571        NA 1.0000000
#> Track_09 1.0000000 1.0000000 0.8593750 0.3928571 1.0000000 1.0000000        NA
#> Track_13 0.3928571 1.0000000 0.6111111 0.3928571 0.8593750 0.8593750 0.6111111
#> Track_15 0.6111111 0.3928571 0.6111111 0.3928571 0.3928571 0.3928571 0.6111111
#> Track_16 0.3928571 0.3928571 0.3928571 0.3928571 0.3928571 0.3928571 1.0000000
#> Track_18 0.6111111 0.3928571 0.6111111 0.3928571 0.3928571 0.6111111 0.6111111
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 0.3928571 0.6111111 0.3928571 0.6111111
#> Track_02 1.0000000 0.3928571 0.3928571 0.3928571
#> Track_03 0.6111111 0.6111111 0.3928571 0.6111111
#> Track_04 0.3928571 0.3928571 0.3928571 0.3928571
#> Track_07 0.8593750 0.3928571 0.3928571 0.3928571
#> Track_08 0.8593750 0.3928571 0.3928571 0.6111111
#> Track_09 0.6111111 0.6111111 1.0000000 0.6111111
#> Track_13        NA 0.3928571 1.0000000 0.3928571
#> Track_15 0.3928571        NA 0.3928571 0.3928571
#> Track_16 1.0000000 0.3928571        NA 0.3928571
#> Track_18 0.3928571 0.3928571 0.3928571        NA
#> 
#> $DTW_metric_p_values_combined
#> [1] 0
#> 
#> $DTW_distance_metric_simulations
#> $DTW_distance_metric_simulations[[1]]
#>          Track_01 Track_02  Track_03 Track_04 Track_07  Track_08 Track_09
#> Track_01       NA 42.47589 40.311066 72.77273 18.23529 41.408291 31.09833
#> Track_02 42.47589       NA 11.311410 48.94278 35.77647 39.209366 11.42958
#> Track_03 40.31107 11.31141        NA 36.50607 22.45184 21.558099 12.84512
#> Track_04 72.77273 48.94278 36.506072       NA 49.59908 35.673796 48.58947
#> Track_07 18.23529 35.77647 22.451842 49.59908       NA 16.118092 22.63527
#> Track_08 41.40829 39.20937 21.558099 35.67380 16.11809        NA 29.19482
#> Track_09 31.09833 11.42958 12.845122 48.58947 22.63527 29.194819       NA
#> Track_13 50.06331 24.98533 10.696201 21.96203 29.97486 22.923788 23.43347
#> Track_15 39.62142 43.90742 27.076781 43.08614 14.38121  7.169873 32.57983
#> Track_16 34.03540 17.39151  7.417121 38.69365 14.68001 15.434633 10.04317
#> Track_18 43.67521 38.98900 20.980624 31.15723 18.57692  2.899171 29.74348
#>          Track_13  Track_15  Track_16  Track_18
#> Track_01 50.06331 39.621425 34.035401 43.675209
#> Track_02 24.98533 43.907421 17.391514 38.989004
#> Track_03 10.69620 27.076781  7.417121 20.980624
#> Track_04 21.96203 43.086141 38.693651 31.157234
#> Track_07 29.97486 14.381213 14.680007 18.576916
#> Track_08 22.92379  7.169873 15.434633  2.899171
#> Track_09 23.43347 32.579832 10.043166 29.743479
#> Track_13       NA 29.787928 15.503752 21.128593
#> Track_15 29.78793        NA 19.882182  9.784580
#> Track_16 15.50375 19.882182        NA 15.744219
#> Track_18 21.12859  9.784580 15.744219        NA
#> 
#> $DTW_distance_metric_simulations[[2]]
#>          Track_01 Track_02  Track_03 Track_04  Track_07  Track_08  Track_09
#> Track_01       NA 25.61680 34.467612 16.81163 43.303719 33.356983 45.231458
#> Track_02 25.61680       NA 36.071323 39.98856 29.715710 14.085289 25.214765
#> Track_03 34.46761 36.07132        NA 28.32518 24.680292 21.825850 34.196734
#> Track_04 16.81163 39.98856 28.325177       NA 49.142973 41.031904 54.415952
#> Track_07 43.30372 29.71571 24.680292 49.14297        NA 10.397980  8.831106
#> Track_08 33.35698 14.08529 21.825850 41.03190 10.397980        NA 11.115288
#> Track_09 45.23146 25.21476 34.196734 54.41595  8.831106 11.115288        NA
#> Track_13 41.13327 22.74470 29.652505 49.76163  8.453755  7.305674  2.988680
#> Track_15 34.65877 37.05261  1.910981 27.02536 26.503780 23.351369 35.902293
#> Track_16 35.74066 31.80623  7.935338 35.64915 16.394261 14.939977 26.258164
#> Track_18 44.84902 34.83599 20.169524 48.06015  7.875074 14.951753 19.289614
#>           Track_13  Track_15  Track_16  Track_18
#> Track_01 41.133269 34.658772 35.740664 44.849020
#> Track_02 22.744698 37.052611 31.806226 34.835989
#> Track_03 29.652505  1.910981  7.935338 20.169524
#> Track_04 49.761629 27.025357 35.649147 48.060148
#> Track_07  8.453755 26.503780 16.394261  7.875074
#> Track_08  7.305674 23.351369 14.939977 14.951753
#> Track_09  2.988680 35.902293 26.258164 19.289614
#> Track_13        NA 31.321012 21.885400 16.647471
#> Track_15 31.321012        NA  9.743024 22.017496
#> Track_16 21.885400  9.743024        NA 12.296557
#> Track_18 16.647471 22.017496 12.296557        NA
#> 
#> $DTW_distance_metric_simulations[[3]]
#>           Track_01  Track_02  Track_03  Track_04  Track_07 Track_08 Track_09
#> Track_01        NA 50.063875  8.404453 57.131271 48.361205 42.79676 46.04793
#> Track_02 50.063875        NA 33.599038  3.016497 33.927255 12.83598 45.55582
#> Track_03  8.404453 33.599038        NA 43.168906 24.868767 19.61443 25.71570
#> Track_04 57.131271  3.016497 43.168906        NA 43.621175 18.83857 55.58202
#> Track_07 48.361205 33.927255 24.868767 43.621175        NA 14.02200 15.57339
#> Track_08 42.796763 12.835983 19.614427 18.838566 14.022005       NA 25.91430
#> Track_09 46.047934 45.555824 25.715696 55.582017 15.573390 25.91430       NA
#> Track_13 52.770640 29.619186 29.113131 38.485064  7.196171 12.98176 22.28757
#> Track_15 12.969407 38.596760  5.441627 48.295251 25.686237 23.26915 23.66086
#> Track_16 35.808576 34.208985 12.909317 44.319729 12.980594 15.09734 11.27135
#> Track_18  8.690007 33.750546  1.242563 43.406014 23.864216 19.14845 24.51367
#>           Track_13  Track_15 Track_16  Track_18
#> Track_01 52.770640 12.969407 35.80858  8.690007
#> Track_02 29.619186 38.596760 34.20898 33.750546
#> Track_03 29.113131  5.441627 12.90932  1.242563
#> Track_04 38.485064 48.295251 44.31973 43.406014
#> Track_07  7.196171 25.686237 12.98059 23.864216
#> Track_08 12.981760 23.269150 15.09734 19.148449
#> Track_09 22.287568 23.660861 11.27135 24.513669
#> Track_13        NA 31.008401 18.95217 28.274447
#> Track_15 31.008401        NA 12.76701  4.996216
#> Track_16 18.952166 12.767006       NA 11.759344
#> Track_18 28.274447  4.996216 11.75934        NA
#> 
#> 
```
