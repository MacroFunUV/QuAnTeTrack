# Test for differences in direction means with pairwise comparisons

`test_direction()` assesses differences in mean direction across tracks
using a selected circular statistical test. It provides two options: the
parametric Watson–Williams test (default), which assumes that directions
follow a von Mises distribution and that dispersion is homogeneous
across tracks—these assumptions are automatically checked before
testing—and the non-parametric Watson–Wheeler test, which does not
require these assumptions. For datasets with more than two tracks, the
function also performs pairwise comparisons to identify which tracks
differ significantly in mean direction.

## Usage

``` r
test_direction(
  data,
  analysis = NULL,
  permutation = NULL,
  B = 1000,
  seed = NULL
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

- analysis:

  A character string specifying the type of circular analysis:
  `"Watson-Williams"` (default; parametric, assumes von Mises and
  similar concentration) or `"Watson-Wheeler"` (non-parametric).

- permutation:

  Logical or `NULL`. For `"Watson-Wheeler"` only, whether to compute a
  permutation (Monte-Carlo) *p*-value. If `NULL` (default), the function
  automatically switches to permutation when some groups have \\n \<
  10\\ or ties are detected; otherwise honors the user choice.

- B:

  Integer. Number of permutations for the Watson–Wheeler permutation
  *p*-value (default `1000`).

- seed:

  Optional integer seed for reproducibility when using permutation
  (default `NULL`).

## Value

A list with the results of the statistical analysis and diagnostic
tests:

- `assumption_results`: A list with per-track Rayleigh test results
  (statistic and *p*-value), estimated \\\kappa\\ by track, and summary
  diagnostics (`kappa_range`, `kappa_ratio`).

- `global_test`: Result of the selected k-sample circular test
  (Watson–Williams or Watson–Wheeler). If permutation was used for
  Watson–Wheeler, the object is an `"htest"` with the permutation
  *p*-value and the method labeled as
  `"Watson-Wheeler (permutation, B=...)"`.

- `pairwise`: Data frame of pairwise comparisons (test statistic, raw
  *p*-value, and Holm-adjusted *p*-value) when more than two tracks are
  present.

## Details

The `test_direction()` function performs circular data analyses using
functions primarily from the circular package. Estimation of the
concentration parameter (\\\kappa\\) relies on the `est.kappa()`
function of the CircStats package. It includes:

- **Condition Testing:**

  - **Non-uniformity within tracks:** Rayleigh test on step directions
    within each track.

  - **Similarity of concentrations:** Per-track concentration
    \\(\kappa)\\ is estimated via `est.kappa`. Large heterogeneity
    (e.g., \\\max(\kappa)/\min(\kappa) \> 2\\) or unstable estimates
    will trigger a warning, since the Watson–Williams test assumes
    approximately similar concentration across tracks.

- **Statistical Analysis:**

  - **Watson–Williams (default):** Parametric comparison of mean
    directions across tracks (assumes von Mises and similar
    concentration).

  - **Watson–Wheeler:** Non-parametric k-sample test for differences in
    direction. When some groups have \\n \< 10\\ or when there are ties
    (identical angles), the asymptotic chi-squared *p*-value can be
    unreliable. To address this, the function can compute a permutation
    (Monte-Carlo) *p*-value for the same Watson–Wheeler statistic (W).
    If `permutation = NULL` (default), this permutation *p*-value is
    used automatically under small samples or ties; otherwise, the
    user's choice is respected. A tiny internal jitter (in radians) is
    applied to angles only when ties are detected to stabilize the test;
    this is not exposed as an argument.

- **Permutation (Monte-Carlo) p-values:** For the Watson–Wheeler test,
  when group sizes are small (\\n \< 10\\) or when ties are present, the
  default chi-squared *p*-value can be unreliable. In these cases, the
  function can compute a permutation (Monte-Carlo) *p*-value. This is
  done by repeatedly shuffling group labels, recalculating the
  Watson–Wheeler statistic each time, and comparing the observed value
  to the resulting null distribution. The *p*-value is then the
  proportion of permuted statistics at least as extreme as the observed
  one. This approach provides a more accurate significance level for
  small or tied samples.

  - **Pairwise comparisons:** For datasets with more than two tracks,
    pairwise tests are performed using the same family as the selected
    global test. *P*-values are adjusted using Holm’s method. For
    Watson–Wheeler, permutation *p*-values are used pairwise whenever
    applicable by the same rule (small n or ties, or when
    `permutation = TRUE`).

- **Direction Measurement:**

  - The direction is measured in degrees. The reference direction, or 0
    degrees, is considered to be along the positive x-axis. Angles are
    measured counterclockwise from the positive x-axis, with 0 degrees
    pointing directly along this axis.

## Logo

![](figures/Logo.png)

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`plot_direction`](https://macrofunuv.github.io/QuAnTeTrack/reference/plot_direction.md)

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
# Example 1: Parametric Circular Test (Watson–Williams) for Differences in Direction Means with Pairwise Comparisons in MountTom dataset
test_direction(MountTom, analysis = "Watson-Williams")
#> Warning: The following tracks were removed from the analysis due to having 3 or fewer footprints: Track 05, Track 06, Track 10, Track 11, Track 12, Track 14, Track 17, Track 19, Track 20, Track 21, Track 22, Track 23.
#> Warning: Estimated concentration (kappa) appears heterogeneous across tracks (ratio > 2). Watson–Williams assumes similar concentration; consider 'Watson-Wheeler'.
#> $assumption_results
#> $assumption_results$rayleigh
#>          statistic     p_value
#> Track 01 0.9976355 0.000000000
#> Track 02 0.9966997 0.000000000
#> Track 03 0.9997300 0.001043937
#> Track 04 0.9966016 0.001161187
#> Track 07 0.9885045 0.000000000
#> Track 08 0.9987424 0.000000000
#> Track 09 0.9987039 0.001081750
#> Track 13 0.9861080 0.001598982
#> Track 15 0.9978616 0.001113260
#> Track 16 0.9928830 0.000000000
#> Track 18 0.9994406 0.001054538
#> 
#> $assumption_results$kappa
#>  Track 01  Track 02  Track 03  Track 04  Track 07  Track 08  Track 09  Track 13 
#> 0.3910399 0.5495348 0.7396528 2.0773420 0.4898489 0.3654361 1.5453896 4.2678523 
#>  Track 15  Track 16  Track 18 
#> 0.9589439 1.5269202 0.3557234 
#> 
#> $assumption_results$kappa_range
#> [1] 3.912129
#> 
#> $assumption_results$kappa_ratio
#> [1] 11.99767
#> 
#> 
#> $global_test
#> 
#>  Watson-Williams test for homogeneity of means
#> 
#> data:  all_circ by groups
#> F = 8.3734, df1 = 10, df2 = 58, p-value = 3.056e-08
#> sample estimates:
#> Circular Data: 
#> Type = angles 
#> Units = radians 
#> Template = none 
#> Modulo = 2pi 
#> Zero = 0 
#> Rotation = counter 
#> mean of Track 01 mean of Track 02 mean of Track 03 mean of Track 04 
#>        3.1366963        1.5386635        0.7346896        0.8368350 
#> mean of Track 07 mean of Track 08 mean of Track 09 mean of Track 13 
#>        1.4885801        3.9810186        0.6642104        5.5461275 
#> mean of Track 15 mean of Track 16 mean of Track 18 
#>        0.7123928        3.3851511        3.4537059 
#> 
#> 
#> $pairwise
#>      track1   track2    statistic      p_value          method      p_adj
#> 1  Track 01 Track 02 2.980766e+00 0.1035129494 Watson-Williams 1.00000000
#> 2  Track 01 Track 03 7.776502e+00 0.0163858701 Watson-Williams 0.73736416
#> 3  Track 01 Track 04 6.833837e+00 0.0226273246 Watson-Williams 0.97297496
#> 4  Track 01 Track 07 2.763332e+00 0.1186649872 Watson-Williams 1.00000000
#> 5  Track 01 Track 08 5.148557e-01 0.4857388086 Watson-Williams 1.00000000
#> 6  Track 01 Track 09 8.280785e+00 0.0138914601 Watson-Williams 0.65177998
#> 7  Track 01 Track 13 7.550041e+00 0.0176763298 Watson-Williams 0.77775851
#> 8  Track 01 Track 15 8.286101e+00 0.0138676592 Watson-Williams 0.65177998
#> 9  Track 01 Track 16 8.173556e-02 0.7788694550 Watson-Williams 1.00000000
#> 10 Track 01 Track 18 6.138353e-02 0.8085110205 Watson-Williams 1.00000000
#> 11 Track 02 Track 03 6.466307e-01 0.4369608130 Watson-Williams 1.00000000
#> 12 Track 02 Track 04 7.400587e-01 0.4065015185 Watson-Williams 1.00000000
#> 13 Track 02 Track 07 2.473063e-03 0.9610401580 Watson-Williams 1.00000000
#> 14 Track 02 Track 08 5.219097e+00 0.0397834262 Watson-Williams 1.00000000
#> 15 Track 02 Track 09 1.054881e+00 0.3246442132 Watson-Williams 1.00000000
#> 16 Track 02 Track 13 9.891294e+00 0.0084523834 Watson-Williams 0.43107155
#> 17 Track 02 Track 15 7.759792e-01 0.3956743468 Watson-Williams 1.00000000
#> 18 Track 02 Track 16 6.154840e+00 0.0254466930 Watson-Williams 1.00000000
#> 19 Track 02 Track 18 2.440906e+00 0.1441822582 Watson-Williams 1.00000000
#> 20 Track 03 Track 04 1.345343e-02 0.9105200633 Watson-Williams 1.00000000
#> 21 Track 03 Track 07 4.678385e-01 0.5095222992 Watson-Williams 1.00000000
#> 22 Track 03 Track 08 9.572190e+00 0.0128484173 Watson-Williams 0.61672403
#> 23 Track 03 Track 09 5.656749e-03 0.9418933046 Watson-Williams 1.00000000
#> 24 Track 03 Track 13 3.423479e+00 0.1014354369 Watson-Williams 1.00000000
#> 25 Track 03 Track 15 4.526982e-04 0.9835460333 Watson-Williams 1.00000000
#> 26 Track 03 Track 16 9.150469e+00 0.0115525818 Watson-Williams 0.56607651
#> 27 Track 03 Track 18 4.979097e+00 0.0561751241 Watson-Williams 1.00000000
#> 28 Track 04 Track 07 4.951183e-01 0.4977092314 Watson-Williams 1.00000000
#> 29 Track 04 Track 08 5.522642e+00 0.0433029414 Watson-Williams 1.00000000
#> 30 Track 04 Track 09 6.963350e-02 0.7985413218 Watson-Williams 1.00000000
#> 31 Track 04 Track 13 1.130779e+01 0.0098919936 Watson-Williams 0.49459968
#> 32 Track 04 Track 15 2.468287e-02 0.8790525990 Watson-Williams 1.00000000
#> 33 Track 04 Track 16 2.520577e+01 0.0003897153 Watson-Williams 0.02143434
#> 34 Track 04 Track 18 3.755804e+00 0.0886209669 Watson-Williams 1.00000000
#> 35 Track 07 Track 08 6.120988e+00 0.0308977509 Watson-Williams 1.00000000
#> 36 Track 07 Track 09 7.301567e-01 0.4128349469 Watson-Williams 1.00000000
#> 37 Track 07 Track 13 6.494097e+00 0.0289456538 Watson-Williams 1.00000000
#> 38 Track 07 Track 15 5.495813e-01 0.4755380098 Watson-Williams 1.00000000
#> 39 Track 07 Track 16 4.647577e+00 0.0504170910 Watson-Williams 1.00000000
#> 40 Track 07 Track 18 2.535878e+00 0.1423692496 Watson-Williams 1.00000000
#> 41 Track 08 Track 09 5.630774e+00 0.0417085367 Watson-Williams 1.00000000
#> 42 Track 08 Track 13 2.163339e+00 0.1754143713 Watson-Williams 1.00000000
#> 43 Track 08 Track 15 7.103917e+00 0.0258198163 Watson-Williams 1.00000000
#> 44 Track 08 Track 16 3.234319e-01 0.5800404171 Watson-Williams 1.00000000
#> 45 Track 08 Track 18 1.379765e-01 0.7188949441 Watson-Williams 1.00000000
#> 46 Track 09 Track 13 6.592936e+00 0.0332477844 Watson-Williams 1.00000000
#> 47 Track 09 Track 15 3.206817e-03 0.9562297400 Watson-Williams 1.00000000
#> 48 Track 09 Track 16 2.348867e+01 0.0005136178 Watson-Williams 0.02773536
#> 49 Track 09 Track 18 4.003953e+00 0.0803930903 Watson-Williams 1.00000000
#> 50 Track 13 Track 15 4.338768e+00 0.0707874930 Watson-Williams 1.00000000
#> 51 Track 13 Track 16 2.147301e+01 0.0007237989 Watson-Williams 0.03836134
#> 52 Track 13 Track 18 2.979469e+00 0.1226030066 Watson-Williams 1.00000000
#> 53 Track 15 Track 16 1.277916e+01 0.0043582591 Watson-Williams 0.22662947
#> 54 Track 15 Track 18 4.384561e+00 0.0695909887 Watson-Williams 1.00000000
#> 55 Track 16 Track 18 3.637189e-03 0.9529912151 Watson-Williams 1.00000000
#> 

# Example 2: Non-Parametric Circular Test (Watson–Wheeler) with automatic permutation under small n/ties
test_direction(MountTom, analysis = "Watson-Wheeler",
               permutation = TRUE, B = 10, seed = 42)
#> Warning: The following tracks were removed from the analysis due to having 3 or fewer footprints: Track 05, Track 06, Track 10, Track 11, Track 12, Track 14, Track 17, Track 19, Track 20, Track 21, Track 22, Track 23.
#> Warning: Estimated concentration (kappa) appears heterogeneous across tracks (ratio > 2). Watson–Williams assumes similar concentration; consider 'Watson-Wheeler'.
#> Warning: Using permutation p-value for Watson–Wheeler because some groups have n < 10.
#> $assumption_results
#> $assumption_results$rayleigh
#>          statistic     p_value
#> Track 01 0.9976355 0.000000000
#> Track 02 0.9966997 0.000000000
#> Track 03 0.9997300 0.001043937
#> Track 04 0.9966016 0.001161187
#> Track 07 0.9885045 0.000000000
#> Track 08 0.9987424 0.000000000
#> Track 09 0.9987039 0.001081750
#> Track 13 0.9861080 0.001598982
#> Track 15 0.9978616 0.001113260
#> Track 16 0.9928830 0.000000000
#> Track 18 0.9994406 0.001054538
#> 
#> $assumption_results$kappa
#>  Track 01  Track 02  Track 03  Track 04  Track 07  Track 08  Track 09  Track 13 
#> 0.3910399 0.5495348 0.7396528 2.0773420 0.4898489 0.3654361 1.5453896 4.2678523 
#>  Track 15  Track 16  Track 18 
#> 0.9589439 1.5269202 0.3557234 
#> 
#> $assumption_results$kappa_range
#> [1] 3.912129
#> 
#> $assumption_results$kappa_ratio
#> [1] 11.99767
#> 
#> 
#> $global_test
#> 
#>  Watson-Wheeler (permutation, B=10)
#> 
#> data:  
#> W = 24.978, p-value = 0.3636
#> 
#> 
#> $pairwise
#>      track1   track2 statistic    p_value                             method
#> 1  Track 01 Track 02 1.9529867 0.18181818 Watson-Wheeler (permutation, B=10)
#> 2  Track 01 Track 03 1.2543066 0.81818182 Watson-Wheeler (permutation, B=10)
#> 3  Track 01 Track 04 5.1647562 0.18181818 Watson-Wheeler (permutation, B=10)
#> 4  Track 01 Track 07 0.8576608 0.81818182 Watson-Wheeler (permutation, B=10)
#> 5  Track 01 Track 08 0.3961129 0.81818182 Watson-Wheeler (permutation, B=10)
#> 6  Track 01 Track 09 3.9935156 0.36363636 Watson-Wheeler (permutation, B=10)
#> 7  Track 01 Track 13 1.7968405 0.54545455 Watson-Wheeler (permutation, B=10)
#> 8  Track 01 Track 15 2.9171522 0.45454545 Watson-Wheeler (permutation, B=10)
#> 9  Track 01 Track 16 1.6129754 0.63636364 Watson-Wheeler (permutation, B=10)
#> 10 Track 01 Track 18 0.5777778 0.72727273 Watson-Wheeler (permutation, B=10)
#> 11 Track 02 Track 03 2.7744529 0.45454545 Watson-Wheeler (permutation, B=10)
#> 12 Track 02 Track 04 3.7520653 0.45454545 Watson-Wheeler (permutation, B=10)
#> 13 Track 02 Track 07 0.8576608 0.72727273 Watson-Wheeler (permutation, B=10)
#> 14 Track 02 Track 08 0.8062315 0.72727273 Watson-Wheeler (permutation, B=10)
#> 15 Track 02 Track 09 1.9904686 0.63636364 Watson-Wheeler (permutation, B=10)
#> 16 Track 02 Track 13 3.7520653 0.36363636 Watson-Wheeler (permutation, B=10)
#> 17 Track 02 Track 15 1.1555556 0.81818182 Watson-Wheeler (permutation, B=10)
#> 18 Track 02 Track 16 3.6113309 0.18181818 Watson-Wheeler (permutation, B=10)
#> 19 Track 02 Track 18 1.0128562 0.90909091 Watson-Wheeler (permutation, B=10)
#> 20 Track 03 Track 04 3.3249845 0.45454545 Watson-Wheeler (permutation, B=10)
#> 21 Track 03 Track 07 0.6285714 0.63636364 Watson-Wheeler (permutation, B=10)
#> 22 Track 03 Track 08 1.4782138 0.72727273 Watson-Wheeler (permutation, B=10)
#> 23 Track 03 Track 09 0.1050466 1.00000000 Watson-Wheeler (permutation, B=10)
#> 24 Track 03 Track 13 3.3249845 0.36363636 Watson-Wheeler (permutation, B=10)
#> 25 Track 03 Track 15 0.7200000 0.90909091 Watson-Wheeler (permutation, B=10)
#> 26 Track 03 Track 16 6.0910205 0.09090909 Watson-Wheeler (permutation, B=10)
#> 27 Track 03 Track 18 0.7200000 0.63636364 Watson-Wheeler (permutation, B=10)
#> 28 Track 04 Track 07 0.1684252 0.81818182 Watson-Wheeler (permutation, B=10)
#> 29 Track 04 Track 08 2.6833945 0.54545455 Watson-Wheeler (permutation, B=10)
#> 30 Track 04 Track 09 3.3249845 0.09090909 Watson-Wheeler (permutation, B=10)
#> 31 Track 04 Track 13 4.9349534 0.18181818 Watson-Wheeler (permutation, B=10)
#> 32 Track 04 Track 15 0.2750155 0.90909091 Watson-Wheeler (permutation, B=10)
#> 33 Track 04 Track 16 6.0910205 0.09090909 Watson-Wheeler (permutation, B=10)
#> 34 Track 04 Track 18 0.7200000 0.45454545 Watson-Wheeler (permutation, B=10)
#> 35 Track 07 Track 08 0.9698425 0.54545455 Watson-Wheeler (permutation, B=10)
#> 36 Track 07 Track 09 0.6285714 0.72727273 Watson-Wheeler (permutation, B=10)
#> 37 Track 07 Track 13 8.7548706 0.09090909 Watson-Wheeler (permutation, B=10)
#> 38 Track 07 Track 15 0.6285714 0.90909091 Watson-Wheeler (permutation, B=10)
#> 39 Track 07 Track 16 5.1788577 0.36363636 Watson-Wheeler (permutation, B=10)
#> 40 Track 07 Track 18 1.4255681 0.81818182 Watson-Wheeler (permutation, B=10)
#> 41 Track 08 Track 09 1.4782138 0.63636364 Watson-Wheeler (permutation, B=10)
#> 42 Track 08 Track 13 8.2290250 0.09090909 Watson-Wheeler (permutation, B=10)
#> 43 Track 08 Track 15 2.0000000 0.36363636 Watson-Wheeler (permutation, B=10)
#> 44 Track 08 Track 16 5.4696605 0.09090909 Watson-Wheeler (permutation, B=10)
#> 45 Track 08 Track 18 0.1667894 0.90909091 Watson-Wheeler (permutation, B=10)
#> 46 Track 09 Track 13 3.3249845 0.27272727 Watson-Wheeler (permutation, B=10)
#> 47 Track 09 Track 15 0.7200000 0.72727273 Watson-Wheeler (permutation, B=10)
#> 48 Track 09 Track 16 6.0910205 0.09090909 Watson-Wheeler (permutation, B=10)
#> 49 Track 09 Track 18 0.7200000 0.90909091 Watson-Wheeler (permutation, B=10)
#> 50 Track 13 Track 15 4.9349534 0.18181818 Watson-Wheeler (permutation, B=10)
#> 51 Track 13 Track 16 9.1589842 0.09090909 Watson-Wheeler (permutation, B=10)
#> 52 Track 13 Track 18 7.5399379 0.09090909 Watson-Wheeler (permutation, B=10)
#> 53 Track 15 Track 16 6.0910205 0.09090909 Watson-Wheeler (permutation, B=10)
#> 54 Track 15 Track 18 0.7200000 0.72727273 Watson-Wheeler (permutation, B=10)
#> 55 Track 16 Track 18 3.1237388 0.27272727 Watson-Wheeler (permutation, B=10)
#>    p_adj
#> 1      1
#> 2      1
#> 3      1
#> 4      1
#> 5      1
#> 6      1
#> 7      1
#> 8      1
#> 9      1
#> 10     1
#> 11     1
#> 12     1
#> 13     1
#> 14     1
#> 15     1
#> 16     1
#> 17     1
#> 18     1
#> 19     1
#> 20     1
#> 21     1
#> 22     1
#> 23     1
#> 24     1
#> 25     1
#> 26     1
#> 27     1
#> 28     1
#> 29     1
#> 30     1
#> 31     1
#> 32     1
#> 33     1
#> 34     1
#> 35     1
#> 36     1
#> 37     1
#> 38     1
#> 39     1
#> 40     1
#> 41     1
#> 42     1
#> 43     1
#> 44     1
#> 45     1
#> 46     1
#> 47     1
#> 48     1
#> 49     1
#> 50     1
#> 51     1
#> 52     1
#> 53     1
#> 54     1
#> 55     1
#> 

# Example 3: Force permutation with more replicates and a seed
test_direction(PaluxyRiver, analysis = "Watson-Wheeler",
               permutation = TRUE, B = 100, seed = 42)
#> $assumption_results
#> $assumption_results$rayleigh
#>         statistic      p_value
#> Track 1 0.9886376 5.798188e-12
#> Track 2 0.9892781 4.254234e-10
#> 
#> $assumption_results$kappa
#>   Track 1   Track 2 
#> 0.3090983 0.1880619 
#> 
#> $assumption_results$kappa_range
#> [1] 0.1210364
#> 
#> $assumption_results$kappa_ratio
#> [1] 1.643599
#> 
#> 
#> $global_test
#> 
#>  Watson-Wheeler (permutation, B=100)
#> 
#> data:  
#> W = 0.0010818, p-value = 1
#> 
#> 
#> $pairwise
#>    track1  track2   statistic p_value                              method p_adj
#> 1 Track 1 Track 2 0.001081754       1 Watson-Wheeler (permutation, B=100)     1
#> 

# Example 4: Keep asymptotic p-values (even if small n or ties), but warn
test_direction(MountTom, analysis = "Watson-Wheeler", permutation = FALSE)
#> Warning: The following tracks were removed from the analysis due to having 3 or fewer footprints: Track 05, Track 06, Track 10, Track 11, Track 12, Track 14, Track 17, Track 19, Track 20, Track 21, Track 22, Track 23.
#> Warning: Estimated concentration (kappa) appears heterogeneous across tracks (ratio > 2). Watson–Williams assumes similar concentration; consider 'Watson-Wheeler'.
#> $assumption_results
#> $assumption_results$rayleigh
#>          statistic     p_value
#> Track 01 0.9976355 0.000000000
#> Track 02 0.9966997 0.000000000
#> Track 03 0.9997300 0.001043937
#> Track 04 0.9966016 0.001161187
#> Track 07 0.9885045 0.000000000
#> Track 08 0.9987424 0.000000000
#> Track 09 0.9987039 0.001081750
#> Track 13 0.9861080 0.001598982
#> Track 15 0.9978616 0.001113260
#> Track 16 0.9928830 0.000000000
#> Track 18 0.9994406 0.001054538
#> 
#> $assumption_results$kappa
#>  Track 01  Track 02  Track 03  Track 04  Track 07  Track 08  Track 09  Track 13 
#> 0.3910399 0.5495348 0.7396528 2.0773420 0.4898489 0.3654361 1.5453896 4.2678523 
#>  Track 15  Track 16  Track 18 
#> 0.9589439 1.5269202 0.3557234 
#> 
#> $assumption_results$kappa_range
#> [1] 3.912129
#> 
#> $assumption_results$kappa_ratio
#> [1] 11.99767
#> 
#> 
#> $global_test
#> 
#>  Watson-Wheeler test for homogeneity of angles
#> 
#> data:  Track 01 and Track 02 and Track 03 and Track 04 and Track 07 and Track 08 and Track 09 and Track 13 and Track 15 and Track 16 and Track 18
#> W = 24.978, df = 20, p-value = 0.2023
#> 
#> 
#> $pairwise
#>      track1   track2 statistic    p_value         method     p_adj
#> 1  Track 01 Track 02 1.9529867 0.37662949 Watson-Wheeler 1.0000000
#> 2  Track 01 Track 03 1.2543066 0.53411010 Watson-Wheeler 1.0000000
#> 3  Track 01 Track 04 5.1647562 0.07559402 Watson-Wheeler 1.0000000
#> 4  Track 01 Track 07 0.8576608 0.65127038 Watson-Wheeler 1.0000000
#> 5  Track 01 Track 08 0.3961129 0.82032354 Watson-Wheeler 1.0000000
#> 6  Track 01 Track 09 3.9935156 0.13577478 Watson-Wheeler 1.0000000
#> 7  Track 01 Track 13 1.7968405 0.40721244 Watson-Wheeler 1.0000000
#> 8  Track 01 Track 15 2.9171522 0.23256719 Watson-Wheeler 1.0000000
#> 9  Track 01 Track 16 1.6129754 0.44642329 Watson-Wheeler 1.0000000
#> 10 Track 01 Track 18 0.5777778 0.74909543 Watson-Wheeler 1.0000000
#> 11 Track 02 Track 03 2.7744529 0.24976708 Watson-Wheeler 1.0000000
#> 12 Track 02 Track 04 3.7520653 0.15319669 Watson-Wheeler 1.0000000
#> 13 Track 02 Track 07 0.8576608 0.65127038 Watson-Wheeler 1.0000000
#> 14 Track 02 Track 08 0.8062315 0.66823474 Watson-Wheeler 1.0000000
#> 15 Track 02 Track 09 1.9904686 0.36963682 Watson-Wheeler 1.0000000
#> 16 Track 02 Track 13 3.7520653 0.15319669 Watson-Wheeler 1.0000000
#> 17 Track 02 Track 15 1.1555556 0.56114397 Watson-Wheeler 1.0000000
#> 18 Track 02 Track 16 3.6113309 0.16436504 Watson-Wheeler 1.0000000
#> 19 Track 02 Track 18 1.0128562 0.60264432 Watson-Wheeler 1.0000000
#> 20 Track 03 Track 04 3.3249845 0.18966570 Watson-Wheeler 1.0000000
#> 21 Track 03 Track 07 0.6285714 0.73031034 Watson-Wheeler 1.0000000
#> 22 Track 03 Track 08 1.4782138 0.47754021 Watson-Wheeler 1.0000000
#> 23 Track 03 Track 09 0.1050466 0.94883222 Watson-Wheeler 1.0000000
#> 24 Track 03 Track 13 3.3249845 0.18966570 Watson-Wheeler 1.0000000
#> 25 Track 03 Track 15 0.7200000 0.69767633 Watson-Wheeler 1.0000000
#> 26 Track 03 Track 16 6.0910205 0.04757203 Watson-Wheeler 1.0000000
#> 27 Track 03 Track 18 0.7200000 0.69767633 Watson-Wheeler 1.0000000
#> 28 Track 04 Track 07 0.1684252 0.91923580 Watson-Wheeler 1.0000000
#> 29 Track 04 Track 08 2.6833945 0.26140162 Watson-Wheeler 1.0000000
#> 30 Track 04 Track 09 3.3249845 0.18966570 Watson-Wheeler 1.0000000
#> 31 Track 04 Track 13 4.9349534 0.08479856 Watson-Wheeler 1.0000000
#> 32 Track 04 Track 15 0.2750155 0.87152758 Watson-Wheeler 1.0000000
#> 33 Track 04 Track 16 6.0910205 0.04757203 Watson-Wheeler 1.0000000
#> 34 Track 04 Track 18 0.7200000 0.69767633 Watson-Wheeler 1.0000000
#> 35 Track 07 Track 08 0.9698425 0.61574569 Watson-Wheeler 1.0000000
#> 36 Track 07 Track 09 0.6285714 0.73031034 Watson-Wheeler 1.0000000
#> 37 Track 07 Track 13 8.7548706 0.01255752 Watson-Wheeler 0.6781063
#> 38 Track 07 Track 15 0.6285714 0.73031034 Watson-Wheeler 1.0000000
#> 39 Track 07 Track 16 5.1788577 0.07506290 Watson-Wheeler 1.0000000
#> 40 Track 07 Track 18 1.4255681 0.49027735 Watson-Wheeler 1.0000000
#> 41 Track 08 Track 09 1.4782138 0.47754021 Watson-Wheeler 1.0000000
#> 42 Track 08 Track 13 8.2290250 0.01633390 Watson-Wheeler 0.8656968
#> 43 Track 08 Track 15 2.0000000 0.36787944 Watson-Wheeler 1.0000000
#> 44 Track 08 Track 16 5.4696605 0.06490503 Watson-Wheeler 1.0000000
#> 45 Track 08 Track 18 0.1667894 0.91998798 Watson-Wheeler 1.0000000
#> 46 Track 09 Track 13 3.3249845 0.18966570 Watson-Wheeler 1.0000000
#> 47 Track 09 Track 15 0.7200000 0.69767633 Watson-Wheeler 1.0000000
#> 48 Track 09 Track 16 6.0910205 0.04757203 Watson-Wheeler 1.0000000
#> 49 Track 09 Track 18 0.7200000 0.69767633 Watson-Wheeler 1.0000000
#> 50 Track 13 Track 15 4.9349534 0.08479856 Watson-Wheeler 1.0000000
#> 51 Track 13 Track 16 9.1589842 0.01026011 Watson-Wheeler 0.5643058
#> 52 Track 13 Track 18 7.5399379 0.02305278 Watson-Wheeler 1.0000000
#> 53 Track 15 Track 16 6.0910205 0.04757203 Watson-Wheeler 1.0000000
#> 54 Track 15 Track 18 0.7200000 0.69767633 Watson-Wheeler 1.0000000
#> 55 Track 16 Track 18 3.1237388 0.20974361 Watson-Wheeler 1.0000000
#> 
```
