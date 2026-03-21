# Calculate velocities and relative stride lengths for tracks

`velocity_track()` calculates the relative stride lengths and velocities
for each step in a series of tracks, based on the step length, height at
the hip, and gravity acceleration.

## Usage

``` r
velocity_track(data, H, G = NULL, method = NULL)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- H:

  A numeric vector representing the height at the hip (in meters) for
  each track maker. The length of this vector should match the number of
  tracks in the data.

- G:

  Gravity acceleration (in meters per second squared). Default is `9.8`.

- method:

  A character vector specifying the method to calculate velocities for
  each track. Method `"A"` follows the approach from Alexander (1976),
  while method `"B"` is based on Ruiz & Torices (2013). If `NULL`,
  method `"A"` will be used for all tracks.

## Value

A `track velocity` R object consisting of a list of lists, where each
sublist contains the computed parameters for a corresponding track. The
parameters included are:

- `Step_velocities`: A vector of velocities for each step in the track
  (in meters per second).

- `Mean_velocity`: The mean velocity across all steps in the track (in
  meters per second).

- `Standard_deviation_velocity`: The standard deviation of velocities
  across all steps in the track (in meters per second).

- `Maximum_velocity`: The maximum velocity among all steps in the track
  (in meters per second).

- `Minimum_velocity`: The minimum velocity among all steps in the track
  (in meters per second).

- `Step_relative_stride`: A vector of relative stride lengths for each
  step in the track (dimensionless).

- `Mean_relative_stride`: The mean relative stride length across all
  steps in the track (dimensionless).

- `Standard_deviation_relative_stride`: The standard deviation of
  relative stride lengths across all steps in the track (dimensionless).

- `Maximum_relative_stride`: The maximum relative stride length among
  all steps in the track (dimensionless).

- `Minimum_relative_stride`: The minimum relative stride length among
  all steps in the track (dimensionless).

## Details

The `velocity_track()` estimates speed from stride length using
classical formulas.

As shown by Prescott et al. (2025), such estimates may be misleading—
particularly when tracks are produced on unconsolidated substrates,
where actual speeds can be substantially lower than calculated.
Moreover, equal stride lengths may correspond to different velocities,
and conversely, different stride lengths can sometimes reflect similar
velocities, depending on gait and substrate conditions. These values
should therefore be treated with caution, both in this function and in
subsequent functions that incorporate velocity estimates. By contrast,
when trackways are formed on firmer or semi-consolidated surfaces, the
formulas are more reliable and can provide useful estimates. Use is best
regarded as comparative or historical unless the substrate context
supports more confident interpretation.

\#' The `velocity_track()` function calculates velocities using two
methods:

**Method A**: Based on Alexander (1976), with the formula: \$\$v = 0.25
\cdot \sqrt{G} \cdot S^{1.67} \cdot H^{-1.17}\$\$

- **v**: Velocity of the track-maker (in meters per second).

- **G**: Acceleration due to gravity (in meters per second squared),
  typically \\9.81\\ \text{m/s}^2\\.

- **S**: Stride length, which is the distance between consecutive
  footprints (in meters).

- **H**: Height at the hip of the track-maker (in meters).

- The coefficients \\0.25\\, \\1.67\\, and \\-1.17\\ are derived from
  empirical studies. These coefficients adjust the formula to account
  for different animal sizes and gaits.

This method applies to a wide range of terrestrial vertebrates and is
used to estimate velocity across different gaits.

**Method B**: Based on Ruiz & Torices (2013), with the formula: \$\$v =
0.226 \cdot \sqrt{G} \cdot S^{1.67} \cdot H^{-1.17}\$\$

- **v**: Velocity of the track-maker (in meters per second).

- **G**: Acceleration due to gravity (in meters per second squared),
  typically \\9.81\\ \text{m/s}^2\\.

- **S**: Stride length (in meters).

- **H**: Height at the hip of the track-maker (in meters).

- The oefficient \\0.226\\ in method B is a refinement based on updated
  data for bipedal locomotion.

Based on Thulborn & Wade (1984), it is possible to identify the gaits of
track-makers on the basis of relative stride length, as follows:

- **Walk**: \\A/H \< 2.0\\; locomotor performance equivalent to walking
  in mammals.

- **Trot**: \\2.0 \leq A/H \leq 2.9\\; locomotor performance equivalent
  to trotting or racking in mammals.

- **Run**: \\A/H \> 2.9\\; locomotor performance equivalent to
  cantering, galloping, or sprinting in mammals.

## Logo

![](figures/Logo.png)

## References

Alexander, R. M. (1976). Estimates of speeds of dinosaurs. Nature,
261(5556), 129-130.

Prescott, T. L., Griffin, B. W., Demuth, O. E., Gatesy, S. M.,
Lallensack, J. N., & Falkingham, P. L. (2025). Speed from fossil
trackways: calculations not validated by extant birds on compliant
substrates. Biology Letters, 21(6), 20250191.

Ruiz, J., & Torices, A. (2013). Humans running at stadiums and beaches
and the accuracy of speed estimations from fossil trackways. Ichnos,
20(1), 31-35.

Thulborn, R. A., & Wade, M. (1984). Dinosaur trackways in the Winton
Formation (mid-Cretaceous) of Queensland. Memoirs of the Queensland
Museum, 21(2), 413-517.

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md)

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
# Example 1: Calculate velocities for the MountTom dataset using default settings.
# H_mounttom contains hip heights for each track in the MountTom dataset.
# The function will use the default method "A" for all tracks.
# Hip heights are inferred as four times the footprint length, following Alexander's approach.
H_mounttom <- c(
  1.380, 1.404, 1.320, 1.736, 1.364, 1.432, 1.508, 1.768, 1.600,
  1.848, 1.532, 1.532, 0.760, 1.532, 1.688, 1.620, 0.636, 1.784, 1.676, 1.872,
  1.648, 1.760, 1.612
)
velocity_track(MountTom, H = H_mounttom)
#> $Track_01
#> $Track_01$Step_velocities
#> [1] 1.917671 1.699122 1.724574 1.690528 1.646310 1.521394 1.326607 1.634729
#> 
#> $Track_01$Mean_velocity
#> [1] 1.645117
#> 
#> $Track_01$Standard_deviation_velocity
#> [1] 0.17006
#> 
#> $Track_01$Maximum_velocity
#> [1] 1.917671
#> 
#> $Track_01$Minimum_velocity
#> [1] 1.326607
#> 
#> $Track_01$Step_relative_stride
#> [1] 1.553054 1.444507 1.457425 1.440128 1.417452 1.352034 1.245544 1.411473
#> 
#> $Track_01$Mean_relative_stride
#> [1] 1.415202
#> 
#> $Track_01$Standard_deviation_relative_stride
#> [1] 0.08868808
#> 
#> $Track_01$Maximum_relative_stride
#> [1] 1.553054
#> 
#> $Track_01$Minimum_relative_stride
#> [1] 1.245544
#> 
#> 
#> $Track_02
#> $Track_02$Step_velocities
#> [1] 1.446657 1.678316 1.833994 1.901630 2.320863 1.891043 1.662219 1.903494
#> 
#> $Track_02$Mean_velocity
#> [1] 1.829777
#> 
#> $Track_02$Standard_deviation_velocity
#> [1] 0.2544829
#> 
#> $Track_02$Maximum_velocity
#> [1] 2.320863
#> 
#> $Track_02$Minimum_velocity
#> [1] 1.446657
#> 
#> $Track_02$Step_relative_stride
#> [1] 1.305107 1.426506 1.504327 1.537305 1.732091 1.532175 1.418298 1.538208
#> 
#> $Track_02$Mean_relative_stride
#> [1] 1.499252
#> 
#> $Track_02$Standard_deviation_relative_stride
#> [1] 0.1241122
#> 
#> $Track_02$Maximum_relative_stride
#> [1] 1.732091
#> 
#> $Track_02$Minimum_relative_stride
#> [1] 1.305107
#> 
#> 
#> $Track_03
#> $Track_03$Step_velocities
#> [1] 2.434079 2.860861 3.171859 2.291902
#> 
#> $Track_03$Mean_velocity
#> [1] 2.689676
#> 
#> $Track_03$Standard_deviation_velocity
#> [1] 0.4022265
#> 
#> $Track_03$Maximum_velocity
#> [1] 3.171859
#> 
#> $Track_03$Minimum_velocity
#> [1] 2.291902
#> 
#> $Track_03$Step_relative_stride
#> [1] 1.815427 1.999825 2.127300 1.751165
#> 
#> $Track_03$Mean_relative_stride
#> [1] 1.923429
#> 
#> $Track_03$Standard_deviation_relative_stride
#> [1] 0.1719874
#> 
#> $Track_03$Maximum_relative_stride
#> [1] 2.1273
#> 
#> $Track_03$Minimum_relative_stride
#> [1] 1.751165
#> 
#> 
#> $Track_04
#> $Track_04$Step_velocities
#> [1] 1.830533 1.623483 1.409200 1.532089
#> 
#> $Track_04$Mean_velocity
#> [1] 1.598826
#> 
#> $Track_04$Standard_deviation_velocity
#> [1] 0.1776776
#> 
#> $Track_04$Maximum_velocity
#> [1] 1.830533
#> 
#> $Track_04$Minimum_velocity
#> [1] 1.4092
#> 
#> $Track_04$Step_relative_stride
#> [1] 1.410104 1.312308 1.205659 1.267558
#> 
#> $Track_04$Mean_relative_stride
#> [1] 1.298907
#> 
#> $Track_04$Standard_deviation_relative_stride
#> [1] 0.08606683
#> 
#> $Track_04$Maximum_relative_stride
#> [1] 1.410104
#> 
#> $Track_04$Minimum_relative_stride
#> [1] 1.205659
#> 
#> 
#> $Track_05
#> $Track_05$Step_velocities
#> [1] 2.175188
#> 
#> $Track_05$Mean_velocity
#> [1] 2.175188
#> 
#> $Track_05$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_05$Maximum_velocity
#> [1] 2.175188
#> 
#> $Track_05$Minimum_velocity
#> [1] 2.175188
#> 
#> $Track_05$Step_relative_stride
#> [1] 1.680626
#> 
#> $Track_05$Mean_relative_stride
#> [1] 1.680626
#> 
#> $Track_05$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_05$Maximum_relative_stride
#> [1] 1.680626
#> 
#> $Track_05$Minimum_relative_stride
#> [1] 1.680626
#> 
#> 
#> $Track_06
#> $Track_06$Step_velocities
#> [1] 1.868553 1.995720
#> 
#> $Track_06$Mean_velocity
#> [1] 1.932136
#> 
#> $Track_06$Standard_deviation_velocity
#> [1] 0.08992089
#> 
#> $Track_06$Maximum_velocity
#> [1] 1.99572
#> 
#> $Track_06$Minimum_velocity
#> [1] 1.868553
#> 
#> $Track_06$Step_relative_stride
#> [1] 1.512270 1.573083
#> 
#> $Track_06$Mean_relative_stride
#> [1] 1.542676
#> 
#> $Track_06$Standard_deviation_relative_stride
#> [1] 0.04300131
#> 
#> $Track_06$Maximum_relative_stride
#> [1] 1.573083
#> 
#> $Track_06$Minimum_relative_stride
#> [1] 1.51227
#> 
#> 
#> $Track_07
#> $Track_07$Step_velocities
#> [1] 1.1580528 0.8799351 1.1565610 1.5947244 1.6218920 1.3654251
#> 
#> $Track_07$Mean_velocity
#> [1] 1.296098
#> 
#> $Track_07$Standard_deviation_velocity
#> [1] 0.2869998
#> 
#> $Track_07$Maximum_velocity
#> [1] 1.621892
#> 
#> $Track_07$Minimum_velocity
#> [1] 0.8799351
#> 
#> $Track_07$Step_relative_stride
#> [1] 1.1181185 0.9485583 1.1172557 1.3542404 1.3680084 1.2340304
#> 
#> $Track_07$Mean_relative_stride
#> [1] 1.190035
#> 
#> $Track_07$Standard_deviation_relative_stride
#> [1] 0.1608436
#> 
#> $Track_07$Maximum_relative_stride
#> [1] 1.368008
#> 
#> $Track_07$Minimum_relative_stride
#> [1] 0.9485583
#> 
#> 
#> $Track_08
#> $Track_08$Step_velocities
#> [1] 1.840475 1.822117 1.615424 1.429871 1.893798
#> 
#> $Track_08$Mean_velocity
#> [1] 1.720337
#> 
#> $Track_08$Standard_deviation_velocity
#> [1] 0.1938159
#> 
#> $Track_08$Maximum_velocity
#> [1] 1.893798
#> 
#> $Track_08$Minimum_velocity
#> [1] 1.429871
#> 
#> $Track_08$Step_relative_stride
#> [1] 1.406970 1.398550 1.301268 1.209585 1.431239
#> 
#> $Track_08$Mean_relative_stride
#> [1] 1.349522
#> 
#> $Track_08$Standard_deviation_relative_stride
#> [1] 0.09259133
#> 
#> $Track_08$Maximum_relative_stride
#> [1] 1.431239
#> 
#> $Track_08$Minimum_relative_stride
#> [1] 1.209585
#> 
#> 
#> $Track_09
#> $Track_09$Step_velocities
#> [1] 2.836024 2.628135 2.426088 1.824772
#> 
#> $Track_09$Mean_velocity
#> [1] 2.428755
#> 
#> $Track_09$Standard_deviation_velocity
#> [1] 0.4360513
#> 
#> $Track_09$Maximum_velocity
#> [1] 2.836024
#> 
#> $Track_09$Minimum_velocity
#> [1] 1.824772
#> 
#> $Track_09$Step_relative_stride
#> [1] 1.878065 1.794374 1.710448 1.442246
#> 
#> $Track_09$Mean_relative_stride
#> [1] 1.706283
#> 
#> $Track_09$Standard_deviation_relative_stride
#> [1] 0.188858
#> 
#> $Track_09$Maximum_relative_stride
#> [1] 1.878065
#> 
#> $Track_09$Minimum_relative_stride
#> [1] 1.442246
#> 
#> 
#> $Track_10
#> $Track_10$Step_velocities
#> [1] 1.531647
#> 
#> $Track_10$Mean_velocity
#> [1] 1.531647
#> 
#> $Track_10$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_10$Maximum_velocity
#> [1] 1.531647
#> 
#> $Track_10$Minimum_velocity
#> [1] 1.531647
#> 
#> $Track_10$Step_relative_stride
#> [1] 1.243837
#> 
#> $Track_10$Mean_relative_stride
#> [1] 1.243837
#> 
#> $Track_10$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_10$Maximum_relative_stride
#> [1] 1.243837
#> 
#> $Track_10$Minimum_relative_stride
#> [1] 1.243837
#> 
#> 
#> $Track_11
#> $Track_11$Step_velocities
#> [1] 2.076045
#> 
#> $Track_11$Mean_velocity
#> [1] 2.076045
#> 
#> $Track_11$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_11$Maximum_velocity
#> [1] 2.076045
#> 
#> $Track_11$Minimum_velocity
#> [1] 2.076045
#> 
#> $Track_11$Step_relative_stride
#> [1] 1.578469
#> 
#> $Track_11$Mean_relative_stride
#> [1] 1.578469
#> 
#> $Track_11$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_11$Maximum_relative_stride
#> [1] 1.578469
#> 
#> $Track_11$Minimum_relative_stride
#> [1] 1.578469
#> 
#> 
#> $Track_12
#> $Track_12$Step_velocities
#> [1] 2.055415
#> 
#> $Track_12$Mean_velocity
#> [1] 2.055415
#> 
#> $Track_12$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_12$Maximum_velocity
#> [1] 2.055415
#> 
#> $Track_12$Minimum_velocity
#> [1] 2.055415
#> 
#> $Track_12$Step_relative_stride
#> [1] 1.569058
#> 
#> $Track_12$Mean_relative_stride
#> [1] 1.569058
#> 
#> $Track_12$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_12$Maximum_relative_stride
#> [1] 1.569058
#> 
#> $Track_12$Minimum_relative_stride
#> [1] 1.569058
#> 
#> 
#> $Track_13
#> $Track_13$Step_velocities
#> [1] 1.0479855 0.9131596 0.7996623 0.9715874
#> 
#> $Track_13$Mean_velocity
#> [1] 0.9330987
#> 
#> $Track_13$Standard_deviation_velocity
#> [1] 0.1046951
#> 
#> $Track_13$Maximum_velocity
#> [1] 1.047985
#> 
#> $Track_13$Minimum_velocity
#> [1] 0.7996623
#> 
#> $Track_13$Step_relative_stride
#> [1] 1.293046 1.190695 1.099729 1.235746
#> 
#> $Track_13$Mean_relative_stride
#> [1] 1.204804
#> 
#> $Track_13$Standard_deviation_relative_stride
#> [1] 0.08161704
#> 
#> $Track_13$Maximum_relative_stride
#> [1] 1.293046
#> 
#> $Track_13$Minimum_relative_stride
#> [1] 1.099729
#> 
#> 
#> $Track_14
#> $Track_14$Step_velocities
#> [1] 1.642322 1.384645
#> 
#> $Track_14$Mean_velocity
#> [1] 1.513483
#> 
#> $Track_14$Standard_deviation_velocity
#> [1] 0.1822052
#> 
#> $Track_14$Maximum_velocity
#> [1] 1.642322
#> 
#> $Track_14$Minimum_velocity
#> [1] 1.384645
#> 
#> $Track_14$Step_relative_stride
#> [1] 1.371800 1.238534
#> 
#> $Track_14$Mean_relative_stride
#> [1] 1.305167
#> 
#> $Track_14$Standard_deviation_relative_stride
#> [1] 0.09423395
#> 
#> $Track_14$Maximum_relative_stride
#> [1] 1.3718
#> 
#> $Track_14$Minimum_relative_stride
#> [1] 1.238534
#> 
#> 
#> $Track_15
#> $Track_15$Step_velocities
#> [1] 1.481571 1.819262 1.714949 1.727258
#> 
#> $Track_15$Mean_velocity
#> [1] 1.68576
#> 
#> $Track_15$Standard_deviation_velocity
#> [1] 0.1438632
#> 
#> $Track_15$Maximum_velocity
#> [1] 1.819262
#> 
#> $Track_15$Minimum_velocity
#> [1] 1.481571
#> 
#> $Track_15$Step_relative_stride
#> [1] 1.252836 1.416743 1.367525 1.373394
#> 
#> $Track_15$Mean_relative_stride
#> [1] 1.352624
#> 
#> $Track_15$Standard_deviation_relative_stride
#> [1] 0.0700531
#> 
#> $Track_15$Maximum_relative_stride
#> [1] 1.416743
#> 
#> $Track_15$Minimum_relative_stride
#> [1] 1.252836
#> 
#> 
#> $Track_16
#> $Track_16$Step_velocities
#> [1] 1.204506 1.319338 1.530753 1.656302 1.486676 1.558870 1.438280
#> 
#> $Track_16$Mean_velocity
#> [1] 1.456389
#> 
#> $Track_16$Standard_deviation_velocity
#> [1] 0.1524638
#> 
#> $Track_16$Maximum_velocity
#> [1] 1.656302
#> 
#> $Track_16$Minimum_velocity
#> [1] 1.204506
#> 
#> $Track_16$Step_relative_stride
#> [1] 1.120470 1.183262 1.293401 1.355917 1.270970 1.307575 1.246030
#> 
#> $Track_16$Mean_relative_stride
#> [1] 1.253947
#> 
#> $Track_16$Standard_deviation_relative_stride
#> [1] 0.07957735
#> 
#> $Track_16$Maximum_relative_stride
#> [1] 1.355917
#> 
#> $Track_16$Minimum_relative_stride
#> [1] 1.12047
#> 
#> 
#> $Track_17
#> $Track_17$Step_velocities
#> [1] 0.5915414
#> 
#> $Track_17$Mean_velocity
#> [1] 0.5915414
#> 
#> $Track_17$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_17$Maximum_velocity
#> [1] 0.5915414
#> 
#> $Track_17$Minimum_velocity
#> [1] 0.5915414
#> 
#> $Track_17$Step_relative_stride
#> [1] 0.9683894
#> 
#> $Track_17$Mean_relative_stride
#> [1] 0.9683894
#> 
#> $Track_17$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_17$Maximum_relative_stride
#> [1] 0.9683894
#> 
#> $Track_17$Minimum_relative_stride
#> [1] 0.9683894
#> 
#> 
#> $Track_18
#> $Track_18$Step_velocities
#> [1] 1.855606 1.676314 1.873988 1.404891
#> 
#> $Track_18$Mean_velocity
#> [1] 1.7027
#> 
#> $Track_18$Standard_deviation_velocity
#> [1] 0.2176438
#> 
#> $Track_18$Maximum_velocity
#> [1] 1.873988
#> 
#> $Track_18$Minimum_velocity
#> [1] 1.404891
#> 
#> $Track_18$Step_relative_stride
#> [1] 1.410077 1.326836 1.418425 1.193662
#> 
#> $Track_18$Mean_relative_stride
#> [1] 1.33725
#> 
#> $Track_18$Standard_deviation_relative_stride
#> [1] 0.1042735
#> 
#> $Track_18$Maximum_relative_stride
#> [1] 1.418425
#> 
#> $Track_18$Minimum_relative_stride
#> [1] 1.193662
#> 
#> 
#> $Track_19
#> $Track_19$Step_velocities
#> [1] 2.650903
#> 
#> $Track_19$Mean_velocity
#> [1] 2.650903
#> 
#> $Track_19$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_19$Maximum_velocity
#> [1] 2.650903
#> 
#> $Track_19$Minimum_velocity
#> [1] 2.650903
#> 
#> $Track_19$Step_relative_stride
#> [1] 1.778779
#> 
#> $Track_19$Mean_relative_stride
#> [1] 1.778779
#> 
#> $Track_19$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_19$Maximum_relative_stride
#> [1] 1.778779
#> 
#> $Track_19$Minimum_relative_stride
#> [1] 1.778779
#> 
#> 
#> $Track_20
#> $Track_20$Step_velocities
#> [1] 1.631206 1.396409
#> 
#> $Track_20$Mean_velocity
#> [1] 1.513807
#> 
#> $Track_20$Standard_deviation_velocity
#> [1] 0.1660266
#> 
#> $Track_20$Maximum_velocity
#> [1] 1.631206
#> 
#> $Track_20$Minimum_velocity
#> [1] 1.396409
#> 
#> $Track_20$Step_relative_stride
#> [1] 1.286657 1.172319
#> 
#> $Track_20$Mean_relative_stride
#> [1] 1.229488
#> 
#> $Track_20$Standard_deviation_relative_stride
#> [1] 0.08084897
#> 
#> $Track_20$Maximum_relative_stride
#> [1] 1.286657
#> 
#> $Track_20$Minimum_relative_stride
#> [1] 1.172319
#> 
#> 
#> $Track_21
#> $Track_21$Step_velocities
#> [1] 1.719009 1.860508
#> 
#> $Track_21$Mean_velocity
#> [1] 1.789759
#> 
#> $Track_21$Standard_deviation_velocity
#> [1] 0.1000548
#> 
#> $Track_21$Maximum_velocity
#> [1] 1.860508
#> 
#> $Track_21$Minimum_velocity
#> [1] 1.719009
#> 
#> $Track_21$Step_relative_stride
#> [1] 1.379331 1.446237
#> 
#> $Track_21$Mean_relative_stride
#> [1] 1.412784
#> 
#> $Track_21$Standard_deviation_relative_stride
#> [1] 0.04730944
#> 
#> $Track_21$Maximum_relative_stride
#> [1] 1.446237
#> 
#> $Track_21$Minimum_relative_stride
#> [1] 1.379331
#> 
#> 
#> $Track_22
#> $Track_22$Step_velocities
#> [1] 0.917559 1.250737
#> 
#> $Track_22$Mean_velocity
#> [1] 1.084148
#> 
#> $Track_22$Standard_deviation_velocity
#> [1] 0.2355923
#> 
#> $Track_22$Maximum_velocity
#> [1] 1.250737
#> 
#> $Track_22$Minimum_velocity
#> [1] 0.917559
#> 
#> $Track_22$Step_relative_stride
#> [1] 0.9286649 1.1179364
#> 
#> $Track_22$Mean_relative_stride
#> [1] 1.023301
#> 
#> $Track_22$Standard_deviation_relative_stride
#> [1] 0.1338351
#> 
#> $Track_22$Maximum_relative_stride
#> [1] 1.117936
#> 
#> $Track_22$Minimum_relative_stride
#> [1] 0.9286649
#> 
#> 
#> $Track_23
#> $Track_23$Step_velocities
#> [1] 3.753784
#> 
#> $Track_23$Mean_velocity
#> [1] 3.753784
#> 
#> $Track_23$Standard_deviation_velocity
#> [1] NA
#> 
#> $Track_23$Maximum_velocity
#> [1] 3.753784
#> 
#> $Track_23$Minimum_velocity
#> [1] 3.753784
#> 
#> $Track_23$Step_relative_stride
#> [1] 2.216404
#> 
#> $Track_23$Mean_relative_stride
#> [1] 2.216404
#> 
#> $Track_23$Standard_deviation_relative_stride
#> [1] NA
#> 
#> $Track_23$Maximum_relative_stride
#> [1] 2.216404
#> 
#> $Track_23$Minimum_relative_stride
#> [1] 2.216404
#> 
#> 

# Example 2: Calculate velocities for the PaluxyRiver dataset using default settings.
# H_paluxyriver contains hip heights for each track in the PaluxyRiver dataset.
# The function will use the default method "A" for all tracks.
# Hip heights are inferred as four times the footprint length, following Alexander's approach.
H_paluxyriver <- c(3.472, 2.200)
velocity_track(PaluxyRiver, H = H_paluxyriver)
#> $Track_1
#> $Track_1$Step_velocities
#>  [1] 0.2884087 0.2667234 0.2246319 0.2038125 0.2297207 0.2673225 0.2495612
#>  [8] 0.2303131 0.2352243 0.1496251 0.1616271 0.2372131 0.2426553 0.2254857
#> [15] 0.2041018 0.2211659 0.2609607 0.2515620 0.3030204 0.2519311 0.2015630
#> [22] 0.1627231 0.2052861 0.2788176 0.2525451 0.2455988 0.2475302 0.3025814
#> 
#> $Track_1$Mean_velocity
#> [1] 0.2357754
#> 
#> $Track_1$Standard_deviation_velocity
#> [1] 0.03865919
#> 
#> $Track_1$Maximum_velocity
#> [1] 0.3030204
#> 
#> $Track_1$Minimum_velocity
#> [1] 0.1496251
#> 
#> $Track_1$Step_relative_stride
#>  [1] 0.3789154 0.3615884 0.3262496 0.3077913 0.3306554 0.3620745 0.3474711
#>  [8] 0.3311657 0.3353764 0.2557886 0.2678840 0.3370715 0.3416810 0.3269916
#> [15] 0.3080528 0.3232259 0.3568899 0.3491365 0.3902964 0.3494432 0.3057525
#> [22] 0.2689704 0.3091219 0.3713188 0.3499529 0.3441569 0.3457750 0.3899578
#> 
#> $Track_1$Mean_relative_stride
#> [1] 0.3347413
#> 
#> $Track_1$Standard_deviation_relative_stride
#> [1] 0.03373255
#> 
#> $Track_1$Maximum_relative_stride
#> [1] 0.3902964
#> 
#> $Track_1$Minimum_relative_stride
#> [1] 0.2557886
#> 
#> 
#> $Track_2
#> $Track_2$Step_velocities
#>  [1] 0.4217631 0.4424676 0.5313636 0.5147757 0.5386046 0.5371984 0.4793255
#>  [8] 0.3964918 0.3818617 0.4191057 0.4207669 0.3792643 0.5615777 0.7805553
#> [15] 0.5455458 0.3713092 0.3961137 0.4110622 0.5113603 0.5568550 0.6203647
#> [22] 0.6723679 0.7106800
#> 
#> $Track_2$Mean_velocity
#> [1] 0.5043818
#> 
#> $Track_2$Standard_deviation_velocity
#> [1] 0.11228
#> 
#> $Track_2$Maximum_velocity
#> [1] 0.7805553
#> 
#> $Track_2$Minimum_velocity
#> [1] 0.3713092
#> 
#> $Track_2$Step_relative_stride
#>  [1] 0.5453924 0.5612701 0.6263007 0.6145187 0.6313974 0.6304098 0.5888162
#>  [8] 0.5255821 0.5138819 0.5433321 0.5446207 0.5117860 0.6473884 0.7884788
#> [15] 0.6362574 0.5053306 0.5252820 0.5370638 0.6120740 0.6441228 0.6871563
#> [22] 0.7210902 0.7454200
#> 
#> $Track_2$Mean_relative_stride
#> [1] 0.6037814
#> 
#> $Track_2$Standard_deviation_relative_stride
#> [1] 0.0789691
#> 
#> $Track_2$Maximum_relative_stride
#> [1] 0.7884788
#> 
#> $Track_2$Minimum_relative_stride
#> [1] 0.5053306
#> 
#> 

# Example 3: Calculate velocities for the PaluxyRiver dataset using different methods
# for velocity calculation. Method "A" is used for sauropods, which is more
# appropriate for quadrupedal dinosaurs. Method "B" is used for theropods, which
# is more appropriate for bipedal dinosaurs. Hip heights are inferred as four times
# the footprint length, following Alexander's approach.
H_paluxyriver <- c(3.472, 2.200)
Method_paluxyriver <- c("A", "B")
velocity_track(PaluxyRiver, H = H_paluxyriver, method = Method_paluxyriver)
#> $Track_1
#> $Track_1$Step_velocities
#>  [1] 0.2884087 0.2667234 0.2246319 0.2038125 0.2297207 0.2673225 0.2495612
#>  [8] 0.2303131 0.2352243 0.1496251 0.1616271 0.2372131 0.2426553 0.2254857
#> [15] 0.2041018 0.2211659 0.2609607 0.2515620 0.3030204 0.2519311 0.2015630
#> [22] 0.1627231 0.2052861 0.2788176 0.2525451 0.2455988 0.2475302 0.3025814
#> 
#> $Track_1$Mean_velocity
#> [1] 0.2357754
#> 
#> $Track_1$Standard_deviation_velocity
#> [1] 0.03865919
#> 
#> $Track_1$Maximum_velocity
#> [1] 0.3030204
#> 
#> $Track_1$Minimum_velocity
#> [1] 0.1496251
#> 
#> $Track_1$Step_relative_stride
#>  [1] 0.3789154 0.3615884 0.3262496 0.3077913 0.3306554 0.3620745 0.3474711
#>  [8] 0.3311657 0.3353764 0.2557886 0.2678840 0.3370715 0.3416810 0.3269916
#> [15] 0.3080528 0.3232259 0.3568899 0.3491365 0.3902964 0.3494432 0.3057525
#> [22] 0.2689704 0.3091219 0.3713188 0.3499529 0.3441569 0.3457750 0.3899578
#> 
#> $Track_1$Mean_relative_stride
#> [1] 0.3347413
#> 
#> $Track_1$Standard_deviation_relative_stride
#> [1] 0.03373255
#> 
#> $Track_1$Maximum_relative_stride
#> [1] 0.3902964
#> 
#> $Track_1$Minimum_relative_stride
#> [1] 0.2557886
#> 
#> 
#> $Track_2
#> $Track_2$Step_velocities
#>  [1] 0.3812738 0.3999907 0.4803527 0.4653572 0.4868985 0.4856273 0.4333103
#>  [8] 0.3584286 0.3452030 0.3788715 0.3803733 0.3428550 0.5076662 0.7056220
#> [15] 0.4931734 0.3356635 0.3580868 0.3716003 0.4622697 0.5033969 0.5608097
#> [22] 0.6078206 0.6424547
#> 
#> $Track_2$Mean_velocity
#> [1] 0.4559611
#> 
#> $Track_2$Standard_deviation_velocity
#> [1] 0.1015011
#> 
#> $Track_2$Maximum_velocity
#> [1] 0.705622
#> 
#> $Track_2$Minimum_velocity
#> [1] 0.3356635
#> 
#> $Track_2$Step_relative_stride
#>  [1] 0.5453924 0.5612701 0.6263007 0.6145187 0.6313974 0.6304098 0.5888162
#>  [8] 0.5255821 0.5138819 0.5433321 0.5446207 0.5117860 0.6473884 0.7884788
#> [15] 0.6362574 0.5053306 0.5252820 0.5370638 0.6120740 0.6441228 0.6871563
#> [22] 0.7210902 0.7454200
#> 
#> $Track_2$Mean_relative_stride
#> [1] 0.6037814
#> 
#> $Track_2$Standard_deviation_relative_stride
#> [1] 0.0789691
#> 
#> $Track_2$Maximum_relative_stride
#> [1] 0.7884788
#> 
#> $Track_2$Minimum_relative_stride
#> [1] 0.5053306
#> 
#> 
```
