# Print track parameters

`track_param()` is a function to compute and print various parameters of
trackways.

## Usage

``` r
track_param(data, gauge_size = NA)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- gauge_size:

  Numeric. Pes/manus length (or width) used to compute Gauge as
  `Trackway_width / gauge_size`. If `NULL` or `NA`, Gauge is returned as
  `NA`.

## Value

A list of lists, where each sublist contains the computed parameters for
a corresponding track. The parameters included are:

- `Turning_angles`:

  A vector of turning angles for the track (in degrees).

- `Mean_turning_angle`:

  The mean of the turning angles (in degrees).

- `Standard_deviation_turning_angle`:

  The standard deviation of the turning angles (in degrees).

- `Path_length`:

  Total path length (in meters), computed as the sum of distances
  between consecutive trajectory points.

- `Beeline_length`:

  Straight-line distance (in meters) between the first and last points
  of the trajectory.

- `Step_lengths`:

  A vector of step lengths for the track (in meters).

- `Mean_step_length`:

  The mean of the step lengths (in meters).

- `Standard_deviation_step_length`:

  The standard deviation of the step lengths (in meters).

- `Stride_length`:

  A vector of stride lengths (in meters), computed as distances between
  consecutive footprints of the same side (L–L and R–R), in footprint
  order.

- `Mean_stride_length`:

  Mean stride length (in meters). Computed as `Mean_step_length * 2`.

- `Pace_length`:

  A vector of pace lengths (in meters), computed as distances between
  consecutive contralateral footprints (L–R or R–L) in footprint order.

- `Mean_pace_length`:

  Mean pace length (in meters).

- `Sinuosity`:

  The sinuosity of the track (dimensionless).

- `Straightness`:

  The straightness of the track (dimensionless).

- `Trackway_width`:

  Mean lateral separation between left and right footprints (in meters),
  measured perpendicular to the inferred trackway axis.

- `Gauge`:

  Trackway gauge (dimensionless), computed as
  `Trackway_width / gauge_size`.

- `Pace_angulation`:

  Mean interior angle (in degrees) computed from alternating triplets
  (L–R–L or R–L–R).

- `Step_angle`:

  Mean step angle (in degrees), computed as the mean angle between each
  pace segment and the inferred stride/trackway axis.

## Details

This function calculates various movement parameters for each track in
the provided data.

It uses the following helper functions:

From the trajr package:

- `TrajAngles()`: Calculates the turning angles of the track.

- `TrajDistance()`: Calculates the total distance covered by the track.

- `TrajLength()`: Calculates the length of the track.

- `TrajStepLengths()`: Calculates the step lengths of the track.

- `TrajSinuosity2()`: Calculates the sinuosity of the track.

- `TrajStraightness()`: Calculates the straightness of the track.

From the circular package:

- `circular()`: Converts raw angles (in radians) into a circular data
  type.

- `mean.circular()`: Computes the circular mean of the turning angles.

- `sd.circular()`: Computes the circular standard deviation of the
  turning angles.

The circular mean and circular standard deviation are returned in
degrees in this function.

The reference direction, or 0 degrees, is considered to be along the
positive x-axis. This means that angles are measured counterclockwise
from the positive x-axis, with 0 degrees pointing directly along this
axis. For a detailed explanation and appropriate methods for analyzing
circular data, refer to Batschelet (1981).

Circular mean of turning angles is computed as:

\$\$\overline{\theta} = atan2\left(\frac{1}{n}\sum\_{i=1}^n \sin
\theta_i, \frac{1}{n}\sum\_{i=1}^n \cos \theta_i\right)\$\$

where:

- \\\theta_i\\:

  is the \\i^{th}\\ turning angle in radians.

- \\n\\:

  is the total number of turning angles.

- \\\sin \theta_i\\, \\\cos \theta_i\\:

  are the sine and cosine components of each turning angle.

- \\\mathrm{atan2}(y,x)\\:

  is the two-argument arctangent that returns the angle in the correct
  quadrant.

Circular standard deviation of turning angles is computed as:

\$\$s_c = \sqrt{-2 \ln(\overline{R})}, \quad \overline{R} =
\sqrt{\left(\frac{1}{n}\sum\_{i=1}^n \cos \theta_i\right)^2 +
\left(\frac{1}{n}\sum\_{i=1}^n \sin \theta_i\right)^2}\$\$

where:

- \\\overline{R}\\:

  is the mean resultant length, measuring concentration of angles around
  the mean direction.

- \\n\\:

  is the total number of turning angles.

- \\\cos \theta_i\\, \\\sin \theta_i\\:

  are the trigonometric components of each angle.

- \\s_c\\:

  is the circular standard deviation in radians (converted to degrees in
  this function).

Sinuosity is calculated according to Benhamou (2004), as defined in
equation 8. The formula used here is a refined version of the sinuosity
index presented by Bovet & Benhamou (1988), which is applicable to a
broader range of turning angle distributions and does not require a
constant step length.

The sinuosity is computed using the formula: \$\$S = 2 \left\[ p \left(
\frac{1 + c}{1 - c} + b^2 \right) \right\]^{-0.5}\$\$ where:

- \\p\\:

  is the mean step length (in meters),

- \\c\\:

  is the mean cosine of turning angles (in radians), and

- \\b\\:

  is the coefficient of variation of the step length (in meters).

The straightness index is defined as the ratio D/L, where:

- \\D\\:

  is the beeline distance between the first and last points in the
  trajectory (in meters), and

- \\L\\:

  is the total path length traveled (in meters).

Straightness index is based on the method described by Batschelet
(1981). According to Benhamou (2004), the straightness index serves as a
reliable measure of the efficiency of a directed walk. However, it is
not suitable for random trajectories, as the index for a random walk
tends towards zero with increasing steps. Thus, it is recommended to use
this measure to compare the tortuosity of random walks only if they
consist of a similar number of steps.

## Logo

![](figures/Logo.png)

## References

Batschelet, E. (1981). Circular statistics in biology. Academic press,
111 Fifth Ave., New York, NY 10003, 1981, 388.

Benhamou, S. (2004). How to reliably estimate the tortuosity of an
animal's path:: straightness, sinuosity, or fractal dimension?. Journal
of theoretical biology, 229(2), 209-220.

Bovet, P., & Benhamou, S. (1988). Spatial analysis of animals' movements
using a correlated random walk model. Journal of theoretical biology,
131(4), 419-433.

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
# Example 1:
track_param(PaluxyRiver)
#> [[1]]
#> [[1]]$Turning_angles
#>  [1] 83.36749 77.21330 80.06835 79.46614 74.64986 73.75563 73.91299 73.54942
#>  [9] 72.64598 79.29432 87.32457 83.82411 82.64091 83.85337 82.06724 78.16842
#> [17] 79.91062 84.65694 88.34770 93.28218 98.22952 98.01709 97.67122 96.57457
#> [25] 94.50862 95.42060 96.01985 96.25904
#> 
#> [[1]]$Mean_turning_angle
#> [1] 85.16105
#> 
#> [[1]]$Standard_deviation_turning_angle
#> [1] 8.716076
#> 
#> [[1]]$Beeline_length
#> [1] 16.09043
#> 
#> [[1]]$Path_length
#> [1] 16.2711
#> 
#> [[1]]$Step_lengths
#>  [1] 0.6577971 0.6277174 0.5663694 0.5343256 0.5740178 0.6285614 0.6032098
#>  [8] 0.5749037 0.5822135 0.4440490 0.4650467 0.5851561 0.5931582 0.5676575
#> [15] 0.5347796 0.5611201 0.6195609 0.6061010 0.6775546 0.6066334 0.5307864
#> [22] 0.4669326 0.5366356 0.6446094 0.6075182 0.5974563 0.6002654 0.6769667
#> 
#> [[1]]$Mean_step_length
#> [1] 0.5811108
#> 
#> [[1]]$Standard_deviation_step_length
#> [1] 0.05855971
#> 
#> [[1]]$Stride_length
#>  [1] 1.2554349 1.0686512 1.2571228 1.1498074 0.8880980 1.1703122 1.1353149
#>  [8] 1.1222403 1.2122021 1.2132667 0.9338651 1.2892188 1.1949126 1.3539334
#> [15] 1.3155942 1.1327387 1.1480356 1.2064196 1.1644269 0.9300933 1.1863164
#> [22] 1.0695592 1.2391217 1.3551093 1.0615727 1.0732712 1.2150365 1.2005308
#> 
#> [[1]]$Mean_stride_length
#> [1] 1.162222
#> 
#> [[1]]$Pace_length
#>  [1] 0.8489478 0.7038175 0.7299520 0.5681387 0.6674331 0.6849247 0.7526702
#>  [8] 0.6399838 0.6993039 0.6620884 0.5574720 0.6183847 0.7240928 0.6448396
#> [15] 0.6916613 0.6109176 0.6990343 0.6937972 0.7244181 0.7805154 0.6604209
#> [22] 0.5999591 0.5526160 0.7225684 0.7236762 0.6795511 0.7156487 0.6777459
#> [29] 0.8547561
#> 
#> [[1]]$Mean_pace_length
#> [1] 0.6858392
#> 
#> [[1]]$Sinuosity
#> [1] 0.07738533
#> 
#> [[1]]$Straightness
#> [1] 0.9888963
#> 
#> [[1]]$Trackway_width
#> [1] 0.3594178
#> 
#> [[1]]$Gauge
#> [1] NA
#> 
#> [[1]]$Pace_angulation
#> [1] 117.293
#> 
#> [[1]]$Step_angle
#> [1] 31.38441
#> 
#> 
#> [[2]]
#> [[2]]$Turning_angles
#>  [1]  77.67357  82.52517  77.62784  73.77659  74.96477  72.69947  74.84593
#>  [8]  80.92672  84.04858  84.37223  83.13249  77.30345  73.16693  86.98721
#> [15]  78.55066  91.34262  97.33861  95.27067  91.84761  81.72095  80.07781
#> [22]  84.97592 100.83397
#> 
#> [[2]]$Mean_turning_angle
#> [1] 82.85258
#> 
#> [[2]]$Standard_deviation_turning_angle
#> [1] 7.941837
#> 
#> [[2]]$Beeline_length
#> [1] 15.13227
#> 
#> [[2]]$Path_length
#> [1] 15.27567
#> 
#> [[2]]$Step_lengths
#>  [1] 0.5999317 0.6173971 0.6889308 0.6759706 0.6945371 0.6934507 0.6476978
#>  [8] 0.5781403 0.5652701 0.5976653 0.5990828 0.5629646 0.7121273 0.8673266
#> [15] 0.6998831 0.5558637 0.5778101 0.5907702 0.6732814 0.7085351 0.7558719
#> [22] 0.7931992 0.8199620
#> 
#> [[2]]$Mean_step_length
#> [1] 0.6641595
#> 
#> [[2]]$Standard_deviation_step_length
#> [1] 0.08686601
#> 
#> [[2]]$Stride_length
#>  [1] 1.199863 1.377862 1.389074 1.295396 1.130540 1.198166 1.424255 1.399766
#>  [9] 1.155620 1.346563 1.511744 1.639924 1.234794 1.351941 1.386901 1.156281
#> [17] 1.195331 1.125929 1.734653 1.111727 1.181540 1.417070 1.586398
#> 
#> [[2]]$Mean_stride_length
#> [1] 1.328319
#> 
#> [[2]]$Pace_length
#>  [1] 0.7269895 0.5475275 0.7337508 0.6608203 0.7245222 0.7035362 0.7358799
#>  [8] 0.6320415 0.5652117 0.6241038 0.6248433 0.6195571 0.5950380 0.8797283
#> [15] 0.8555496 0.5766510 0.6583341 0.5675579 0.6707853 0.7478466 0.6816282
#> [22] 0.8471698 0.7816255 0.8594073
#> 
#> [[2]]$Mean_pace_length
#> [1] 0.6925044
#> 
#> [[2]]$Sinuosity
#> [1] 0.1466281
#> 
#> [[2]]$Straightness
#> [1] 0.9906128
#> 
#> [[2]]$Trackway_width
#> [1] 0.0876221
#> 
#> [[2]]$Gauge
#> [1] NA
#> 
#> [[2]]$Pace_angulation
#> [1] 150.9086
#> 
#> [[2]]$Step_angle
#> [1] 15.62501
#> 
#> 

# Example 2:
track_param(MountTom)
#> [[1]]
#> [[1]]$Turning_angles
#> [1] 135.4170 134.9103 139.1806 134.1905 128.4984 139.7955 143.4536 137.0192
#> 
#> [[1]]$Mean_turning_angle
#> [1] 136.5592
#> 
#> [[1]]$Standard_deviation_turning_angle
#> [1] 4.46736
#> 
#> [[1]]$Beeline_length
#> [1] 7.792014
#> 
#> [[1]]$Path_length
#> [1] 7.811916
#> 
#> [[1]]$Step_lengths
#> [1] 1.0716072 0.9967099 1.0056236 0.9936882 0.9780421 0.9329032 0.8594255
#> [8] 0.9739165
#> 
#> [[1]]$Mean_step_length
#> [1] 0.9764895
#> 
#> [[1]]$Standard_deviation_step_length
#> [1] 0.06119478
#> 
#> [[1]]$Stride_length
#> [1] 1.993420 1.987376 1.865806 1.947833 2.143214 2.011247 1.956084 1.718851
#> 
#> [[1]]$Mean_stride_length
#> [1] 1.952979
#> 
#> [[1]]$Pace_length
#> [1] 1.1465582 1.0363986 0.9960140 1.0306737 0.9579627 1.0153771 0.9604282
#> [8] 0.8395926 1.1404306
#> 
#> [[1]]$Mean_pace_length
#> [1] 1.013715
#> 
#> [[1]]$Sinuosity
#> [1] 0.1073126
#> 
#> [[1]]$Straightness
#> [1] 0.9974524
#> 
#> [[1]]$Trackway_width
#> [1] 0.1898719
#> 
#> [[1]]$Gauge
#> [1] NA
#> 
#> [[1]]$Pace_angulation
#> [1] 158.5226
#> 
#> [[1]]$Step_angle
#> [1] 10.62683
#> 
#> 
#> [[2]]
#> [[2]]$Turning_angles
#> [1] 132.9511 134.6430 134.4076 127.5257 121.4985 123.2012 131.3173 125.6057
#> 
#> [[2]]$Mean_turning_angle
#> [1] 128.895
#> 
#> [[2]]$Standard_deviation_turning_angle
#> [1] 5.148772
#> 
#> [[2]]$Beeline_length
#> [1] 8.389291
#> 
#> [[2]]$Path_length
#> [1] 8.419799
#> 
#> [[2]]$Step_lengths
#> [1] 0.9161851 1.0014075 1.0560373 1.0791884 1.2159276 1.0755866 0.9956450
#> [8] 1.0798217
#> 
#> [[2]]$Mean_step_length
#> [1] 1.052475
#> 
#> [[2]]$Standard_deviation_step_length
#> [1] 0.08712676
#> 
#> [[2]]$Stride_length
#> [1] 1.832370 2.112075 2.431855 1.991290 2.002815 2.158377 2.151173 2.159643
#> 
#> [[2]]$Mean_stride_length
#> [1] 2.10495
#> 
#> [[2]]$Pace_length
#> [1] 0.9615624 0.8708180 1.1325182 0.9802660 1.1857229 1.2471023 0.9089200
#> [8] 1.0856887 1.0989901
#> 
#> [[2]]$Mean_pace_length
#> [1] 1.052399
#> 
#> [[2]]$Sinuosity
#> [1] 0.08820136
#> 
#> [[2]]$Straightness
#> [1] 0.9963766
#> 
#> [[2]]$Trackway_width
#> [1] 0.01763514
#> 
#> [[2]]$Gauge
#> [1] NA
#> 
#> [[2]]$Pace_angulation
#> [1] 173.6657
#> 
#> [[2]]$Step_angle
#> [1] 5.892303
#> 
#> 
#> [[3]]
#> [[3]]$Turning_angles
#> [1] 133.8065 133.7811 130.2205 132.0604
#> 
#> [[3]]$Mean_turning_angle
#> [1] 132.4672
#> 
#> [[3]]$Standard_deviation_turning_angle
#> [1] 1.706197
#> 
#> [[3]]$Beeline_length
#> [1] 5.076073
#> 
#> [[3]]$Path_length
#> [1] 5.077853
#> 
#> [[3]]$Step_lengths
#> [1] 1.198182 1.319885 1.404018 1.155769
#> 
#> [[3]]$Mean_step_length
#> [1] 1.269463
#> 
#> [[3]]$Standard_deviation_step_length
#> [1] 0.1135117
#> 
#> [[3]]$Stride_length
#> [1] 2.396364 2.808036 2.639770 2.311538
#> 
#> [[3]]$Mean_stride_length
#> [1] 2.538927
#> 
#> [[3]]$Pace_length
#> [1] 1.296939 1.106465 1.539902 1.268143 1.045104
#> 
#> [[3]]$Mean_pace_length
#> [1] 1.251311
#> 
#> [[3]]$Sinuosity
#> [1] 0.03584756
#> 
#> [[3]]$Straightness
#> [1] 0.9996494
#> 
#> [[3]]$Trackway_width
#> [1] 0.02343501
#> 
#> [[3]]$Gauge
#> [1] NA
#> 
#> [[3]]$Pace_angulation
#> [1] 174.5695
#> 
#> [[3]]$Step_angle
#> [1] 3.035714
#> 
#> 
#> [[4]]
#> [[4]]$Turning_angles
#> [1] 125.9818 121.1350 125.1305 133.1314
#> 
#> [[4]]$Mean_turning_angle
#> [1] 126.3425
#> 
#> [[4]]$Standard_deviation_turning_angle
#> [1] 4.993546
#> 
#> [[4]]$Beeline_length
#> [1] 4.497133
#> 
#> [[4]]$Path_length
#> [1] 4.509806
#> 
#> [[4]]$Step_lengths
#> [1] 1.223971 1.139084 1.046512 1.100240
#> 
#> [[4]]$Mean_step_length
#> [1] 1.127452
#> 
#> [[4]]$Standard_deviation_step_length
#> [1] 0.07470601
#> 
#> [[4]]$Stride_length
#> [1] 2.447941 2.093023 2.278167 2.200480
#> 
#> [[4]]$Mean_stride_length
#> [1] 2.254903
#> 
#> [[4]]$Pace_length
#> [1] 1.226479 1.221637 1.063455 1.054587 1.146117
#> 
#> [[4]]$Mean_pace_length
#> [1] 1.142455
#> 
#> [[4]]$Sinuosity
#> [1] 0.09660109
#> 
#> [[4]]$Straightness
#> [1] 0.9971898
#> 
#> [[4]]$Trackway_width
#> [1] 0.005275705
#> 
#> [[4]]$Gauge
#> [1] NA
#> 
#> [[4]]$Pace_angulation
#> [1] 172.6057
#> 
#> [[4]]$Step_angle
#> [1] 5.141067
#> 
#> 
#> [[5]]
#> [[5]]$Turning_angles
#> [1] 126.6273
#> 
#> [[5]]$Mean_turning_angle
#> [1] 126.6273
#> 
#> [[5]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[5]]$Beeline_length
#> [1] 1.146187
#> 
#> [[5]]$Path_length
#> [1] 1.146187
#> 
#> [[5]]$Step_lengths
#> [1] 1.146187
#> 
#> [[5]]$Mean_step_length
#> [1] 1.146187
#> 
#> [[5]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[5]]$Stride_length
#> [1] 2.292374
#> 
#> [[5]]$Mean_stride_length
#> [1] 2.292374
#> 
#> [[5]]$Pace_length
#> [1] 1.064845 1.233165
#> 
#> [[5]]$Mean_pace_length
#> [1] 1.149005
#> 
#> [[5]]$Sinuosity
#> [1] NaN
#> 
#> [[5]]$Straightness
#> [1] 1
#> 
#> [[5]]$Trackway_width
#> [1] 0.08021087
#> 
#> [[5]]$Gauge
#> [1] NA
#> 
#> [[5]]$Pace_angulation
#> [1] 171.9506
#> 
#> [[5]]$Step_angle
#> [1] 4.024696
#> 
#> 
#> [[6]]
#> [[6]]$Turning_angles
#> [1] 133.6793 134.0478
#> 
#> [[6]]$Mean_turning_angle
#> [1] 133.8635
#> 
#> [[6]]$Standard_deviation_turning_angle
#> [1] 0.2605755
#> 
#> [[6]]$Beeline_length
#> [1] 2.209101
#> 
#> [[6]]$Path_length
#> [1] 2.209112
#> 
#> [[6]]$Step_lengths
#> [1] 1.082785 1.126327
#> 
#> [[6]]$Mean_step_length
#> [1] 1.104556
#> 
#> [[6]]$Standard_deviation_step_length
#> [1] 0.03078894
#> 
#> [[6]]$Stride_length
#> [1] 2.165570 2.252654
#> 
#> [[6]]$Mean_stride_length
#> [1] 2.209112
#> 
#> [[6]]$Pace_length
#> [1] 1.040503 1.130394 1.126444
#> 
#> [[6]]$Mean_pace_length
#> [1] 1.099114
#> 
#> [[6]]$Sinuosity
#> [1] 0.006119749
#> 
#> [[6]]$Straightness
#> [1] 0.9999948
#> 
#> [[6]]$Trackway_width
#> [1] 0.07225275
#> 
#> [[6]]$Gauge
#> [1] NA
#> 
#> [[6]]$Pace_angulation
#> [1] 172.4928
#> 
#> [[6]]$Step_angle
#> [1] 3.784005
#> 
#> 
#> [[7]]
#> [[7]]$Turning_angles
#> [1] 104.5351 109.8321 126.4855 123.7244 127.0905 126.5167
#> 
#> [[7]]$Mean_turning_angle
#> [1] 119.726
#> 
#> [[7]]$Standard_deviation_turning_angle
#> [1] 9.906308
#> 
#> [[7]]$Beeline_length
#> [1] 5.320692
#> 
#> [[7]]$Path_length
#> [1] 5.38372
#> 
#> [[7]]$Step_lengths
#> [1] 0.8430613 0.7152130 0.8424108 1.0210973 1.0314783 0.9304589
#> 
#> [[7]]$Mean_step_length
#> [1] 0.8972866
#> 
#> [[7]]$Standard_deviation_step_length
#> [1] 0.1212761
#> 
#> [[7]]$Stride_length
#> [1] 1.430426 2.042195 1.860918 1.686123 1.684822 2.062957
#> 
#> [[7]]$Mean_stride_length
#> [1] 1.794573
#> 
#> [[7]]$Pace_length
#> [1] 0.9642508 0.7395032 0.7395032 0.9460639 1.1010072 0.9795509 0.9037121
#> 
#> [[7]]$Mean_pace_length
#> [1] 0.9105131
#> 
#> [[7]]$Sinuosity
#> [1] 0.1483819
#> 
#> [[7]]$Straightness
#> [1] 0.9882929
#> 
#> [[7]]$Trackway_width
#> [1] 0.05675529
#> 
#> [[7]]$Gauge
#> [1] NA
#> 
#> [[7]]$Pace_angulation
#> [1] 164.9846
#> 
#> [[7]]$Step_angle
#> [1] 9.280859
#> 
#> 
#> [[8]]
#> [[8]]$Turning_angles
#> [1] 126.1992 120.6832 117.1547 122.0179 123.9117
#> 
#> [[8]]$Mean_turning_angle
#> [1] 121.9937
#> 
#> [[8]]$Standard_deviation_turning_angle
#> [1] 3.409506
#> 
#> [[8]]$Beeline_length
#> [1] 5.956414
#> 
#> [[8]]$Path_length
#> [1] 5.964889
#> 
#> [[8]]$Step_lengths
#> [1] 1.243762 1.236318 1.150321 1.069273 1.265215
#> 
#> [[8]]$Mean_step_length
#> [1] 1.192978
#> 
#> [[8]]$Standard_deviation_step_length
#> [1] 0.08185074
#> 
#> [[8]]$Stride_length
#> [1] 2.472636 2.138546 2.487524 2.300641 2.530431
#> 
#> [[8]]$Mean_stride_length
#> [1] 2.385956
#> 
#> [[8]]$Pace_length
#> [1] 1.272930 1.250889 1.228849 1.096570 1.104405 1.466419
#> 
#> [[8]]$Mean_pace_length
#> [1] 1.236677
#> 
#> [[8]]$Sinuosity
#> [1] 0.06692552
#> 
#> [[8]]$Straightness
#> [1] 0.9985792
#> 
#> [[8]]$Trackway_width
#> [1] 0.2059027
#> 
#> [[8]]$Gauge
#> [1] NA
#> 
#> [[8]]$Pace_angulation
#> [1] 161.4382
#> 
#> [[8]]$Step_angle
#> [1] 9.462611
#> 
#> 
#> [[9]]
#> [[9]]$Turning_angles
#> [1] -26.79077 -20.88580 -20.77225 -18.85063
#> 
#> [[9]]$Mean_turning_angle
#> [1] 338.1763
#> 
#> [[9]]$Standard_deviation_turning_angle
#> [1] 3.439774
#> 
#> [[9]]$Beeline_length
#> [1] 5.452531
#> 
#> [[9]]$Path_length
#> [1] 5.460107
#> 
#> [[9]]$Step_lengths
#> [1] 1.502452 1.435499 1.368359 1.153797
#> 
#> [[9]]$Mean_step_length
#> [1] 1.365027
#> 
#> [[9]]$Standard_deviation_step_length
#> [1] 0.1510864
#> 
#> [[9]]$Stride_length
#> [1] 3.004904 2.736717 2.870998 2.307593
#> 
#> [[9]]$Mean_stride_length
#> [1] 2.730053
#> 
#> [[9]]$Pace_length
#> [1] 1.637710 1.383333 1.487786 1.249122 1.059053
#> 
#> [[9]]$Mean_pace_length
#> [1] 1.363401
#> 
#> [[9]]$Sinuosity
#> [1] 0.0535709
#> 
#> [[9]]$Straightness
#> [1] 0.9986127
#> 
#> [[9]]$Trackway_width
#> [1] 0.06713357
#> 
#> [[9]]$Gauge
#> [1] NA
#> 
#> [[9]]$Pace_angulation
#> [1] 175.7803
#> 
#> [[9]]$Step_angle
#> [1] 3.718311
#> 
#> 
#> [[10]]
#> [[10]]$Turning_angles
#> [1] -59.30028
#> 
#> [[10]]$Mean_turning_angle
#> [1] 300.6997
#> 
#> [[10]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[10]]$Beeline_length
#> [1] 1.149305
#> 
#> [[10]]$Path_length
#> [1] 1.149305
#> 
#> [[10]]$Step_lengths
#> [1] 1.149305
#> 
#> [[10]]$Mean_step_length
#> [1] 1.149305
#> 
#> [[10]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[10]]$Stride_length
#> [1] 2.29861
#> 
#> [[10]]$Mean_stride_length
#> [1] 2.29861
#> 
#> [[10]]$Pace_length
#> [1] 1.263684 1.036464
#> 
#> [[10]]$Mean_pace_length
#> [1] 1.150074
#> 
#> [[10]]$Sinuosity
#> [1] NaN
#> 
#> [[10]]$Straightness
#> [1] 1
#> 
#> [[10]]$Trackway_width
#> [1] 0.04184679
#> 
#> [[10]]$Gauge
#> [1] NA
#> 
#> [[10]]$Pace_angulation
#> [1] 175.7884
#> 
#> [[10]]$Step_angle
#> [1] 2.105807
#> 
#> 
#> [[11]]
#> [[11]]$Turning_angles
#> [1] -65.30517
#> 
#> [[11]]$Mean_turning_angle
#> [1] 294.6948
#> 
#> [[11]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[11]]$Beeline_length
#> [1] 1.209107
#> 
#> [[11]]$Path_length
#> [1] 1.209107
#> 
#> [[11]]$Step_lengths
#> [1] 1.209107
#> 
#> [[11]]$Mean_step_length
#> [1] 1.209107
#> 
#> [[11]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[11]]$Stride_length
#> [1] 2.418215
#> 
#> [[11]]$Mean_stride_length
#> [1] 2.418215
#> 
#> [[11]]$Pace_length
#> [1] 1.205510 1.214254
#> 
#> [[11]]$Mean_pace_length
#> [1] 1.209882
#> 
#> [[11]]$Sinuosity
#> [1] NaN
#> 
#> [[11]]$Straightness
#> [1] 1
#> 
#> [[11]]$Trackway_width
#> [1] 0.04328632
#> 
#> [[11]]$Gauge
#> [1] NA
#> 
#> [[11]]$Pace_angulation
#> [1] 175.8993
#> 
#> [[11]]$Step_angle
#> [1] 2.050353
#> 
#> 
#> [[12]]
#> [[12]]$Turning_angles
#> [1] 111.1961
#> 
#> [[12]]$Mean_turning_angle
#> [1] 111.1961
#> 
#> [[12]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[12]]$Beeline_length
#> [1] 1.201898
#> 
#> [[12]]$Path_length
#> [1] 1.201898
#> 
#> [[12]]$Step_lengths
#> [1] 1.201898
#> 
#> [[12]]$Mean_step_length
#> [1] 1.201898
#> 
#> [[12]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[12]]$Stride_length
#> [1] 2.403797
#> 
#> [[12]]$Mean_stride_length
#> [1] 2.403797
#> 
#> [[12]]$Pace_length
#> [1] 1.102385 1.312083
#> 
#> [[12]]$Mean_pace_length
#> [1] 1.207234
#> 
#> [[12]]$Sinuosity
#> [1] NaN
#> 
#> [[12]]$Straightness
#> [1] 1
#> 
#> [[12]]$Trackway_width
#> [1] 0.1129458
#> 
#> [[12]]$Gauge
#> [1] NA
#> 
#> [[12]]$Pace_angulation
#> [1] 169.1812
#> 
#> [[12]]$Step_angle
#> [1] 5.409406
#> 
#> 
#> [[13]]
#> [[13]]$Turning_angles
#> [1] -134.4543 -139.9441 -128.1402 -116.2039
#> 
#> [[13]]$Mean_turning_angle
#> [1] 230.2982
#> 
#> [[13]]$Standard_deviation_turning_angle
#> [1] 10.2
#> 
#> [[13]]$Beeline_length
#> [1] 1.809277
#> 
#> [[13]]$Path_length
#> [1] 1.831302
#> 
#> [[13]]$Step_lengths
#> [1] 0.4913575 0.4524641 0.4178969 0.4695836
#> 
#> [[13]]$Mean_step_length
#> [1] 0.4578255
#> 
#> [[13]]$Standard_deviation_step_length
#> [1] 0.03101447
#> 
#> [[13]]$Stride_length
#> [1] 0.9827151 0.8357937 0.9049282 0.9391673
#> 
#> [[13]]$Mean_stride_length
#> [1] 0.9156511
#> 
#> [[13]]$Pace_length
#> [1] 0.4212485 0.5618379 0.3535198 0.5277364 0.4131542
#> 
#> [[13]]$Mean_pace_length
#> [1] 0.4554993
#> 
#> [[13]]$Sinuosity
#> [1] 0.2636325
#> 
#> [[13]]$Straightness
#> [1] 0.9879729
#> 
#> [[13]]$Trackway_width
#> [1] 0.02941368
#> 
#> [[13]]$Gauge
#> [1] NA
#> 
#> [[13]]$Pace_angulation
#> [1] 163.576
#> 
#> [[13]]$Step_angle
#> [1] 11.26813
#> 
#> 
#> [[14]]
#> [[14]]$Turning_angles
#> [1] 104.5905 101.1270
#> 
#> [[14]]$Mean_turning_angle
#> [1] 102.8588
#> 
#> [[14]]$Standard_deviation_turning_angle
#> [1] 2.449041
#> 
#> [[14]]$Beeline_length
#> [1] 1.998605
#> 
#> [[14]]$Path_length
#> [1] 1.999516
#> 
#> [[14]]$Step_lengths
#> [1] 1.0507992 0.9487167
#> 
#> [[14]]$Mean_step_length
#> [1] 0.9997579
#> 
#> [[14]]$Standard_deviation_step_length
#> [1] 0.0721832
#> 
#> [[14]]$Stride_length
#> [1] 1.897433 2.101598
#> 
#> [[14]]$Mean_stride_length
#> [1] 1.999516
#> 
#> [[14]]$Pace_length
#> [1] 1.2329758 0.8936302 1.0158945
#> 
#> [[14]]$Mean_pace_length
#> [1] 1.0475
#> 
#> [[14]]$Sinuosity
#> [1] 0.06047449
#> 
#> [[14]]$Straightness
#> [1] 0.9995445
#> 
#> [[14]]$Trackway_width
#> [1] 0.1346045
#> 
#> [[14]]$Gauge
#> [1] NA
#> 
#> [[14]]$Pace_angulation
#> [1] 164.6249
#> 
#> [[14]]$Step_angle
#> [1] 7.34011
#> 
#> 
#> [[15]]
#> [[15]]$Turning_angles
#> [1] 114.2659 111.3107 118.9163 120.6646
#> 
#> [[15]]$Mean_turning_angle
#> [1] 116.2898
#> 
#> [[15]]$Standard_deviation_turning_angle
#> [1] 4.278817
#> 
#> [[15]]$Beeline_length
#> [1] 4.556698
#> 
#> [[15]]$Path_length
#> [1] 4.56646
#> 
#> [[15]]$Step_lengths
#> [1] 1.057393 1.195731 1.154191 1.159145
#> 
#> [[15]]$Mean_step_length
#> [1] 1.141615
#> 
#> [[15]]$Standard_deviation_step_length
#> [1] 0.05912482
#> 
#> [[15]]$Stride_length
#> [1] 2.391462 2.318289 2.114787 2.308382
#> 
#> [[15]]$Mean_stride_length
#> [1] 2.28323
#> 
#> [[15]]$Pace_length
#> [1] 0.9269852 1.1959790 1.1962719 1.1262364 1.1996107
#> 
#> [[15]]$Mean_pace_length
#> [1] 1.129017
#> 
#> [[15]]$Sinuosity
#> [1] 0.07871941
#> 
#> [[15]]$Straightness
#> [1] 0.9978623
#> 
#> [[15]]$Trackway_width
#> [1] 0.0005184741
#> 
#> [[15]]$Gauge
#> [1] NA
#> 
#> [[15]]$Pace_angulation
#> [1] 171.2543
#> 
#> [[15]]$Step_angle
#> [1] 4.587023
#> 
#> 
#> [[16]]
#> [[16]]$Turning_angles
#> [1] 116.5651 122.1172 120.0741 121.6150 129.5226 129.5917 135.7970
#> 
#> [[16]]$Mean_turning_angle
#> [1] 125.0358
#> 
#> [[16]]$Standard_deviation_turning_angle
#> [1] 6.748851
#> 
#> [[16]]$Beeline_length
#> [1] 7.068614
#> 
#> [[16]]$Path_length
#> [1] 7.109877
#> 
#> [[16]]$Step_lengths
#> [1] 0.9075806 0.9584426 1.0476549 1.0982925 1.0294856 1.0591361 1.0092847
#> 
#> [[16]]$Mean_step_length
#> [1] 1.015697
#> 
#> [[16]]$Standard_deviation_step_length
#> [1] 0.06445766
#> 
#> [[16]]$Stride_length
#> [1] 1.815161 2.095310 2.058971 2.018569 1.916885 2.196585 2.118272
#> 
#> [[16]]$Mean_stride_length
#> [1] 2.031393
#> 
#> [[16]]$Pace_length
#> [1] 0.9103429 0.9142684 1.0026361 1.0936195 1.1063331 0.9608334 1.1652406
#> [8] 0.8996648
#> 
#> [[16]]$Mean_pace_length
#> [1] 1.006617
#> 
#> [[16]]$Sinuosity
#> [1] 0.08320294
#> 
#> [[16]]$Straightness
#> [1] 0.9941965
#> 
#> [[16]]$Trackway_width
#> [1] 0.1388154
#> 
#> [[16]]$Gauge
#> [1] NA
#> 
#> [[16]]$Pace_angulation
#> [1] 170.4862
#> 
#> [[16]]$Step_angle
#> [1] 7.654888
#> 
#> 
#> [[17]]
#> [[17]]$Turning_angles
#> [1] -131.2245
#> 
#> [[17]]$Mean_turning_angle
#> [1] 228.7755
#> 
#> [[17]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[17]]$Beeline_length
#> [1] 0.3079478
#> 
#> [[17]]$Path_length
#> [1] 0.3079478
#> 
#> [[17]]$Step_lengths
#> [1] 0.3079478
#> 
#> [[17]]$Mean_step_length
#> [1] 0.3079478
#> 
#> [[17]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[17]]$Stride_length
#> [1] 0.6158957
#> 
#> [[17]]$Mean_stride_length
#> [1] 0.6158957
#> 
#> [[17]]$Pace_length
#> [1] 0.3010686 0.3184427
#> 
#> [[17]]$Mean_pace_length
#> [1] 0.3097556
#> 
#> [[17]]$Sinuosity
#> [1] NaN
#> 
#> [[17]]$Straightness
#> [1] 1
#> 
#> [[17]]$Trackway_width
#> [1] 0.03340355
#> 
#> [[17]]$Gauge
#> [1] NA
#> 
#> [[17]]$Pace_angulation
#> [1] 167.6087
#> 
#> [[17]]$Step_angle
#> [1] 6.195643
#> 
#> 
#> [[18]]
#> [[18]]$Turning_angles
#> [1] 107.7676 110.0561 113.7495 110.8756
#> 
#> [[18]]$Mean_turning_angle
#> [1] 110.6121
#> 
#> [[18]]$Standard_deviation_turning_angle
#> [1] 2.470679
#> 
#> [[18]]$Beeline_length
#> [1] 4.767794
#> 
#> [[18]]$Path_length
#> [1] 4.771307
#> 
#> [[18]]$Step_lengths
#> [1] 1.257788 1.183538 1.265235 1.064747
#> 
#> [[18]]$Mean_step_length
#> [1] 1.192827
#> 
#> [[18]]$Standard_deviation_step_length
#> [1] 0.09301194
#> 
#> [[18]]$Stride_length
#> [1] 2.367075 2.129494 2.515577 2.530469
#> 
#> [[18]]$Mean_stride_length
#> [1] 2.385654
#> 
#> [[18]]$Pace_length
#> [1] 1.3872949 1.1555288 1.2559978 1.2922481 0.8444469
#> 
#> [[18]]$Mean_pace_length
#> [1] 1.187103
#> 
#> [[18]]$Sinuosity
#> [1] 0.04807342
#> 
#> [[18]]$Straightness
#> [1] 0.9992637
#> 
#> [[18]]$Trackway_width
#> [1] 0.1359546
#> 
#> [[18]]$Gauge
#> [1] NA
#> 
#> [[18]]$Pace_angulation
#> [1] 164.4785
#> 
#> [[18]]$Step_angle
#> [1] 7.24451
#> 
#> 
#> [[19]]
#> [[19]]$Turning_angles
#> [1] 111.8959
#> 
#> [[19]]$Mean_turning_angle
#> [1] 111.8959
#> 
#> [[19]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[19]]$Beeline_length
#> [1] 1.490617
#> 
#> [[19]]$Path_length
#> [1] 1.490617
#> 
#> [[19]]$Step_lengths
#> [1] 1.490617
#> 
#> [[19]]$Mean_step_length
#> [1] 1.490617
#> 
#> [[19]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[19]]$Stride_length
#> [1] 2.981234
#> 
#> [[19]]$Mean_stride_length
#> [1] 2.981234
#> 
#> [[19]]$Pace_length
#> [1] 1.731892 1.249348
#> 
#> [[19]]$Mean_pace_length
#> [1] 1.49062
#> 
#> [[19]]$Sinuosity
#> [1] NaN
#> 
#> [[19]]$Straightness
#> [1] 1
#> 
#> [[19]]$Trackway_width
#> [1] 0.00287917
#> 
#> [[19]]$Gauge
#> [1] NA
#> 
#> [[19]]$Pace_angulation
#> [1] 179.7727
#> 
#> [[19]]$Step_angle
#> [1] 0.1136457
#> 
#> 
#> [[20]]
#> [[20]]$Turning_angles
#> [1] 104.3162 108.0343
#> 
#> [[20]]$Mean_turning_angle
#> [1] 106.1753
#> 
#> [[20]]$Standard_deviation_turning_angle
#> [1] 2.629063
#> 
#> [[20]]$Beeline_length
#> [1] 2.300393
#> 
#> [[20]]$Path_length
#> [1] 2.301602
#> 
#> [[20]]$Step_lengths
#> [1] 1.204311 1.097291
#> 
#> [[20]]$Mean_step_length
#> [1] 1.150801
#> 
#> [[20]]$Standard_deviation_step_length
#> [1] 0.07567463
#> 
#> [[20]]$Stride_length
#> [1] 2.194581 2.408622
#> 
#> [[20]]$Mean_stride_length
#> [1] 2.301602
#> 
#> [[20]]$Pace_length
#> [1] 1.237467 1.172650 1.022968
#> 
#> [[20]]$Mean_pace_length
#> [1] 1.144362
#> 
#> [[20]]$Sinuosity
#> [1] 0.06051244
#> 
#> [[20]]$Straightness
#> [1] 0.9994748
#> 
#> [[20]]$Trackway_width
#> [1] 0.005588892
#> 
#> [[20]]$Gauge
#> [1] NA
#> 
#> [[20]]$Pace_angulation
#> [1] 176.2159
#> 
#> [[20]]$Step_angle
#> [1] 2.632583
#> 
#> 
#> [[21]]
#> [[21]]$Turning_angles
#> [1] 140.7497 134.7000
#> 
#> [[21]]$Mean_turning_angle
#> [1] 137.7249
#> 
#> [[21]]$Standard_deviation_turning_angle
#> [1] 4.277776
#> 
#> [[21]]$Beeline_length
#> [1] 2.325026
#> 
#> [[21]]$Path_length
#> [1] 2.328268
#> 
#> [[21]]$Step_lengths
#> [1] 1.136569 1.191699
#> 
#> [[21]]$Mean_step_length
#> [1] 1.164134
#> 
#> [[21]]$Standard_deviation_step_length
#> [1] 0.03898298
#> 
#> [[21]]$Stride_length
#> [1] 2.383399 2.273138
#> 
#> [[21]]$Mean_stride_length
#> [1] 2.328268
#> 
#> [[21]]$Pace_length
#> [1] 0.9802163 1.2987085 1.0874083
#> 
#> [[21]]$Mean_pace_length
#> [1] 1.122111
#> 
#> [[21]]$Sinuosity
#> [1] 0.09795169
#> 
#> [[21]]$Straightness
#> [1] 0.9986075
#> 
#> [[21]]$Trackway_width
#> [1] 0.01064242
#> 
#> [[21]]$Gauge
#> [1] NA
#> 
#> [[21]]$Pace_angulation
#> [1] 173.1296
#> 
#> [[21]]$Step_angle
#> [1] 4.722496
#> 
#> 
#> [[22]]
#> [[22]]$Turning_angles
#> [1] 118.3636 109.6538
#> 
#> [[22]]$Mean_turning_angle
#> [1] 114.0087
#> 
#> [[22]]$Standard_deviation_turning_angle
#> [1] 6.158744
#> 
#> [[22]]$Beeline_length
#> [1] 1.795854
#> 
#> [[22]]$Path_length
#> [1] 1.801009
#> 
#> [[22]]$Step_lengths
#> [1] 0.8172251 0.9837840
#> 
#> [[22]]$Mean_step_length
#> [1] 0.9005046
#> 
#> [[22]]$Standard_deviation_step_length
#> [1] 0.1177749
#> 
#> [[22]]$Stride_length
#> [1] 1.967568 1.634450
#> 
#> [[22]]$Mean_stride_length
#> [1] 1.801009
#> 
#> [[22]]$Pace_length
#> [1] 0.7882257 0.9617344 1.0348481
#> 
#> [[22]]$Mean_pace_length
#> [1] 0.9282694
#> 
#> [[22]]$Sinuosity
#> [1] 0.1604935
#> 
#> [[22]]$Straightness
#> [1] 0.9971376
#> 
#> [[22]]$Trackway_width
#> [1] 0.2363052
#> 
#> [[22]]$Gauge
#> [1] NA
#> 
#> [[22]]$Pace_angulation
#> [1] 149.1711
#> 
#> [[22]]$Step_angle
#> [1] 15.8974
#> 
#> 
#> [[23]]
#> [[23]]$Turning_angles
#> [1] -9.811257
#> 
#> [[23]]$Mean_turning_angle
#> [1] 350.1887
#> 
#> [[23]]$Standard_deviation_turning_angle
#> [1] NA
#> 
#> [[23]]$Beeline_length
#> [1] 1.786422
#> 
#> [[23]]$Path_length
#> [1] 1.786422
#> 
#> [[23]]$Step_lengths
#> [1] 1.786422
#> 
#> [[23]]$Mean_step_length
#> [1] 1.786422
#> 
#> [[23]]$Standard_deviation_step_length
#> [1] NA
#> 
#> [[23]]$Stride_length
#> [1] 3.572843
#> 
#> [[23]]$Mean_stride_length
#> [1] 3.572843
#> 
#> [[23]]$Pace_length
#> [1] 2.053557 1.539428
#> 
#> [[23]]$Mean_pace_length
#> [1] 1.796492
#> 
#> [[23]]$Sinuosity
#> [1] NaN
#> 
#> [[23]]$Straightness
#> [1] 1
#> 
#> [[23]]$Trackway_width
#> [1] 0.1879773
#> 
#> [[23]]$Gauge
#> [1] NA
#> 
#> [[23]]$Pace_angulation
#> [1] 167.7341
#> 
#> [[23]]$Step_angle
#> [1] 6.132934
#> 
#> 
```
