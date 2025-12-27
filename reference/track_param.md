# Print track parameters

`track_param()` is a function to compute and print various parameters of
tracks from a list of track data.

## Usage

``` r
track_param(data)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

## Value

A list of lists, where each sublist contains the computed parameters for
a corresponding track. The parameters included are:

- `Turning_angles`: A vector of turning angles for the track (in
  degrees).

- `Mean_turning_angle`: The mean of the turning angles (in degrees).

- `Standard_deviation_turning_angle`: The standard deviation of the
  turning angles (in degrees).

- `Distance`: The total distance covered by the track (in meters).

- `Length`: The length of the track (in meters).

- `Step_lengths`: A vector of step lengths for the track (in meters).

- `Mean_step_length`: The mean of the step lengths (in meters).

- `Standard_deviation_step_length`: The standard deviation of the step
  lengths (in meters).

- `Sinuosity`: The sinuosity of the track (dimensionless).

- `Straightness`: The straightness of the track (dimensionless).

- `Trackway_width`: Mean lateral separation between left and right
  footprints (in meters), measured perpendicular to the inferred
  trackway axis.

- `Pace_angulation`: Mean interior angle (in degrees) computed from
  alternating triplets (L–R–L or R–L–R).

The reference direction, or 0 degrees, is considered to be along the
positive x-axis. This means that angles are measured counterclockwise
from the positive x-axis, with 0 degrees (or 0 degrees) pointing
directly along this axis. For a detailed explanation and appropriate
methods for analyzing circular data, refer to Batschelet (1981).

Circular mean of turning angles is computed as:

\$\$\bar{θ} = atan2\left(\frac{1}{n}\sum\_{i=1}^n \sin θ_i,
\frac{1}{n}\sum\_{i=1}^n \cos θ_i\right)\$\$

where:

- \\θ_i\\ is the \\i^{th}\\ turning angle in radians.

- \\n\\ is the total number of turning angles.

- \\\sin θ_i\\, \\\cos θ_i\\ are the sine and cosine components of each
  turning angle.

- \\atan2(y, x)\\ is the two-argument arctangent that returns the angle
  in the correct quadrant.

Circular standard deviation of turning angles is computed as:

\$\$s_c = \sqrt{-2 \ln(\bar R)}, \quad \bar R =
\sqrt{\left(\frac{1}{n}\sum\_{i=1}^n \cos θ_i\right)^2 +
\left(\frac{1}{n}\sum\_{i=1}^n \sin θ_i\right)^2}\$\$

where:

- \\\bar R\\ is the mean resultant length, measuring concentration of
  angles around the mean direction.

- \\n\\ is the total number of turning angles.

- \\\cos θ_i\\, \\\sin θ_i\\ are the trigonometric components of each
  angle.

- \\s_c\\ is the circular standard deviation in radians (converted to
  degrees in this function).

Sinuosity is calculated according to Benhamou (2004), as defined in
equation 8. The formula used here is a refined version of the sinuosity
index presented by Bovet & Benhamou (1988), which is applicable to a
broader range of turning angle distributions and does not require a
constant step length.

The sinuosity is computed using the formula: \$\$S = 2 \left\[ p \left(
\frac{1 + c}{1 - c} + b^2 \right) \right\]^{-0.5}\$\$ where:

- \\p\\ is the mean step length (in meters),

- \\c\\ is the mean cosine of turning angles (in radians), and

- \\b\\ is the coefficient of variation of the step length (in meters).

The straightness index is defined as the ratio D/L, where:

- \\D\\ is the beeline distance between the first and last points in the
  trajectory (in meters), and

- \\L\\ is the total path length traveled (in meters).

Straightness index is based on the method described by Batschelet
(1981). According to Benhamou (2004), the straightness index serves as a
reliable measure of the efficiency of a directed walk. However, it is
not suitable for random trajectories, as the index for a random walk
tends towards zero with increasing steps. Thus, it is recommended to use
this measure to compare the tortuosity of random walks only if they
consist of a similar number of steps.

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
#> $Track_1
#> $Track_1$Turning_angles
#>  [1] 83.36749 77.21330 80.06835 79.46614 74.64986 73.75563 73.91299 73.54942
#>  [9] 72.64598 79.29432 87.32457 83.82411 82.64091 83.85337 82.06724 78.16842
#> [17] 79.91062 84.65694 88.34770 93.28218 98.22952 98.01709 97.67122 96.57457
#> [25] 94.50862 95.42060 96.01985 96.25904
#> 
#> $Track_1$Mean_turning_angle
#> [1] 85.16105
#> 
#> $Track_1$Standard_deviation_turning_angle
#> [1] 8.716076
#> 
#> $Track_1$Distance
#> [1] 16.09043
#> 
#> $Track_1$Length
#> [1] 16.2711
#> 
#> $Track_1$Step_lengths
#>  [1] 0.6577971 0.6277174 0.5663694 0.5343256 0.5740178 0.6285614 0.6032098
#>  [8] 0.5749037 0.5822135 0.4440490 0.4650467 0.5851561 0.5931582 0.5676575
#> [15] 0.5347796 0.5611201 0.6195609 0.6061010 0.6775546 0.6066334 0.5307864
#> [22] 0.4669326 0.5366356 0.6446094 0.6075182 0.5974563 0.6002654 0.6769667
#> 
#> $Track_1$Mean_step_length
#> [1] 0.5811108
#> 
#> $Track_1$Standard_deviation_step_length
#> [1] 0.05855971
#> 
#> $Track_1$Sinuosity
#> [1] 0.07738533
#> 
#> $Track_1$Straightness
#> [1] 0.9888963
#> 
#> $Track_1$Trackway_width
#> [1] 0.3594178
#> 
#> $Track_1$Pace_angulation
#> [1] 117.293
#> 
#> 
#> $Track_2
#> $Track_2$Turning_angles
#>  [1]  77.67357  82.52517  77.62784  73.77659  74.96477  72.69947  74.84593
#>  [8]  80.92672  84.04858  84.37223  83.13249  77.30345  73.16693  86.98721
#> [15]  78.55066  91.34262  97.33861  95.27067  91.84761  81.72095  80.07781
#> [22]  84.97592 100.83397
#> 
#> $Track_2$Mean_turning_angle
#> [1] 82.85258
#> 
#> $Track_2$Standard_deviation_turning_angle
#> [1] 7.941837
#> 
#> $Track_2$Distance
#> [1] 15.13227
#> 
#> $Track_2$Length
#> [1] 15.27567
#> 
#> $Track_2$Step_lengths
#>  [1] 0.5999317 0.6173971 0.6889308 0.6759706 0.6945371 0.6934507 0.6476978
#>  [8] 0.5781403 0.5652701 0.5976653 0.5990828 0.5629646 0.7121273 0.8673266
#> [15] 0.6998831 0.5558637 0.5778101 0.5907702 0.6732814 0.7085351 0.7558719
#> [22] 0.7931992 0.8199620
#> 
#> $Track_2$Mean_step_length
#> [1] 0.6641595
#> 
#> $Track_2$Standard_deviation_step_length
#> [1] 0.08686601
#> 
#> $Track_2$Sinuosity
#> [1] 0.1466281
#> 
#> $Track_2$Straightness
#> [1] 0.9906128
#> 
#> $Track_2$Trackway_width
#> [1] 0.0876221
#> 
#> $Track_2$Pace_angulation
#> [1] 150.9086
#> 
#> 

# Example 2:
track_param(MountTom)
#> $Track_01
#> $Track_01$Turning_angles
#> [1] 135.4170 134.9103 139.1806 134.1905 128.4984 139.7955 143.4536 137.0192
#> 
#> $Track_01$Mean_turning_angle
#> [1] 136.5592
#> 
#> $Track_01$Standard_deviation_turning_angle
#> [1] 4.46736
#> 
#> $Track_01$Distance
#> [1] 7.792014
#> 
#> $Track_01$Length
#> [1] 7.811916
#> 
#> $Track_01$Step_lengths
#> [1] 1.0716072 0.9967099 1.0056236 0.9936882 0.9780421 0.9329032 0.8594255
#> [8] 0.9739165
#> 
#> $Track_01$Mean_step_length
#> [1] 0.9764895
#> 
#> $Track_01$Standard_deviation_step_length
#> [1] 0.06119478
#> 
#> $Track_01$Sinuosity
#> [1] 0.1073126
#> 
#> $Track_01$Straightness
#> [1] 0.9974524
#> 
#> $Track_01$Trackway_width
#> [1] 0.1898719
#> 
#> $Track_01$Pace_angulation
#> [1] 158.5226
#> 
#> 
#> $Track_02
#> $Track_02$Turning_angles
#> [1] 132.9511 134.6430 134.4076 127.5257 121.4985 123.2012 131.3173 125.6057
#> 
#> $Track_02$Mean_turning_angle
#> [1] 128.895
#> 
#> $Track_02$Standard_deviation_turning_angle
#> [1] 5.148772
#> 
#> $Track_02$Distance
#> [1] 8.389291
#> 
#> $Track_02$Length
#> [1] 8.419799
#> 
#> $Track_02$Step_lengths
#> [1] 0.9161851 1.0014075 1.0560373 1.0791884 1.2159276 1.0755866 0.9956450
#> [8] 1.0798217
#> 
#> $Track_02$Mean_step_length
#> [1] 1.052475
#> 
#> $Track_02$Standard_deviation_step_length
#> [1] 0.08712676
#> 
#> $Track_02$Sinuosity
#> [1] 0.08820136
#> 
#> $Track_02$Straightness
#> [1] 0.9963766
#> 
#> $Track_02$Trackway_width
#> [1] 0.01763514
#> 
#> $Track_02$Pace_angulation
#> [1] 173.6657
#> 
#> 
#> $Track_03
#> $Track_03$Turning_angles
#> [1] 133.8065 133.7811 130.2205 132.0604
#> 
#> $Track_03$Mean_turning_angle
#> [1] 132.4672
#> 
#> $Track_03$Standard_deviation_turning_angle
#> [1] 1.706197
#> 
#> $Track_03$Distance
#> [1] 5.076073
#> 
#> $Track_03$Length
#> [1] 5.077853
#> 
#> $Track_03$Step_lengths
#> [1] 1.198182 1.319885 1.404018 1.155769
#> 
#> $Track_03$Mean_step_length
#> [1] 1.269463
#> 
#> $Track_03$Standard_deviation_step_length
#> [1] 0.1135117
#> 
#> $Track_03$Sinuosity
#> [1] 0.03584756
#> 
#> $Track_03$Straightness
#> [1] 0.9996494
#> 
#> $Track_03$Trackway_width
#> [1] 0.02343501
#> 
#> $Track_03$Pace_angulation
#> [1] 174.5695
#> 
#> 
#> $Track_04
#> $Track_04$Turning_angles
#> [1] 125.9818 121.1350 125.1305 133.1314
#> 
#> $Track_04$Mean_turning_angle
#> [1] 126.3425
#> 
#> $Track_04$Standard_deviation_turning_angle
#> [1] 4.993546
#> 
#> $Track_04$Distance
#> [1] 4.497133
#> 
#> $Track_04$Length
#> [1] 4.509806
#> 
#> $Track_04$Step_lengths
#> [1] 1.223971 1.139084 1.046512 1.100240
#> 
#> $Track_04$Mean_step_length
#> [1] 1.127452
#> 
#> $Track_04$Standard_deviation_step_length
#> [1] 0.07470601
#> 
#> $Track_04$Sinuosity
#> [1] 0.09660109
#> 
#> $Track_04$Straightness
#> [1] 0.9971898
#> 
#> $Track_04$Trackway_width
#> [1] 0.005275705
#> 
#> $Track_04$Pace_angulation
#> [1] 172.6057
#> 
#> 
#> $Track_05
#> $Track_05$Turning_angles
#> [1] 126.6273
#> 
#> $Track_05$Mean_turning_angle
#> [1] 126.6273
#> 
#> $Track_05$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_05$Distance
#> [1] 1.146187
#> 
#> $Track_05$Length
#> [1] 1.146187
#> 
#> $Track_05$Step_lengths
#> [1] 1.146187
#> 
#> $Track_05$Mean_step_length
#> [1] 1.146187
#> 
#> $Track_05$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_05$Sinuosity
#> [1] NaN
#> 
#> $Track_05$Straightness
#> [1] 1
#> 
#> $Track_05$Trackway_width
#> [1] 0.08021087
#> 
#> $Track_05$Pace_angulation
#> [1] 171.9506
#> 
#> 
#> $Track_06
#> $Track_06$Turning_angles
#> [1] 133.6793 134.0478
#> 
#> $Track_06$Mean_turning_angle
#> [1] 133.8635
#> 
#> $Track_06$Standard_deviation_turning_angle
#> [1] 0.2605755
#> 
#> $Track_06$Distance
#> [1] 2.209101
#> 
#> $Track_06$Length
#> [1] 2.209112
#> 
#> $Track_06$Step_lengths
#> [1] 1.082785 1.126327
#> 
#> $Track_06$Mean_step_length
#> [1] 1.104556
#> 
#> $Track_06$Standard_deviation_step_length
#> [1] 0.03078894
#> 
#> $Track_06$Sinuosity
#> [1] 0.006119749
#> 
#> $Track_06$Straightness
#> [1] 0.9999948
#> 
#> $Track_06$Trackway_width
#> [1] 0.07225275
#> 
#> $Track_06$Pace_angulation
#> [1] 172.4928
#> 
#> 
#> $Track_07
#> $Track_07$Turning_angles
#> [1] 104.5351 109.8321 126.4855 123.7244 127.0905 126.5167
#> 
#> $Track_07$Mean_turning_angle
#> [1] 119.726
#> 
#> $Track_07$Standard_deviation_turning_angle
#> [1] 9.906308
#> 
#> $Track_07$Distance
#> [1] 5.320692
#> 
#> $Track_07$Length
#> [1] 5.38372
#> 
#> $Track_07$Step_lengths
#> [1] 0.8430613 0.7152130 0.8424108 1.0210973 1.0314783 0.9304589
#> 
#> $Track_07$Mean_step_length
#> [1] 0.8972866
#> 
#> $Track_07$Standard_deviation_step_length
#> [1] 0.1212761
#> 
#> $Track_07$Sinuosity
#> [1] 0.1483819
#> 
#> $Track_07$Straightness
#> [1] 0.9882929
#> 
#> $Track_07$Trackway_width
#> [1] 0.05675529
#> 
#> $Track_07$Pace_angulation
#> [1] 164.9846
#> 
#> 
#> $Track_08
#> $Track_08$Turning_angles
#> [1] 126.1992 120.6832 117.1547 122.0179 123.9117
#> 
#> $Track_08$Mean_turning_angle
#> [1] 121.9937
#> 
#> $Track_08$Standard_deviation_turning_angle
#> [1] 3.409506
#> 
#> $Track_08$Distance
#> [1] 5.956414
#> 
#> $Track_08$Length
#> [1] 5.964889
#> 
#> $Track_08$Step_lengths
#> [1] 1.243762 1.236318 1.150321 1.069273 1.265215
#> 
#> $Track_08$Mean_step_length
#> [1] 1.192978
#> 
#> $Track_08$Standard_deviation_step_length
#> [1] 0.08185074
#> 
#> $Track_08$Sinuosity
#> [1] 0.06692552
#> 
#> $Track_08$Straightness
#> [1] 0.9985792
#> 
#> $Track_08$Trackway_width
#> [1] 0.2059027
#> 
#> $Track_08$Pace_angulation
#> [1] 161.4382
#> 
#> 
#> $Track_09
#> $Track_09$Turning_angles
#> [1] -26.79077 -20.88580 -20.77225 -18.85063
#> 
#> $Track_09$Mean_turning_angle
#> [1] 338.1763
#> 
#> $Track_09$Standard_deviation_turning_angle
#> [1] 3.439774
#> 
#> $Track_09$Distance
#> [1] 5.452531
#> 
#> $Track_09$Length
#> [1] 5.460107
#> 
#> $Track_09$Step_lengths
#> [1] 1.502452 1.435499 1.368359 1.153797
#> 
#> $Track_09$Mean_step_length
#> [1] 1.365027
#> 
#> $Track_09$Standard_deviation_step_length
#> [1] 0.1510864
#> 
#> $Track_09$Sinuosity
#> [1] 0.0535709
#> 
#> $Track_09$Straightness
#> [1] 0.9986127
#> 
#> $Track_09$Trackway_width
#> [1] 0.06713357
#> 
#> $Track_09$Pace_angulation
#> [1] 175.7803
#> 
#> 
#> $Track_10
#> $Track_10$Turning_angles
#> [1] -59.30028
#> 
#> $Track_10$Mean_turning_angle
#> [1] 300.6997
#> 
#> $Track_10$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_10$Distance
#> [1] 1.149305
#> 
#> $Track_10$Length
#> [1] 1.149305
#> 
#> $Track_10$Step_lengths
#> [1] 1.149305
#> 
#> $Track_10$Mean_step_length
#> [1] 1.149305
#> 
#> $Track_10$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_10$Sinuosity
#> [1] NaN
#> 
#> $Track_10$Straightness
#> [1] 1
#> 
#> $Track_10$Trackway_width
#> [1] 0.04184679
#> 
#> $Track_10$Pace_angulation
#> [1] 175.7884
#> 
#> 
#> $Track_11
#> $Track_11$Turning_angles
#> [1] -65.30517
#> 
#> $Track_11$Mean_turning_angle
#> [1] 294.6948
#> 
#> $Track_11$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_11$Distance
#> [1] 1.209107
#> 
#> $Track_11$Length
#> [1] 1.209107
#> 
#> $Track_11$Step_lengths
#> [1] 1.209107
#> 
#> $Track_11$Mean_step_length
#> [1] 1.209107
#> 
#> $Track_11$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_11$Sinuosity
#> [1] NaN
#> 
#> $Track_11$Straightness
#> [1] 1
#> 
#> $Track_11$Trackway_width
#> [1] 0.04328632
#> 
#> $Track_11$Pace_angulation
#> [1] 175.8993
#> 
#> 
#> $Track_12
#> $Track_12$Turning_angles
#> [1] 111.1961
#> 
#> $Track_12$Mean_turning_angle
#> [1] 111.1961
#> 
#> $Track_12$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_12$Distance
#> [1] 1.201898
#> 
#> $Track_12$Length
#> [1] 1.201898
#> 
#> $Track_12$Step_lengths
#> [1] 1.201898
#> 
#> $Track_12$Mean_step_length
#> [1] 1.201898
#> 
#> $Track_12$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_12$Sinuosity
#> [1] NaN
#> 
#> $Track_12$Straightness
#> [1] 1
#> 
#> $Track_12$Trackway_width
#> [1] 0.1129458
#> 
#> $Track_12$Pace_angulation
#> [1] 169.1812
#> 
#> 
#> $Track_13
#> $Track_13$Turning_angles
#> [1] -134.4543 -139.9441 -128.1402 -116.2039
#> 
#> $Track_13$Mean_turning_angle
#> [1] 230.2982
#> 
#> $Track_13$Standard_deviation_turning_angle
#> [1] 10.2
#> 
#> $Track_13$Distance
#> [1] 1.809277
#> 
#> $Track_13$Length
#> [1] 1.831302
#> 
#> $Track_13$Step_lengths
#> [1] 0.4913575 0.4524641 0.4178969 0.4695836
#> 
#> $Track_13$Mean_step_length
#> [1] 0.4578255
#> 
#> $Track_13$Standard_deviation_step_length
#> [1] 0.03101447
#> 
#> $Track_13$Sinuosity
#> [1] 0.2636325
#> 
#> $Track_13$Straightness
#> [1] 0.9879729
#> 
#> $Track_13$Trackway_width
#> [1] 0.02941368
#> 
#> $Track_13$Pace_angulation
#> [1] 163.576
#> 
#> 
#> $Track_14
#> $Track_14$Turning_angles
#> [1] 104.5905 101.1270
#> 
#> $Track_14$Mean_turning_angle
#> [1] 102.8588
#> 
#> $Track_14$Standard_deviation_turning_angle
#> [1] 2.449041
#> 
#> $Track_14$Distance
#> [1] 1.998605
#> 
#> $Track_14$Length
#> [1] 1.999516
#> 
#> $Track_14$Step_lengths
#> [1] 1.0507992 0.9487167
#> 
#> $Track_14$Mean_step_length
#> [1] 0.9997579
#> 
#> $Track_14$Standard_deviation_step_length
#> [1] 0.0721832
#> 
#> $Track_14$Sinuosity
#> [1] 0.06047449
#> 
#> $Track_14$Straightness
#> [1] 0.9995445
#> 
#> $Track_14$Trackway_width
#> [1] 0.1346045
#> 
#> $Track_14$Pace_angulation
#> [1] 164.6249
#> 
#> 
#> $Track_15
#> $Track_15$Turning_angles
#> [1] 114.2659 111.3107 118.9163 120.6646
#> 
#> $Track_15$Mean_turning_angle
#> [1] 116.2898
#> 
#> $Track_15$Standard_deviation_turning_angle
#> [1] 4.278817
#> 
#> $Track_15$Distance
#> [1] 4.556698
#> 
#> $Track_15$Length
#> [1] 4.56646
#> 
#> $Track_15$Step_lengths
#> [1] 1.057393 1.195731 1.154191 1.159145
#> 
#> $Track_15$Mean_step_length
#> [1] 1.141615
#> 
#> $Track_15$Standard_deviation_step_length
#> [1] 0.05912482
#> 
#> $Track_15$Sinuosity
#> [1] 0.07871941
#> 
#> $Track_15$Straightness
#> [1] 0.9978623
#> 
#> $Track_15$Trackway_width
#> [1] 0.0005184741
#> 
#> $Track_15$Pace_angulation
#> [1] 171.2543
#> 
#> 
#> $Track_16
#> $Track_16$Turning_angles
#> [1] 116.5651 122.1172 120.0741 121.6150 129.5226 129.5917 135.7970
#> 
#> $Track_16$Mean_turning_angle
#> [1] 125.0358
#> 
#> $Track_16$Standard_deviation_turning_angle
#> [1] 6.748851
#> 
#> $Track_16$Distance
#> [1] 7.068614
#> 
#> $Track_16$Length
#> [1] 7.109877
#> 
#> $Track_16$Step_lengths
#> [1] 0.9075806 0.9584426 1.0476549 1.0982925 1.0294856 1.0591361 1.0092847
#> 
#> $Track_16$Mean_step_length
#> [1] 1.015697
#> 
#> $Track_16$Standard_deviation_step_length
#> [1] 0.06445766
#> 
#> $Track_16$Sinuosity
#> [1] 0.08320294
#> 
#> $Track_16$Straightness
#> [1] 0.9941965
#> 
#> $Track_16$Trackway_width
#> [1] 0.1388154
#> 
#> $Track_16$Pace_angulation
#> [1] 170.4862
#> 
#> 
#> $Track_17
#> $Track_17$Turning_angles
#> [1] -131.2245
#> 
#> $Track_17$Mean_turning_angle
#> [1] 228.7755
#> 
#> $Track_17$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_17$Distance
#> [1] 0.3079478
#> 
#> $Track_17$Length
#> [1] 0.3079478
#> 
#> $Track_17$Step_lengths
#> [1] 0.3079478
#> 
#> $Track_17$Mean_step_length
#> [1] 0.3079478
#> 
#> $Track_17$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_17$Sinuosity
#> [1] NaN
#> 
#> $Track_17$Straightness
#> [1] 1
#> 
#> $Track_17$Trackway_width
#> [1] 0.03340355
#> 
#> $Track_17$Pace_angulation
#> [1] 167.6087
#> 
#> 
#> $Track_18
#> $Track_18$Turning_angles
#> [1] 107.7676 110.0561 113.7495 110.8756
#> 
#> $Track_18$Mean_turning_angle
#> [1] 110.6121
#> 
#> $Track_18$Standard_deviation_turning_angle
#> [1] 2.470679
#> 
#> $Track_18$Distance
#> [1] 4.767794
#> 
#> $Track_18$Length
#> [1] 4.771307
#> 
#> $Track_18$Step_lengths
#> [1] 1.257788 1.183538 1.265235 1.064747
#> 
#> $Track_18$Mean_step_length
#> [1] 1.192827
#> 
#> $Track_18$Standard_deviation_step_length
#> [1] 0.09301194
#> 
#> $Track_18$Sinuosity
#> [1] 0.04807342
#> 
#> $Track_18$Straightness
#> [1] 0.9992637
#> 
#> $Track_18$Trackway_width
#> [1] 0.1359546
#> 
#> $Track_18$Pace_angulation
#> [1] 164.4785
#> 
#> 
#> $Track_19
#> $Track_19$Turning_angles
#> [1] 111.8959
#> 
#> $Track_19$Mean_turning_angle
#> [1] 111.8959
#> 
#> $Track_19$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_19$Distance
#> [1] 1.490617
#> 
#> $Track_19$Length
#> [1] 1.490617
#> 
#> $Track_19$Step_lengths
#> [1] 1.490617
#> 
#> $Track_19$Mean_step_length
#> [1] 1.490617
#> 
#> $Track_19$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_19$Sinuosity
#> [1] NaN
#> 
#> $Track_19$Straightness
#> [1] 1
#> 
#> $Track_19$Trackway_width
#> [1] 0.00287917
#> 
#> $Track_19$Pace_angulation
#> [1] 179.7727
#> 
#> 
#> $Track_20
#> $Track_20$Turning_angles
#> [1] 104.3162 108.0343
#> 
#> $Track_20$Mean_turning_angle
#> [1] 106.1753
#> 
#> $Track_20$Standard_deviation_turning_angle
#> [1] 2.629063
#> 
#> $Track_20$Distance
#> [1] 2.300393
#> 
#> $Track_20$Length
#> [1] 2.301602
#> 
#> $Track_20$Step_lengths
#> [1] 1.204311 1.097291
#> 
#> $Track_20$Mean_step_length
#> [1] 1.150801
#> 
#> $Track_20$Standard_deviation_step_length
#> [1] 0.07567463
#> 
#> $Track_20$Sinuosity
#> [1] 0.06051244
#> 
#> $Track_20$Straightness
#> [1] 0.9994748
#> 
#> $Track_20$Trackway_width
#> [1] 0.005588892
#> 
#> $Track_20$Pace_angulation
#> [1] 176.2159
#> 
#> 
#> $Track_21
#> $Track_21$Turning_angles
#> [1] 140.7497 134.7000
#> 
#> $Track_21$Mean_turning_angle
#> [1] 137.7249
#> 
#> $Track_21$Standard_deviation_turning_angle
#> [1] 4.277776
#> 
#> $Track_21$Distance
#> [1] 2.325026
#> 
#> $Track_21$Length
#> [1] 2.328268
#> 
#> $Track_21$Step_lengths
#> [1] 1.136569 1.191699
#> 
#> $Track_21$Mean_step_length
#> [1] 1.164134
#> 
#> $Track_21$Standard_deviation_step_length
#> [1] 0.03898298
#> 
#> $Track_21$Sinuosity
#> [1] 0.09795169
#> 
#> $Track_21$Straightness
#> [1] 0.9986075
#> 
#> $Track_21$Trackway_width
#> [1] 0.01064242
#> 
#> $Track_21$Pace_angulation
#> [1] 173.1296
#> 
#> 
#> $Track_22
#> $Track_22$Turning_angles
#> [1] 118.3636 109.6538
#> 
#> $Track_22$Mean_turning_angle
#> [1] 114.0087
#> 
#> $Track_22$Standard_deviation_turning_angle
#> [1] 6.158744
#> 
#> $Track_22$Distance
#> [1] 1.795854
#> 
#> $Track_22$Length
#> [1] 1.801009
#> 
#> $Track_22$Step_lengths
#> [1] 0.8172251 0.9837840
#> 
#> $Track_22$Mean_step_length
#> [1] 0.9005046
#> 
#> $Track_22$Standard_deviation_step_length
#> [1] 0.1177749
#> 
#> $Track_22$Sinuosity
#> [1] 0.1604935
#> 
#> $Track_22$Straightness
#> [1] 0.9971376
#> 
#> $Track_22$Trackway_width
#> [1] 0.2363052
#> 
#> $Track_22$Pace_angulation
#> [1] 149.1711
#> 
#> 
#> $Track_23
#> $Track_23$Turning_angles
#> [1] -9.811257
#> 
#> $Track_23$Mean_turning_angle
#> [1] 350.1887
#> 
#> $Track_23$Standard_deviation_turning_angle
#> [1] NA
#> 
#> $Track_23$Distance
#> [1] 1.786422
#> 
#> $Track_23$Length
#> [1] 1.786422
#> 
#> $Track_23$Step_lengths
#> [1] 1.786422
#> 
#> $Track_23$Mean_step_length
#> [1] 1.786422
#> 
#> $Track_23$Standard_deviation_step_length
#> [1] NA
#> 
#> $Track_23$Sinuosity
#> [1] NaN
#> 
#> $Track_23$Straightness
#> [1] 1
#> 
#> $Track_23$Trackway_width
#> [1] 0.1879773
#> 
#> $Track_23$Pace_angulation
#> [1] 167.7341
#> 
#> 
```
