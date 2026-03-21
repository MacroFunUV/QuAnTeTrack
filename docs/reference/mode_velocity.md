# Test for steady, acceleration, or deceleration along trajectories

`mode_velocity()` evaluates the trend in velocity along each trajectory
by applying Spearman's rank correlation test. The function classifies
the trend into "acceleration", "deceleration", or "steady" based on the
correlation and the *p*-value.

## Usage

``` r
mode_velocity(trackvel)
```

## Arguments

- trackvel:

  A `track velocity` object where each element corresponds to a track
  and contains a vector of velocity or relative stride length data.

## Value

A list where each element corresponds to a trajectory from the input
`trackvel` and contains:

- **correlation:** The result of the Spearman correlation test,
  including the correlation coefficient and *p*-value.

- **trend:** A classification of the trend as "Acceleration",
  "Deceleration", or "Steady" based on the *p*-value and the correlation
  coefficient.

- If a trajectory has fewer than 3 steps, the entry contains the message
  "Less than three steps."

## Details

The `mode_velocity()` function performs the following operations:

- **Spearman's Rank Correlation Test:**

  - This non-parametric test assesses the strength and direction of a
    monotonic relationship between two variables. It does not require
    assumptions about the normality of data or a linear relationship
    between velocity and step number.

  - It uses ranks rather than raw values, making it robust to outliers
    and suitable for detecting general trends (acceleration or
    deceleration) in velocity data.

- **Function Operation:**

  - For each trajectory in the `trackvel` list, the function calculates
    the Spearman correlation coefficient and the associated *p*-value
    between velocity and step number.

  - Based on the *p*-value and correlation coefficient, it classifies
    the trend as "acceleration", "deceleration", or "steady".

  - If a trajectory contains fewer than 3 steps, the function returns a
    message indicating insufficient data for correlation analysis.

- **Advantages:**

  - The non-parametric nature allows flexibility with data distributions
    and reduced sensitivity to outliers compared to parametric tests.

  - Effective for detecting monotonic trends (either increasing or
    decreasing) when the correlation is statistically significant.

- **Limitations:**

  - May be unreliable with very small sample sizes (e.g., fewer than 3
    steps), providing potentially non-informative results.

  - Does not capture the magnitude of change or provide detailed
    insights into the rate of acceleration or deceleration.

  - Identifies monotonic trends based on statistical significance but
    does not distinguish between different types of monotonic
    relationships (e.g., steady acceleration vs. abrupt changes).

**Interpretation of Results:**

- **Acceleration:** If the *p*-value is less than 0.05 and the Spearman
  correlation coefficient is positive.

- **Deceleration:** If the *p*-value is less than 0.05 and the Spearman
  correlation coefficient is negative.

- **Steady:** If the *p*-value is greater than or equal to 0.05,
  indicating no significant monotonic relationship.

**Usage Considerations:**

- Ensure that each trajectory in `trackvel` has a sufficient number of
  steps for meaningful analysis.

- For more detailed analysis of velocity trends, consider complementary
  methods such as linear or non-linear regression, or specialized change
  point detection techniques.

## Logo

![](figures/Logo.png)

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`velocity_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/velocity_track.md),
[`plot_velocity`](https://macrofunuv.github.io/QuAnTeTrack/reference/plot_velocity.md)

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
# Example 1: Test for Steady, Acceleration, or Deceleration in MountTom dataset.

# Hip heights for each track in the MountTom dataset
H_mounttom <- c(
  1.380, 1.404, 1.320, 1.736, 1.364, 1.432, 1.508, 1.768, 1.600,
  1.848, 1.532, 1.532, 0.760, 1.532, 1.688, 1.620, 0.636, 1.784,
  1.676, 1.872, 1.648, 1.760, 1.612
)

# Calculate velocities using the default Method "A"
V_mounttom <- velocity_track(MountTom, H = H_mounttom)

# Test for Steady, Acceleration, or Deceleration
mode_velocity(V_mounttom)
#> $Track_01
#> $Track_01$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 160, p-value = 0.002008
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>        rho 
#> -0.9047619 
#> 
#> 
#> $Track_01$trend
#> [1] "Deceleration"
#> 
#> 
#> $Track_02
#> $Track_02$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 42, p-value = 0.207
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#> rho 
#> 0.5 
#> 
#> 
#> $Track_02$trend
#> [1] "Steady"
#> 
#> 
#> $Track_03
#> $Track_03$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 12, p-value = 0.8
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>  rho 
#> -0.2 
#> 
#> 
#> $Track_03$trend
#> [1] "Steady"
#> 
#> 
#> $Track_04
#> $Track_04$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 18, p-value = 0.2
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>  rho 
#> -0.8 
#> 
#> 
#> $Track_04$trend
#> [1] "Steady"
#> 
#> 
#> $Track_05
#> [1] "Less than three steps"
#> 
#> $Track_06
#> [1] "Less than three steps"
#> 
#> $Track_07
#> $Track_07$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 12, p-value = 0.1562
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>       rho 
#> 0.6571429 
#> 
#> 
#> $Track_07$trend
#> [1] "Steady"
#> 
#> 
#> $Track_08
#> $Track_08$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 20, p-value = 1
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#> rho 
#>   0 
#> 
#> 
#> $Track_08$trend
#> [1] "Steady"
#> 
#> 
#> $Track_09
#> $Track_09$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 20, p-value < 2.2e-16
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#> rho 
#>  -1 
#> 
#> 
#> $Track_09$trend
#> [1] "Deceleration"
#> 
#> 
#> $Track_10
#> [1] "Less than three steps"
#> 
#> $Track_11
#> [1] "Less than three steps"
#> 
#> $Track_12
#> [1] "Less than three steps"
#> 
#> $Track_13
#> $Track_13$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 14, p-value = 0.6
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>  rho 
#> -0.4 
#> 
#> 
#> $Track_13$trend
#> [1] "Steady"
#> 
#> 
#> $Track_14
#> [1] "Less than three steps"
#> 
#> $Track_15
#> $Track_15$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 6, p-value = 0.6
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#> rho 
#> 0.4 
#> 
#> 
#> $Track_15$trend
#> [1] "Steady"
#> 
#> 
#> $Track_16
#> $Track_16$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 30, p-value = 0.2939
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>       rho 
#> 0.4642857 
#> 
#> 
#> $Track_16$trend
#> [1] "Steady"
#> 
#> 
#> $Track_17
#> [1] "Less than three steps"
#> 
#> $Track_18
#> $Track_18$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 14, p-value = 0.6
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>  rho 
#> -0.4 
#> 
#> 
#> $Track_18$trend
#> [1] "Steady"
#> 
#> 
#> $Track_19
#> [1] "Less than three steps"
#> 
#> $Track_20
#> [1] "Less than three steps"
#> 
#> $Track_21
#> [1] "Less than three steps"
#> 
#> $Track_22
#> [1] "Less than three steps"
#> 
#> $Track_23
#> [1] "Less than three steps"
#> 

# Example 2: Test for Steady, Acceleration, or Deceleration in PaluxyRiver dataset.

# Hip heights for each track in the PaluxyRiver dataset
H_paluxyriver <- c(3.472, 2.200)

# Specify different methods for different tracks
Method_paluxyriver <- c("A", "B")

# Calculate velocities using specified methods
V_paluxyriver <- velocity_track(PaluxyRiver,
  H = H_paluxyriver,
  method = Method_paluxyriver
)

# Test for Steady, Acceleration, or Deceleration
mode_velocity(V_paluxyriver)
#> $Track_1
#> $Track_1$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 3178, p-value = 0.5088
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>       rho 
#> 0.1302682 
#> 
#> 
#> $Track_1$trend
#> [1] "Steady"
#> 
#> 
#> $Track_2
#> $Track_2$correlation
#> 
#>  Spearman's rank correlation rho
#> 
#> data:  velocity and steps
#> S = 1426, p-value = 0.1711
#> alternative hypothesis: true rho is not equal to 0
#> sample estimates:
#>       rho 
#> 0.2954545 
#> 
#> 
#> $Track_2$trend
#> [1] "Steady"
#> 
#> 
```
