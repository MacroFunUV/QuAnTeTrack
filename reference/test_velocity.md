# Test for differences in velocity means with pairwise comparisons

`test_velocity()` evaluates differences in velocity means across
different tracks using a specified statistical test. It includes options
for ANOVA, Kruskal-Wallis test, and Generalized Linear Models (GLM), and
checks for assumptions such as normality and homogeneity of variances.
For datasets with more than two tracks, it performs pairwise comparisons
to identify specific differences between tracks.

## Usage

``` r
test_velocity(data, trackvel, plot = FALSE, analysis = NULL)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- trackvel:

  A `track velocity` R object consisting of a list where each element
  corresponds to a track and contains velocity or relative stride length
  data.

- plot:

  A logical value indicating whether to plot a boxplot of velocities by
  track (default is `FALSE`).

- analysis:

  A character string specifying the type of analysis: `"ANOVA"`,
  `"Kruskal-Wallis"`, or `"GLM"`. Default is `"ANOVA"`.

## Value

A list with the results of the statistical analysis and diagnostic
tests:

- `normality_results`: A matrix of test statistics and *p*-values from
  the Shapiro-Wilk test for each track, with rows for the test statistic
  and *p*-value, and columns for each track.

- `homogeneity_test`: The result of Levene's test, including the
  *p*-value for homogeneity of variances.

- `ANOVA` (If `analysis` is `"ANOVA"`): A list containing the ANOVA
  table and Tukey HSD post-hoc test results.

- `Kruskal_Wallis` (If `analysis` is `"Kruskal-Wallis"`): A list
  containing the Kruskal-Wallis test result and Dunn's test post-hoc
  results.

- `GLM` (If `analysis` is `"GLM"`): A summary of the GLM fit and
  pairwise comparisons.

- `plot` (If `plot` is `TRUE`): A boxplot of velocities by track is
  generated and displayed.

## Details

The `test_velocity` function performs the following operations:

- **Condition Testing:**

  - **Normality:** Shapiro-Wilk test for normality on velocity data
    within each track.

  - **Homogeneity of Variances:** Levene's test for equal variances
    across tracks.

- **Statistical Analysis:**

  - **ANOVA:** Compares mean velocities across tracks, assuming
    normality and homogeneity of variances. Includes Tukey's HSD
    post-hoc test for pairwise comparisons.

  - **Kruskal-Wallis Test:** Non-parametric alternative to ANOVA for
    comparing median velocities across tracks when assumptions are
    violated. Includes Dunn's test for pairwise comparisons.

  - **GLM:** Generalized Linear Model with a Gaussian family for
    comparing means if ANOVA assumptions are not met. Pairwise
    comparisons in the GLM are conducted using estimated marginal means
    (least-squares means) with the emmeans package, which computes
    differences between group means while adjusting for multiple
    comparisons using Tukey’s method.

- **Plotting:**

  - If `plot` is `TRUE`, a boxplot of velocities by track is generated.

## Logo

![](figures/Logo.png)

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`velocity_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/velocity_track.md)

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
# Example 1: Test for Differences in Velocity Means with Pairwise Comparisons in Trajectories
# in MountTom dataset.

# Hip heights for each track in the MountTom dataset
H_mounttom <- c(
  1.380, 1.404, 1.320, 1.736, 1.364, 1.432, 1.508, 1.768, 1.600, 1.848,
  1.532, 1.532, 0.760, 1.532, 1.688, 1.620, 0.636, 1.784, 1.676, 1.872,
  1.648, 1.760, 1.612
)

# Calculate velocities using the default Method "A"
V_mounttom <- velocity_track(MountTom, H = H_mounttom)

# Test for Differences in Velocity Means with Pairwise Comparisons
test_velocity(MountTom, V_mounttom)
#> Warning: The following tracks were removed from the analysis due to having 3 or fewer footprints: Track 05, Track 06, Track 10, Track 11, Track 12, Track 14, Track 17, Track 19, Track 20, Track 21, Track 22, Track 23.
#> $normality_results
#>              Track 01  Track 02  Track 03  Track 04  Track 07   Track 08
#> statistic.W 0.9209702 0.9049639 0.8522781 0.9286605 0.9371515 0.81516913
#> p_value     0.4003329 0.2821994 0.2018271 0.5872649 0.6132123 0.08011557
#>              Track 09  Track 13  Track 15  Track 16  Track 18
#> statistic.W 0.8697921 0.9420435 0.8148552 0.9692820 0.8301712
#> p_value     0.2655830 0.6804196 0.1064974 0.8923149 0.1395325
#> 
#> $homogeneity_test
#> Levene's Test for Homogeneity of Variance (center = median)
#>       Df F value  Pr(>F)  
#> group 10  1.8073 0.07947 .
#>       58                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $ANOVA
#> $ANOVA$ANOVA
#>             Df Sum Sq Mean Sq F value   Pr(>F)    
#> track       10 10.729  1.0729      19 7.35e-15 ***
#> Residuals   58  3.275  0.0565                     
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $ANOVA$Tukey
#>   Tukey multiple comparisons of means
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = vel ~ track, data = M_analysis)
#> 
#> $track
#>                           diff           lwr          upr     p adj
#> Track 02-Track 01  0.194005116 -0.1812699957  0.569280228 0.8128384
#> Track 03-Track 01  0.966158204  0.5221267034  1.410189704 0.0000001
#> Track 04-Track 01 -0.058483846 -0.5025153465  0.385547654 0.9999966
#> Track 07-Track 01 -0.337960523 -0.7391464928  0.063225447 0.1756512
#> Track 08-Track 01  0.105284346 -0.3142859841  0.524854676 0.9987899
#> Track 09-Track 01  0.663995566  0.2199640654  1.108027066 0.0002666
#> Track 13-Track 01 -0.703166313 -1.1471978130 -0.259134812 0.0000926
#> Track 15-Track 01  0.050096645 -0.3939348551  0.494128146 0.9999992
#> Track 16-Track 01 -0.189837124 -0.5766618555  0.196987607 0.8563548
#> Track 18-Track 01 -0.000824843 -0.4448563435  0.443206657 1.0000000
#> Track 03-Track 02  0.772153088  0.3281215872  1.216184588 0.0000136
#> Track 04-Track 02 -0.252488962 -0.6965204627  0.191542538 0.7114764
#> Track 07-Track 02 -0.531965639 -0.9331516090 -0.130779669 0.0018707
#> Track 08-Track 02 -0.088720770 -0.5082911003  0.330849560 0.9997273
#> Track 09-Track 02  0.469990450  0.0259589492  0.914021950 0.0294239
#> Track 13-Track 02 -0.897171429 -1.3412029292 -0.453139928 0.0000004
#> Track 15-Track 02 -0.143908471 -0.5879399713  0.300123030 0.9904498
#> Track 16-Track 02 -0.383842241 -0.7706669717  0.002982491 0.0534952
#> Track 18-Track 02 -0.194829959 -0.6388614597  0.249201541 0.9236910
#> Track 04-Track 03 -1.024642050 -1.5281264461 -0.521157654 0.0000003
#> Track 07-Track 03 -1.304118727 -1.7702547008 -0.837982753 0.0000000
#> Track 08-Track 03 -0.860873858 -1.3429234671 -0.378824249 0.0000076
#> Track 09-Track 03 -0.302162638 -0.8056470343  0.201321758 0.6425087
#> Track 13-Track 03 -1.669324516 -2.1728089127 -1.165840120 0.0000000
#> Track 15-Track 03 -0.916061559 -1.4195459547 -0.412577162 0.0000050
#> Track 16-Track 03 -1.155995328 -1.6098300300 -0.702160626 0.0000000
#> Track 18-Track 03 -0.966983047 -1.4704674431 -0.463498651 0.0000014
#> Track 07-Track 04 -0.279476677 -0.7456126509  0.186659297 0.6438142
#> Track 08-Track 04  0.163768192 -0.3182814172  0.645817801 0.9864022
#> Track 09-Track 04  0.722479412  0.2189950157  1.225963808 0.0005416
#> Track 13-Track 04 -0.644682467 -1.1481668627 -0.141198070 0.0030869
#> Track 15-Track 04  0.108580491 -0.3949039048  0.612064888 0.9996753
#> Track 16-Track 04 -0.131353278 -0.5851879800  0.322481423 0.9960595
#> Track 18-Track 04  0.057659003 -0.4458253932  0.561143399 0.9999991
#> Track 08-Track 07  0.443244869  0.0003475505  0.886142188 0.0496557
#> Track 09-Track 07  1.001956089  0.5358201149  1.468092063 0.0000001
#> Track 13-Track 07 -0.365205790 -0.8313417635  0.100930184 0.2604157
#> Track 15-Track 07  0.388057168 -0.0780788056  0.854193142 0.1882672
#> Track 16-Track 07  0.148123399 -0.2638864865  0.560133284 0.9794893
#> Track 18-Track 07  0.337135680 -0.1290002940  0.803271654 0.3706410
#> Track 09-Track 08  0.558711220  0.0766616105  1.040760829 0.0110562
#> Track 13-Track 08 -0.808450659 -1.2905002679 -0.326401049 0.0000294
#> Track 15-Track 08 -0.055187701 -0.5372373100  0.426861909 0.9999991
#> Track 16-Track 08 -0.295121470 -0.7250531220  0.134810181 0.4479344
#> Track 18-Track 08 -0.106109189 -0.5881587984  0.375940420 0.9996108
#> Track 13-Track 09 -1.367161878 -1.8706462746 -0.863677482 0.0000000
#> Track 15-Track 09 -0.613898920 -1.1173833167 -0.110414524 0.0059373
#> Track 16-Track 09 -0.853832690 -1.3076673919 -0.399997988 0.0000023
#> Track 18-Track 09 -0.664820409 -1.1683048051 -0.161336013 0.0019895
#> Track 15-Track 13  0.753262958  0.2497785617  1.256747354 0.0002643
#> Track 16-Track 13  0.513329188  0.0594944865  0.967163890 0.0146246
#> Track 18-Track 13  0.702341470  0.1988570733  1.205825866 0.0008589
#> Track 16-Track 15 -0.239933770 -0.6937684714  0.213900932 0.7915111
#> Track 18-Track 15 -0.050921488 -0.5544058846  0.452562908 0.9999997
#> Track 18-Track 16  0.189012281 -0.2648224205  0.642846983 0.9447031
#> 
#> 
#> 

# Example 2: Test for Differences in Velocity Means with Pairwise Comparisons in Trajectories
# in PaluxyRiver dataset.

# Hip heights for each track in the PaluxyRiver dataset
H_paluxyriver <- c(3.472, 2.200)

# Specify different methods for different tracks
Method_paluxyriver <- c("A", "B")

# Calculate velocities using specified methods
V_paluxyriver <- velocity_track(PaluxyRiver, H = H_paluxyriver, method = Method_paluxyriver)

# Test for Differences in Velocity Means with Pairwise Comparisons
test_velocity(PaluxyRiver, V_paluxyriver)
#> Warning: One or more tracks do not follow a normal distribution (p-value <= 0.05). Assumptions for ANOVA are not met. Consider using 'Kruskal-Wallis' or 'GLM'.
#> Warning: Homogeneity of variances assumption is violated (Levene's test p-value <= 0.05). Assumptions for ANOVA are not met. Consider using 'Kruskal-Wallis' or 'GLM'.
#> $normality_results
#>               Track 1    Track 2
#> statistic.W 0.9586486 0.91134871
#> p_value     0.3043396 0.03773134
#> 
#> $homogeneity_test
#> Levene's Test for Homogeneity of Variance (center = median)
#>       Df F value    Pr(>F)    
#> group  1  18.698 7.121e-05 ***
#>       51                      
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $ANOVA
#> $ANOVA$ANOVA
#>             Df Sum Sq Mean Sq F value   Pr(>F)    
#> track        1 0.6687  0.6687   111.9 1.82e-14 ***
#> Residuals   51 0.3046  0.0060                     
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> $ANOVA$Tukey
#>   Tukey multiple comparisons of means
#>     95% family-wise confidence level
#> 
#> Fit: aov(formula = vel ~ track, data = M_analysis)
#> 
#> $track
#>                      diff       lwr     upr p adj
#> Track 2-Track 1 0.2256526 0.1828352 0.26847     0
#> 
#> 
#> 
```
