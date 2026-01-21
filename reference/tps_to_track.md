# Transform a \*.tps file into a `trackway` R object

`tps_to_track()` reads a \*.tps file containing footprint coordinates of
one or several trackways and transforms it into a `trackway` R object.

## Usage

``` r
tps_to_track(file, scale = NULL, R.L.side, missing = FALSE, NAs = NULL)
```

## Arguments

- file:

  A \*.tps file containing (x,y) coordinates of footprints in trackways.

- scale:

  A numeric value specifying the scale in meters per pixel.

- R.L.side:

  A character vector specifying the laterality of the first footprint of
  each trackway. The length of the vector must be equal to the total
  number of trackways in the sample.

  - `"L"`: first footprint corresponds to the left foot.

  - `"R"`: first footprint corresponds to the right foot.

- missing:

  A logical value indicating whether there are missing footprints in any
  trackway to be interpolated: `TRUE`, or `FALSE` (the default).

- NAs:

  A matrix with two columns indicating which missing footprints will be
  interpolated. The first column gives the number of the trackway
  containing missing footprints, and the second column gives the number
  of the footprint that is missing within this trackway. The number of
  rows is equal to the total number of missing footprints in the sample.

## Value

A `trackway` R object, which is a list consisting of two elements:

- **`Trajectories`**: A list of trajectories representing trackway
  midlines, interpolated by connecting the midpoints of successive
  left–right footprint pairs (i.e., footprints linked by pace lines).
  Includes columns `X`, `Y`, `IMAGE`, `ID`, and `Side` (set to
  `"Medial"`).

- **`Footprints`**: A list of data frames containing footprint
  coordinates and associated metadata, with a `Side` column (`"R"` or
  `"L"`) and a `missing` marker (`"Actual"` or `"Inferred"`).

## Details

It is highly recommended that the \*.tps file is built using the TPS
software (Rohlf 2008, 2009). Trackways with a different number of
footprints (i.e., landmarks) are allowed. This function transforms the
coordinates of the footprints of each trackway into a set of trajectory
coordinates. Each point of the trajectory is calculated as:
\$\$Point_i(x,y)= (Footprint_i(x,y) + Footprint\_{i+1}(x,y)/2\$\$

The number of points of the resulting trajectory is \\n\_{footprints} -
1\\.

If `missing` is set to `TRUE`, missing footprints can be interpolated.
The interpolated footprint is then placed at the midpoint between the
two adjacent footprints and shifted laterally along the direction
perpendicular to their connecting segment. The magnitude of this lateral
shift is estimated from the mean trackway width of the corresponding
trackway, and its sign is determined by the inferred anatomical side
(left or right) of the missing footprint. Inferred footprints are
flagged as such in the output.

## Logo

![](figures/Logo.png)

## References

Farlow, J. O., O’Brien, M., Kuban, G. J., Dattilo, B. F., Bates, K. T.,
Falkingham, P. L., & Piñuela, L. (2012). Dinosaur Tracksites of the
Paluxy River Valley (Glen Rose Formation, Lower Cretaceous), Dinosaur
Valley State Park, Somervell County, Texas. In Proceedings of the V
International Symposium about Dinosaur Palaeontology and their
Environment (pp. 41-69). Burgos: Salas de los Infantes.

Ostrom, J. H. (1972). Were some dinosaurs gregarious?. Palaeogeography,
Palaeoclimatology, Palaeoecology, 11(4), 287-301.

Rohlf, F. J. 2008. TPSUTIL. Version 1.40. Department of Ecology and
Evolution, State University of New York. Available at
<https://sbmorphometrics.org/>.

Rohlf, F. J. 2009. tpsDig. Version 2.14. Department of Ecology and
Evolution, State University of New York. Available at
<https://sbmorphometrics.org/>.

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
# Example 1: Trackways without missing footprints.
# Based on the Paluxy River dinosaur chase sequence (Farlow et al., 2011).

# Load the example TPS file provided in the QuAnTeTrack package.
tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")

# Convert the TPS data into a trackway object.
# The 'scale' argument sets the scaling factor,
# 'R.L.side' specifies the starting side for each trackway,
# and 'missing = FALSE' indicates no footprints are missing.
tps_to_track(
  tpsPaluxyRiver,
  scale = 0.004341493,
  R.L.side = c("R", "L"),
  missing = FALSE,
  NAs = NULL
)
#> $Trajectories
#> $Trajectories$Trackway_1
#>            x         y        IMAGE ID   Side time displacementTime
#> 1  0.7554198  1.026763 Sauropod.png  0 Medial 0.00             0.00
#> 2  0.8313959  1.680158 Sauropod.png  0 Medial 0.02             0.02
#> 3  0.9703237  2.292308 Sauropod.png  0 Medial 0.04             0.04
#> 4  1.0680073  2.850190 Sauropod.png  0 Medial 0.06             0.06
#> 5  1.1656909  3.375511 Sauropod.png  0 Medial 0.08             0.08
#> 6  1.3176431  3.929051 Sauropod.png  0 Medial 0.10             0.10
#> 7  1.4934736  4.532519 Sauropod.png  0 Medial 0.12             0.12
#> 8  1.6606211  5.112108 Sauropod.png  0 Medial 0.14             0.14
#> 9  1.8234271  5.663478 Sauropod.png  0 Medial 0.16             0.16
#> 10 1.9970868  6.219189 Sauropod.png  0 Medial 0.18             0.18
#> 11 2.0795751  6.655509 Sauropod.png  0 Medial 0.20             0.20
#> 12 2.1012826  7.120049 Sauropod.png  0 Medial 0.22             0.22
#> 13 2.1642343  7.701809 Sauropod.png  0 Medial 0.24             0.24
#> 14 2.2402104  8.290081 Sauropod.png  0 Medial 0.26             0.26
#> 15 2.3009913  8.854475 Sauropod.png  0 Medial 0.28             0.28
#> 16 2.3747967  9.384137 Sauropod.png  0 Medial 0.30             0.30
#> 17 2.4898462  9.933336 Sauropod.png  0 Medial 0.32             0.32
#> 18 2.5983836 10.543316 Sauropod.png  0 Medial 0.34             0.34
#> 19 2.6548230 11.146783 Sauropod.png  0 Medial 0.36             0.36
#> 20 2.6743597 11.824056 Sauropod.png  0 Medial 0.38             0.38
#> 21 2.6396277 12.429694 Sauropod.png  0 Medial 0.40             0.40
#> 22 2.5636516 12.955015 Sauropod.png  0 Medial 0.42             0.42
#> 23 2.4985292 13.417384 Sauropod.png  0 Medial 0.44             0.44
#> 24 2.4268946 13.949217 Sauropod.png  0 Medial 0.46             0.46
#> 25 2.3530892 14.589587 Sauropod.png  0 Medial 0.48             0.48
#> 26 2.3053328 15.195226 Sauropod.png  0 Medial 0.50             0.50
#> 27 2.2488934 15.790010 Sauropod.png  0 Medial 0.52             0.52
#> 28 2.1859417 16.386965 Sauropod.png  0 Medial 0.54             0.54
#> 29 2.1121363 17.059897 Sauropod.png  0 Medial 0.56             0.56
#>                   polar           displacement
#> 1  0.7554198+ 1.026763i  0.00000000+0.0000000i
#> 2  0.8313959+ 1.680158i  0.07597613+0.6533947i
#> 3  0.9703237+ 2.292308i  0.13892778+0.6121505i
#> 4  1.0680073+ 2.850190i  0.09768359+0.5578819i
#> 5  1.1656909+ 3.375511i  0.09768359+0.5253207i
#> 6  1.3176431+ 3.929051i  0.15195226+0.5535404i
#> 7  1.4934736+ 4.532519i  0.17583047+0.6034675i
#> 8  1.6606211+ 5.112108i  0.16714748+0.5795893i
#> 9  1.8234271+ 5.663478i  0.16280599+0.5513696i
#> 10 1.9970868+ 6.219189i  0.17365972+0.5557111i
#> 11 2.0795751+ 6.655509i  0.08248837+0.4363200i
#> 12 2.1012826+ 7.120049i  0.02170746+0.4645398i
#> 13 2.1642343+ 7.701809i  0.06295165+0.5817601i
#> 14 2.2402104+ 8.290081i  0.07597613+0.5882723i
#> 15 2.3009913+ 8.854475i  0.06078090+0.5643941i
#> 16 2.3747967+ 9.384137i  0.07380538+0.5296621i
#> 17 2.4898462+ 9.933336i  0.11504956+0.5491989i
#> 18 2.5983836+10.543316i  0.10853732+0.6099798i
#> 19 2.6548230+11.146783i  0.05643941+0.6034675i
#> 20 2.6743597+11.824056i  0.01953672+0.6772729i
#> 21 2.6396277+12.429694i -0.03473194+0.6056383i
#> 22 2.5636516+12.955015i -0.07597613+0.5253207i
#> 23 2.4985292+13.417384i -0.06512239+0.4623690i
#> 24 2.4268946+13.949217i -0.07163463+0.5318329i
#> 25 2.3530892+14.589587i -0.07380538+0.6403702i
#> 26 2.3053328+15.195226i -0.04775642+0.6056383i
#> 27 2.2488934+15.790010i -0.05643941+0.5947845i
#> 28 2.1859417+16.386965i -0.06295165+0.5969553i
#> 29 2.1121363+17.059897i -0.07380538+0.6729314i
#> 
#> $Trajectories$Trackway_2
#>            x         y        IMAGE ID   Side time displacementTime
#> 1  0.3646854  1.693182 Theropod.png  1 Medial 0.00             0.00
#> 2  0.4927595  2.279284 Theropod.png  1 Medial 0.02             0.02
#> 3  0.5730771  2.891434 Theropod.png  1 Medial 0.04             0.04
#> 4  0.7206878  3.564366 Theropod.png  1 Medial 0.06             0.06
#> 5  0.9095428  4.213419 Theropod.png  1 Medial 0.08             0.08
#> 6  1.0897147  4.884180 Theropod.png  1 Medial 0.10             0.10
#> 7  1.2959357  5.546257 Theropod.png  1 Medial 0.12             0.12
#> 8  1.4652539  6.171432 Theropod.png  1 Medial 0.14             0.14
#> 9  1.5564252  6.742339 Theropod.png  1 Medial 0.16             0.16
#> 10 1.6150354  7.304562 Theropod.png  1 Medial 0.18             0.18
#> 11 1.6736456  7.899347 Theropod.png  1 Medial 0.20             0.20
#> 12 1.7452802  8.494131 Theropod.png  1 Medial 0.22             0.22
#> 13 1.8690127  9.043330 Theropod.png  1 Medial 0.24             0.24
#> 14 2.0752337  9.724944 Theropod.png  1 Medial 0.26             0.26
#> 15 2.1208193 10.591072 Theropod.png  1 Medial 0.28             0.28
#> 16 2.2597471 11.277028 Theropod.png  1 Medial 0.30             0.30
#> 17 2.2467226 11.832739 Theropod.png  1 Medial 0.32             0.32
#> 18 2.1729172 12.405816 Theropod.png  1 Medial 0.34             0.34
#> 19 2.1186486 12.994089 Theropod.png  1 Medial 0.36             0.36
#> 20 2.0969411 13.667020 Theropod.png  1 Medial 0.38             0.38
#> 21 2.1989662 14.368171 Theropod.png  1 Medial 0.40             0.40
#> 22 2.3292110 15.112737 Theropod.png  1 Medial 0.42             0.42
#> 23 2.3986749 15.902889 Theropod.png  1 Medial 0.44             0.44
#> 24 2.2445519 16.708236 Theropod.png  1 Medial 0.46             0.46
#>                   polar           displacement
#> 1  0.3646854+ 1.693182i  0.00000000+0.0000000i
#> 2  0.4927595+ 2.279284i  0.12807404+0.5861016i
#> 3  0.5730771+ 2.891434i  0.08031762+0.6121505i
#> 4  0.7206878+ 3.564366i  0.14761076+0.6729314i
#> 5  0.9095428+ 4.213419i  0.18885495+0.6490532i
#> 6  1.0897147+ 4.884180i  0.18017196+0.6707607i
#> 7  1.2959357+ 5.546257i  0.20622092+0.6620777i
#> 8  1.4652539+ 6.171432i  0.16931823+0.6251750i
#> 9  1.5564252+ 6.742339i  0.09117135+0.5709063i
#> 10 1.6150354+ 7.304562i  0.05861016+0.5622233i
#> 11 1.6736456+ 7.899347i  0.05861016+0.5947845i
#> 12 1.7452802+ 8.494131i  0.07163463+0.5947845i
#> 13 1.8690127+ 9.043330i  0.12373255+0.5491989i
#> 14 2.0752337+ 9.724944i  0.20622092+0.6816144i
#> 15 2.1208193+10.591072i  0.04558568+0.8661279i
#> 16 2.2597471+11.277028i  0.13892778+0.6859559i
#> 17 2.2467226+11.832739i -0.01302448+0.5557111i
#> 18 2.1729172+12.405816i -0.07380538+0.5730771i
#> 19 2.1186486+12.994089i -0.05426866+0.5882723i
#> 20 2.0969411+13.667020i -0.02170746+0.6729314i
#> 21 2.1989662+14.368171i  0.10202509+0.7011511i
#> 22 2.3292110+15.112737i  0.13024479+0.7445660i
#> 23 2.3986749+15.902889i  0.06946389+0.7901517i
#> 24 2.2445519+16.708236i -0.15412300+0.8053470i
#> 
#> 
#> $Footprints
#> $Footprints[[1]]
#>            X         Y        IMAGE ID Side missing
#> 1  0.9160550  0.633858 Sauropod.png  0    R  Actual
#> 2  0.5947845  1.419668 Sauropod.png  0    L  Actual
#> 3  1.0680073  1.940647 Sauropod.png  0    R  Actual
#> 4  0.8726401  2.643969 Sauropod.png  0    L  Actual
#> 5  1.2633745  3.056411 Sauropod.png  0    R  Actual
#> 6  1.0680073  3.694611 Sauropod.png  0    L  Actual
#> 7  1.5672790  4.163492 Sauropod.png  0    R  Actual
#> 8  1.4196682  4.901546 Sauropod.png  0    L  Actual
#> 9  1.9015739  5.322670 Sauropod.png  0    R  Actual
#> 10 1.7452802  6.004285 Sauropod.png  0    L  Actual
#> 11 2.2488934  6.434093 Sauropod.png  0    R  Actual
#> 12 1.9102569  6.876925 Sauropod.png  0    L  Actual
#> 13 2.2923083  7.363172 Sauropod.png  0    R  Actual
#> 14 2.0361602  8.040445 Sauropod.png  0    L  Actual
#> 15 2.4442606  8.539717 Sauropod.png  0    R  Actual
#> 16 2.1577220  9.169233 Sauropod.png  0    L  Actual
#> 17 2.5918713  9.599041 Sauropod.png  0    R  Actual
#> 18 2.3878212 10.267631 Sauropod.png  0    L  Actual
#> 19 2.8089460 10.819001 Sauropod.png  0    R  Actual
#> 20 2.5007000 11.474566 Sauropod.png  0    L  Actual
#> 21 2.8480194 12.173546 Sauropod.png  0    R  Actual
#> 22 2.4312361 12.685843 Sauropod.png  0    L  Actual
#> 23 2.6960672 13.224188 Sauropod.png  0    R  Actual
#> 24 2.3009913 13.610581 Sauropod.png  0    L  Actual
#> 25 2.5527979 14.287853 Sauropod.png  0    R  Actual
#> 26 2.1533805 14.891321 Sauropod.png  0    L  Actual
#> 27 2.4572850 15.499130 Sauropod.png  0    R  Actual
#> 28 2.0405017 16.080890 Sauropod.png  0    L  Actual
#> 29 2.3313817 16.693041 Sauropod.png  0    R  Actual
#> 30 1.8928909 17.426753 Sauropod.png  0    L  Actual
#> 
#> $Footprints[[2]]
#>            X         Y        IMAGE ID Side missing
#> 1  0.1866842  1.376253 Theropod.png  1    L  Actual
#> 2  0.5426866  2.010111 Theropod.png  1    R  Actual
#> 3  0.4428323  2.548456 Theropod.png  1    L  Actual
#> 4  0.7033219  3.234412 Theropod.png  1    R  Actual
#> 5  0.7380538  3.894319 Theropod.png  1    L  Actual
#> 6  1.0810318  4.532519 Theropod.png  1    R  Actual
#> 7  1.0983977  5.235841 Theropod.png  1    L  Actual
#> 8  1.4934736  5.856674 Theropod.png  1    R  Actual
#> 9  1.4370342  6.486191 Theropod.png  1    L  Actual
#> 10 1.6758163  6.998487 Theropod.png  1    R  Actual
#> 11 1.5542545  7.610637 Theropod.png  1    L  Actual
#> 12 1.7930366  8.188056 Theropod.png  1    R  Actual
#> 13 1.6975238  8.800206 Theropod.png  1    L  Actual
#> 14 2.0405017  9.286454 Theropod.png  1    R  Actual
#> 15 2.1099656 10.163435 Theropod.png  1    L  Actual
#> 16 2.1316731 11.018709 Theropod.png  1    R  Actual
#> 17 2.3878212 11.535347 Theropod.png  1    L  Actual
#> 18 2.1056241 12.130131 Theropod.png  1    R  Actual
#> 19 2.2402104 12.681501 Theropod.png  1    L  Actual
#> 20 1.9970868 13.306676 Theropod.png  1    R  Actual
#> 21 2.1967955 14.027364 Theropod.png  1    L  Actual
#> 22 2.2011370 14.708978 Theropod.png  1    R  Actual
#> 23 2.4572850 15.516496 Theropod.png  1    L  Actual
#> 24 2.3400647 16.289282 Theropod.png  1    R  Actual
#> 25 2.1490390 17.127190 Theropod.png  1    L  Actual
#> 
#> 


# Example 2: Trackways with missing footprints.
# Based on dinosaur trackways from Mount Tom (Ostrom, 1972).

# Load the example TPS file.
tpsMountTom <- system.file("extdata", "MountTom.tps", package = "QuAnTeTrack")

# Define a matrix representing the missing footprints.
# Here, footprint 7 is missing in trackway 3.
NAs <- matrix(c(3, 7), nrow = 1, ncol = 2)

# Convert the TPS data into a trackway object.
# The 'missing = TRUE' flag activates interpolation for missing footprints,
# 'NAs' specifies which footprints are missing,
# and 'R.L.side' indicates the starting side for each trackway.
tps_to_track(
  tpsMountTom,
  scale = 0.004411765,
  R.L.side = c(
    "R", "L", "L", "L", "R", "L", "R", "R", "L", "L", "L",
    "L", "L", "R", "R", "L", "R", "R", "L", "R", "R",
    "R", "R"
  ),
  missing = TRUE,
  NAs = NAs
)
#> $Trajectories
#> $Trajectories$Trackway_01
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 40.67868 15.50294 Track 1.png  0 Medial 0.00             0.00
#> 2 39.91544 16.25515 Track 1.png  0 Medial 0.02             0.02
#> 3 39.21177 16.96103 Track 1.png  0 Medial 0.04             0.04
#> 4 38.45074 17.61838 Track 1.png  0 Medial 0.06             0.06
#> 5 37.75809 18.33088 Track 1.png  0 Medial 0.08             0.08
#> 6 37.14927 19.09632 Track 1.png  0 Medial 0.10             0.10
#> 7 36.43677 19.69853 Track 1.png  0 Medial 0.12             0.12
#> 8 35.74633 20.21030 Track 1.png  0 Medial 0.14             0.14
#> 9 35.03383 20.87427 Track 1.png  0 Medial 0.16             0.16
#>                polar          displacement
#> 1 40.67868+15.50294i  0.0000000+0.0000000i
#> 2 39.91544+16.25515i -0.7632353+0.7522059i
#> 3 39.21177+16.96103i -0.7036765+0.7058824i
#> 4 38.45074+17.61838i -0.7610295+0.6573530i
#> 5 37.75809+18.33088i -0.6926471+0.7125000i
#> 6 37.14927+19.09632i -0.6088236+0.7654412i
#> 7 36.43677+19.69853i -0.7125000+0.6022059i
#> 8 35.74633+20.21030i -0.6904412+0.5117647i
#> 9 35.03383+20.87427i -0.7125000+0.6639706i
#> 
#> $Trajectories$Trackway_02
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 39.18971 14.90735 Track 2.png  1 Medial 0.00             0.00
#> 2 38.56544 15.57794 Track 2.png  1 Medial 0.02             0.02
#> 3 37.86177 16.29044 Track 2.png  1 Medial 0.04             0.04
#> 4 37.12280 17.04485 Track 2.png  1 Medial 0.06             0.06
#> 5 36.46544 17.90074 Track 2.png  1 Medial 0.08             0.08
#> 6 35.83015 18.93750 Track 2.png  1 Medial 0.10             0.10
#> 7 35.24118 19.83750 Track 2.png  1 Medial 0.12             0.12
#> 8 34.58383 20.58530 Track 2.png  1 Medial 0.14             0.14
#> 9 33.95515 21.46324 Track 2.png  1 Medial 0.16             0.16
#>                polar          displacement
#> 1 39.18971+14.90735i  0.0000000+0.0000000i
#> 2 38.56544+15.57794i -0.6242647+0.6705883i
#> 3 37.86177+16.29044i -0.7036765+0.7125000i
#> 4 37.12280+17.04485i -0.7389706+0.7544118i
#> 5 36.46544+17.90074i -0.6573530+0.8558824i
#> 6 35.83015+18.93750i -0.6352942+1.0367648i
#> 7 35.24118+19.83750i -0.5889706+0.9000001i
#> 8 34.58383+20.58530i -0.6573530+0.7477942i
#> 9 33.95515+21.46324i -0.6286765+0.8779412i
#> 
#> $Trajectories$Trackway_03
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 38.10441 15.55147 Track 3.png  2 Medial 0.00             0.00
#> 2 37.27500 16.41618 Track 3.png  2 Medial 0.02             0.02
#> 3 36.36177 17.36912 Track 3.png  2 Medial 0.04             0.04
#> 4 35.45515 18.44118 Track 3.png  2 Medial 0.06             0.06
#> 5 34.68088 19.29927 Track 3.png  2 Medial 0.08             0.08
#>                polar          displacement
#> 1 38.10441+15.55147i  0.0000000+0.0000000i
#> 2 37.27500+16.41618i -0.8294118+0.8647059i
#> 3 36.36177+17.36912i -0.9132354+0.9529412i
#> 4 35.45515+18.44118i -0.9066177+1.0720589i
#> 5 34.68088+19.29927i -0.7742648+0.8580883i
#> 
#> $Trajectories$Trackway_04
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 35.58088 16.01250 Track 4.png  3 Medial 0.00             0.00
#> 2 34.86177 17.00294 Track 4.png  3 Medial 0.02             0.02
#> 3 34.27280 17.97794 Track 4.png  3 Medial 0.04             0.04
#> 4 33.67059 18.83382 Track 4.png  3 Medial 0.06             0.06
#> 5 32.91838 19.63677 Track 4.png  3 Medial 0.08             0.08
#>                polar          displacement
#> 1 35.58088+16.01250i  0.0000000+0.0000000i
#> 2 34.86177+17.00294i -0.7191177+0.9904412i
#> 3 34.27280+17.97794i -0.5889706+0.9750001i
#> 4 33.67059+18.83382i -0.6022059+0.8558824i
#> 5 32.91838+19.63677i -0.7522059+0.8029412i
#> 
#> $Trajectories$Trackway_05
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 32.28750 14.39559 Track 5.png  4 Medial 0.00             0.00
#> 2 31.60368 15.31544 Track 5.png  4 Medial 0.02             0.02
#>                polar         displacement
#> 1 32.28750+14.39559i  0.0000000+0.000000i
#> 2 31.60368+15.31544i -0.6838236+0.919853i
#> 
#> $Trajectories$Trackway_06
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 30.91324 14.48382 Track 6.png  5 Medial 0.00             0.00
#> 2 30.16544 15.26691 Track 6.png  5 Medial 0.02             0.02
#> 3 29.38235 16.07647 Track 6.png  5 Medial 0.04             0.04
#>                polar          displacement
#> 1 30.91324+14.48382i  0.0000000+0.0000000i
#> 2 30.16544+15.26691i -0.7477942+0.7830883i
#> 3 29.38235+16.07647i -0.7830883+0.8095589i
#> 
#> $Trajectories$Trackway_07
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 31.20221 14.97574 Track 7.png  6 Medial 0.00             0.00
#> 2 30.78088 16.09632 Track 7.png  6 Medial 0.02             0.02
#> 3 30.24706 17.14191 Track 7.png  6 Medial 0.04             0.04
#> 4 29.68015 17.99118 Track 7.png  6 Medial 0.06             0.06
#> 5 29.05809 18.81397 Track 7.png  6 Medial 0.08             0.08
#> 6 28.50441 19.56177 Track 7.png  6 Medial 0.10             0.10
#>                polar          displacement
#> 1 31.20221+14.97574i  0.0000000+0.0000000i
#> 2 30.78088+16.09632i -0.4213236+1.1205883i
#> 3 30.24706+17.14191i -0.5338236+1.0455883i
#> 4 29.68015+17.99118i -0.5669118+0.8492648i
#> 5 29.05809+18.81397i -0.6220589+0.8227942i
#> 6 28.50441+19.56177i -0.5536765+0.7477942i
#> 
#> $Trajectories$Trackway_08
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 28.89485 15.10809 Track 8.png  7 Medial 0.00             0.00
#> 2 28.16030 16.11177 Track 8.png  7 Medial 0.02             0.02
#> 3 27.52941 17.17500 Track 8.png  7 Medial 0.04             0.04
#> 4 27.00441 18.19853 Track 8.png  7 Medial 0.06             0.06
#> 5 26.43750 19.10515 Track 8.png  7 Medial 0.08             0.08
#> 6 25.73162 20.15515 Track 8.png  7 Medial 0.10             0.10
#>                polar          displacement
#> 1 28.89485+15.10809i  0.0000000+0.0000000i
#> 2 28.16030+16.11177i -0.7345589+1.0036765i
#> 3 27.52941+17.17500i -0.6308824+1.0632354i
#> 4 27.00441+18.19853i -0.5250000+1.0235295i
#> 5 26.43750+19.10515i -0.5669118+0.9066177i
#> 6 25.73162+20.15515i -0.7058824+1.0500001i
#> 
#> $Trajectories$Trackway_09
#>          x        y       IMAGE ID   Side time displacementTime
#> 1 29.19706 17.09118 Track 9.png  8 Medial 0.00             0.00
#> 2 30.53824 16.41397 Track 9.png  8 Medial 0.02             0.02
#> 3 31.87941 15.90221 Track 9.png  8 Medial 0.04             0.04
#> 4 33.15883 15.41691 Track 9.png  8 Medial 0.06             0.06
#> 5 34.25074 15.04412 Track 9.png  8 Medial 0.08             0.08
#>                polar        displacement
#> 1 29.19706+17.09118i 0.000000+0.0000000i
#> 2 30.53824+16.41397i 1.341177-0.6772059i
#> 3 31.87941+15.90221i 1.341177-0.5117647i
#> 4 33.15883+15.41691i 1.279412-0.4852941i
#> 5 34.25074+15.04412i 1.091912-0.3727941i
#> 
#> $Trajectories$Trackway_10
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 24.15000 16.14485 Track 10.png  9 Medial 0.00             0.00
#> 2 24.73677 15.15662 Track 10.png  9 Medial 0.02             0.02
#>                polar         displacement
#> 1 24.15000+16.14485i 0.0000000+0.0000000i
#> 2 24.73677+15.15662i 0.5867647-0.9882354i
#> 
#> $Trajectories$Trackway_11
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 23.13309 16.49118 Track 11.png 10 Medial 0.00             0.00
#> 2 23.63824 15.39265 Track 11.png 10 Medial 0.02             0.02
#>                polar        displacement
#> 1 23.13309+16.49118i 0.0000000+0.000000i
#> 2 23.63824+15.39265i 0.5051471-1.098529i
#> 
#> $Trajectories$Trackway_12
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 23.42206 14.58750 Track 12.png 11 Medial 0.00             0.00
#> 2 22.98750 15.70809 Track 12.png 11 Medial 0.02             0.02
#>                polar         displacement
#> 1 23.42206+14.58750i  0.0000000+0.000000i
#> 2 22.98750+15.70809i -0.4345589+1.120588i
#> 
#> $Trajectories$Trackway_13
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 14.57868 9.652942 Track 13.png 12 Medial 0.00             0.00
#> 2 14.23456 9.302207 Track 13.png 12 Medial 0.02             0.02
#> 3 13.88824 9.011030 Track 13.png 12 Medial 0.04             0.04
#> 4 13.63015 8.682354 Track 13.png 12 Medial 0.06             0.06
#> 5 13.42280 8.261030 Track 13.png 12 Medial 0.08             0.08
#>                polar          displacement
#> 1 14.57868+9.652942i  0.0000000+0.0000000i
#> 2 14.23456+9.302207i -0.3441177-0.3507353i
#> 3 13.88824+9.011030i -0.3463236-0.2911765i
#> 4 13.63015+8.682354i -0.2580883-0.3286765i
#> 5 13.42280+8.261030i -0.2073530-0.4213236i
#> 
#> $Trajectories$Trackway_14
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 19.77574 7.859559 Track 14.png 13 Medial 0.00             0.00
#> 2 19.51103 8.876471 Track 14.png 13 Medial 0.02             0.02
#> 3 19.32794 9.807354 Track 14.png 13 Medial 0.04             0.04
#>                polar          displacement
#> 1 19.77574+7.859559i  0.0000000+0.0000000i
#> 2 19.51103+8.876471i -0.2647059+1.0169118i
#> 3 19.32794+9.807354i -0.1830882+0.9308824i
#> 
#> $Trajectories$Trackway_15
#>          x         y        IMAGE ID   Side time displacementTime
#> 1 17.24780  9.004412 Track 15.png 14 Medial 0.00             0.00
#> 2 16.81324  9.968383 Track 15.png 14 Medial 0.02             0.02
#> 3 16.37868 11.082354 Track 15.png 14 Medial 0.04             0.04
#> 4 15.82059 12.092648 Track 15.png 14 Medial 0.06             0.06
#> 5 15.22941 13.089707 Track 15.png 14 Medial 0.08             0.08
#>                 polar          displacement
#> 1 17.24780+ 9.004412i  0.0000000+0.0000000i
#> 2 16.81324+ 9.968383i -0.4345589+0.9639707i
#> 3 16.37868+11.082354i -0.4345589+1.1139707i
#> 4 15.82059+12.092648i -0.5580883+1.0102942i
#> 5 15.22941+13.089707i -0.5911765+0.9970589i
#> 
#> $Trajectories$Trackway_16
#>          x         y        IMAGE ID   Side time displacementTime
#> 1 16.21103  7.497795 Track 16.png 15 Medial 0.00             0.00
#> 2 15.80515  8.309559 Track 16.png 15 Medial 0.02             0.02
#> 3 15.29559  9.121324 Track 16.png 15 Medial 0.04             0.04
#> 4 14.77059 10.027942 Track 16.png 15 Medial 0.06             0.06
#> 5 14.19485 10.963236 Track 16.png 15 Medial 0.08             0.08
#> 6 13.53971 11.757354 Track 16.png 15 Medial 0.10             0.10
#> 7 12.86471 12.573530 Track 16.png 15 Medial 0.12             0.12
#> 8 12.14118 13.277207 Track 16.png 15 Medial 0.14             0.14
#>                 polar          displacement
#> 1 16.21103+ 7.497795i  0.0000000+0.0000000i
#> 2 15.80515+ 8.309559i -0.4058824+0.8117648i
#> 3 15.29559+ 9.121324i -0.5095589+0.8117648i
#> 4 14.77059+10.027942i -0.5250000+0.9066177i
#> 5 14.19485+10.963236i -0.5757353+0.9352942i
#> 6 13.53971+11.757354i -0.6551471+0.7941177i
#> 7 12.86471+12.573530i -0.6750000+0.8161765i
#> 8 12.14118+13.277207i -0.7235295+0.7036765i
#> 
#> $Trajectories$Trackway_17
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 10.37427 12.34412 Track 17.png 16 Medial 0.00             0.00
#> 2 10.17132 12.11250 Track 17.png 16 Medial 0.02             0.02
#>                polar          displacement
#> 1 10.37427+12.34412i  0.0000000+0.0000000i
#> 2 10.17132+12.11250i -0.2029412-0.2316177i
#> 
#> $Trajectories$Trackway_18
#>          x         y        IMAGE ID   Side time displacementTime
#> 1 13.24191  8.069118 Track 18.png 17 Medial 0.00             0.00
#> 2 12.85809  9.266912 Track 18.png 17 Medial 0.02             0.02
#> 3 12.45221 10.378677 Track 18.png 17 Medial 0.04             0.04
#> 4 11.94265 11.536765 Track 18.png 17 Medial 0.06             0.06
#> 5 11.56324 12.531618 Track 18.png 17 Medial 0.08             0.08
#>                 polar         displacement
#> 1 13.24191+ 8.069118i  0.0000000+0.000000i
#> 2 12.85809+ 9.266912i -0.3838236+1.197794i
#> 3 12.45221+10.378677i -0.4058824+1.111765i
#> 4 11.94265+11.536765i -0.5095589+1.158088i
#> 5 11.56324+12.531618i -0.3794118+0.994853i
#> 
#> $Trajectories$Trackway_19
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 9.443383 5.647059 Track 19.png 18 Medial 0.00             0.00
#> 2 8.887501 7.030148 Track 19.png 18 Medial 0.02             0.02
#>                polar         displacement
#> 1 9.443383+5.647059i  0.0000000+0.000000i
#> 2 8.887501+7.030148i -0.5558824+1.383088i
#> 
#> $Trajectories$Trackway_20
#>          x         y        IMAGE ID   Side time displacementTime
#> 1 7.058824  8.338236 Track 20.png 19 Medial 0.00             0.00
#> 2 6.761030  9.505148 Track 20.png 19 Medial 0.02             0.02
#> 3 6.421324 10.548530 Track 20.png 19 Medial 0.04             0.04
#>                 polar         displacement
#> 1 7.058824+ 8.338236i  0.0000000+0.000000i
#> 2 6.761030+ 9.505148i -0.2977941+1.166912i
#> 3 6.421324+10.548530i -0.3397059+1.043382i
#> 
#> $Trajectories$Trackway_21
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 4.008089 6.498530 Track 21.png 20 Medial 0.00             0.00
#> 2 3.127941 7.217648 Track 21.png 20 Medial 0.02             0.02
#> 3 2.289706 8.064706 Track 21.png 20 Medial 0.04             0.04
#>                polar          displacement
#> 1 4.008089+6.498530i  0.0000000+0.0000000i
#> 2 3.127941+7.217648i -0.8801471+0.7191177i
#> 3 2.289706+8.064706i -0.8382353+0.8470589i
#> 
#> $Trajectories$Trackway_22
#>          x        y        IMAGE ID   Side time displacementTime
#> 1 5.272059 6.134559 Track 22.png 21 Medial 0.00             0.00
#> 2 4.883824 6.853677 Track 22.png 21 Medial 0.02             0.02
#> 3 4.552941 7.780148 Track 22.png 21 Medial 0.04             0.04
#>                polar          displacement
#> 1 5.272059+6.134559i  0.0000000+0.0000000i
#> 2 4.883824+6.853677i -0.3882353+0.7191177i
#> 3 4.552941+7.780148i -0.3308824+0.9264706i
#> 
#> $Trajectories$Trackway_23
#>           x        y        IMAGE ID   Side time displacementTime
#> 1  9.529412 12.59779 Track 23.png 22 Medial 0.00             0.00
#> 2 11.289707 12.29338 Track 23.png 22 Medial 0.02             0.02
#>                 polar        displacement
#> 1  9.529412+12.59779i 0.000000+0.0000000i
#> 2 11.289707+12.29338i 1.760294-0.3044118i
#> 
#> 
#> $Footprints
#> $Footprints[[1]]
#>           X        Y       IMAGE ID Side missing
#> 1  41.15294 15.18088 Track 1.png  0    R  Actual
#> 2  40.20441 15.82500 Track 1.png  0    L  Actual
#> 3  39.62647 16.68530 Track 1.png  0    R  Actual
#> 4  38.79706 17.23677 Track 1.png  0    L  Actual
#> 5  38.10441 18.00000 Track 1.png  0    R  Actual
#> 6  37.41177 18.66177 Track 1.png  0    L  Actual
#> 7  36.88677 19.53088 Track 1.png  0    R  Actual
#> 8  35.98677 19.86618 Track 1.png  0    L  Actual
#> 9  35.50588 20.55441 Track 1.png  0    R  Actual
#> 10 34.56177 21.19412 Track 1.png  0    L  Actual
#> 
#> $Footprints[[2]]
#>           X        Y       IMAGE ID Side missing
#> 1  39.51618 14.55441 Track 2.png  1    L  Actual
#> 2  38.86324 15.26030 Track 2.png  1    R  Actual
#> 3  38.26765 15.89559 Track 2.png  1    L  Actual
#> 4  37.45588 16.68530 Track 2.png  1    R  Actual
#> 5  36.78971 17.40441 Track 2.png  1    L  Actual
#> 6  36.14118 18.39706 Track 2.png  1    R  Actual
#> 7  35.51912 19.47794 Track 2.png  1    L  Actual
#> 8  34.96324 20.19706 Track 2.png  1    R  Actual
#> 9  34.20441 20.97353 Track 2.png  1    L  Actual
#> 10 33.70588 21.95294 Track 2.png  1    R  Actual
#> 
#> $Footprints[[3]]
#>           X        Y       IMAGE   ID Side  missing
#> 1  38.51912 15.05294 Track 3.png    2    L   Actual
#> 2  37.68971 16.05000 Track 3.png    2    R   Actual
#> 3  36.86030 16.78235 Track 3.png    2    L   Actual
#> 4  35.86324 17.95588 Track 3.png    2    R   Actual
#> 5  35.04706 18.92647 Track 3.png    2    L   Actual
#> 6  34.31471 19.67206 Track 3.png    2    R   Actual
#> NA       NA       NA Track 3.png    2    L Inferred
#> 8        NA       NA        <NA> <NA>    R   Actual
#> 
#> $Footprints[[4]]
#>          X        Y       IMAGE ID Side missing
#> 1 35.94706 15.52059 Track 4.png  3    L  Actual
#> 2 35.21471 16.50441 Track 4.png  3    R  Actual
#> 3 34.50883 17.50147 Track 4.png  3    L  Actual
#> 4 34.03677 18.45441 Track 4.png  3    R  Actual
#> 5 33.30441 19.21324 Track 4.png  3    L  Actual
#> 6 32.53236 20.06030 Track 4.png  3    R  Actual
#> 
#> $Footprints[[5]]
#>          X        Y       IMAGE ID Side missing
#> 1 32.57206 13.94559 Track 5.png  4    R  Actual
#> 2 32.00294 14.84559 Track 5.png  4    L  Actual
#> 3 31.20441 15.78530 Track 5.png  4    R  Actual
#> 
#> $Footprints[[6]]
#>          X        Y       IMAGE ID Side missing
#> 1 31.24412 14.08235 Track 6.png  5    L  Actual
#> 2 30.58235 14.88530 Track 6.png  5    R  Actual
#> 3 29.74853 15.64853 Track 6.png  5    L  Actual
#> 4 29.01618 16.50441 Track 6.png  5    R  Actual
#> 
#> $Footprints[[7]]
#>          X        Y       IMAGE ID Side missing
#> 1 31.38088 14.52794 Track 7.png  6    R  Actual
#> 2 31.02353 15.42353 Track 7.png  6    L  Actual
#> 3 30.53824 16.76912 Track 7.png  6    R  Actual
#> 4 29.95588 17.51471 Track 7.png  6    L  Actual
#> 5 29.40441 18.46765 Track 7.png  6    R  Actual
#> 6 28.71177 19.16030 Track 7.png  6    L  Actual
#> 7 28.29706 19.96324 Track 7.png  6    R  Actual
#> 
#> $Footprints[[8]]
#>          X        Y       IMAGE ID Side missing
#> 1 29.35147 14.66471 Track 8.png  7    R  Actual
#> 2 28.43824 15.55147 Track 8.png  7    L  Actual
#> 3 27.88235 16.67206 Track 8.png  7    R  Actual
#> 4 27.17647 17.67794 Track 8.png  7    L  Actual
#> 5 26.83235 18.71912 Track 8.png  7    R  Actual
#> 6 26.04265 19.49118 Track 8.png  7    L  Actual
#> 7 25.42059 20.81912 Track 8.png  7    R  Actual
#> 
#> $Footprints[[9]]
#>          X        Y       IMAGE ID Side missing
#> 1 28.50441 17.52794 Track 9.png  8    L  Actual
#> 2 29.88971 16.65441 Track 9.png  8    R  Actual
#> 3 31.18677 16.17353 Track 9.png  8    L  Actual
#> 4 32.57206 15.63088 Track 9.png  8    R  Actual
#> 5 33.74559 15.20294 Track 9.png  8    L  Actual
#> 6 34.75588 14.88530 Track 9.png  8    R  Actual
#> 
#> $Footprints[[10]]
#>          X        Y        IMAGE ID Side missing
#> 1 23.84559 16.69853 Track 10.png  9    L  Actual
#> 2 24.45441 15.59118 Track 10.png  9    R  Actual
#> 3 25.01912 14.72206 Track 10.png  9    L  Actual
#> 
#> $Footprints[[11]]
#>          X        Y        IMAGE ID Side missing
#> 1 22.86177 17.02941 Track 11.png 10    L  Actual
#> 2 23.40441 15.95294 Track 11.png 10    R  Actual
#> 3 23.87206 14.83235 Track 11.png 10    L  Actual
#> 
#> $Footprints[[12]]
#>          X        Y        IMAGE ID Side missing
#> 1 23.56765 14.05588 Track 12.png 11    L  Actual
#> 2 23.27647 15.11912 Track 12.png 11    R  Actual
#> 3 22.69853 16.29706 Track 12.png 11    L  Actual
#> 
#> $Footprints[[13]]
#>          X        Y        IMAGE ID Side missing
#> 1 14.73088 9.798530 Track 13.png 12    L  Actual
#> 2 14.42647 9.507354 Track 13.png 12    R  Actual
#> 3 14.04265 9.097059 Track 13.png 12    L  Actual
#> 4 13.73382 8.925001 Track 13.png 12    R  Actual
#> 5 13.52647 8.439706 Track 13.png 12    L  Actual
#> 6 13.31912 8.082353 Track 13.png 12    R  Actual
#> 
#> $Footprints[[14]]
#>          X         Y        IMAGE ID Side missing
#> 1 20.00735  7.288236 Track 14.png 13    R  Actual
#> 2 19.54412  8.430883 Track 14.png 13    L  Actual
#> 3 19.47794  9.322059 Track 14.png 13    R  Actual
#> 4 19.17794 10.292648 Track 14.png 13    L  Actual
#> 
#> $Footprints[[15]]
#>          X         Y        IMAGE ID Side missing
#> 1 17.47941  8.602942 Track 15.png 14    R  Actual
#> 2 17.01618  9.405883 Track 15.png 14    L  Actual
#> 3 16.61030 10.530883 Track 15.png 14    R  Actual
#> 4 16.14706 11.633824 Track 15.png 14    L  Actual
#> 5 15.49412 12.551471 Track 15.png 14    R  Actual
#> 6 14.96471 13.627942 Track 15.png 14    L  Actual
#> 
#> $Footprints[[16]]
#>          X         Y        IMAGE ID Side missing
#> 1 16.37206  7.072059 Track 16.png 15    L  Actual
#> 2 16.05000  7.923530 Track 16.png 15    R  Actual
#> 3 15.56030  8.695589 Track 16.png 15    L  Actual
#> 4 15.03088  9.547059 Track 16.png 15    R  Actual
#> 5 14.51030 10.508824 Track 16.png 15    L  Actual
#> 6 13.87941 11.417648 Track 16.png 15    R  Actual
#> 7 13.20000 12.097060 Track 16.png 15    L  Actual
#> 8 12.52941 13.050001 Track 16.png 15    R  Actual
#> 9 11.75294 13.504413 Track 16.png 15    L  Actual
#> 
#> $Footprints[[17]]
#>          X        Y        IMAGE ID Side missing
#> 1 10.46029 12.46765 Track 17.png 16    R  Actual
#> 2 10.28824 12.22059 Track 17.png 16    L  Actual
#> 3 10.05441 12.00441 Track 17.png 16    R  Actual
#> 
#> $Footprints[[18]]
#>          X         Y        IMAGE ID Side missing
#> 1 13.53971  7.442648 Track 18.png 17    R  Actual
#> 2 12.94412  8.695589 Track 18.png 17    L  Actual
#> 3 12.77206  9.838236 Track 18.png 17    R  Actual
#> 4 12.13235 10.919118 Track 18.png 17    L  Actual
#> 5 11.75294 12.154413 Track 18.png 17    R  Actual
#> 6 11.37353 12.908824 Track 18.png 17    L  Actual
#> 
#> $Footprints[[19]]
#>          X        Y        IMAGE ID Side missing
#> 1 9.767648 4.844118 Track 19.png 18    L  Actual
#> 2 9.119118 6.450000 Track 19.png 18    R  Actual
#> 3 8.655883 7.610295 Track 19.png 18    L  Actual
#> 
#> $Footprints[[20]]
#>          X         Y        IMAGE ID Side missing
#> 1 7.191177  7.733824 Track 20.png 19    R  Actual
#> 2 6.926471  8.942648 Track 20.png 19    L  Actual
#> 3 6.595589 10.067648 Track 20.png 19    R  Actual
#> 4 6.247059 11.029412 Track 20.png 19    L  Actual
#> 
#> $Footprints[[21]]
#>          X        Y        IMAGE ID Side missing
#> 1 4.411765 6.220589 Track 21.png 20    R  Actual
#> 2 3.604412 6.776471 Track 21.png 20    L  Actual
#> 3 2.651471 7.658824 Track 21.png 20    R  Actual
#> 4 1.927941 8.470589 Track 21.png 20    L  Actual
#> 
#> $Footprints[[22]]
#>          X        Y        IMAGE ID Side missing
#> 1 5.580883 5.889706 Track 22.png 21    R  Actual
#> 2 4.963236 6.379412 Track 22.png 21    L  Actual
#> 3 4.804412 7.327942 Track 22.png 21    R  Actual
#> 4 4.301471 8.232353 Track 22.png 21    L  Actual
#> 
#> $Footprints[[23]]
#>           X        Y        IMAGE ID Side missing
#> 1  8.505883 12.67941 Track 23.png 22    R  Actual
#> 2 10.552942 12.51618 Track 23.png 22    L  Actual
#> 3 12.026471 12.07059 Track 23.png 22    R  Actual
#> 
#> 
```
