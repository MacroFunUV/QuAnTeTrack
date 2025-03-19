
<!-- README.md is generated from README.Rmd. Please edit that file -->

# QuAnTeTrack

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/QuAnTeTrack)](https://CRAN.R-project.org/package=QuAnTeTrack)
[![R-CMD-check](https://github.com/MacroFunUV/QuAnTeTrack/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MacroFunUV/QuAnTeTrack/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Studying fossil vertebrate trackways is crucial as they provide valuable
insights into the behavior, locomotion, and environmental interactions
of extinct species, offering a dynamic perspective on prehistoric life
that skeletal remains alone cannot. Here, I present **QuAnTeTrack**, an
integrated tool designed to **semi-automatically extract
palaeobiological information** by digitizing footprint and trackway
coordinates and analyzing them within the R environment.

QuAnTeTrack includes functions to read footprint coordinates and
trackways stored as **.TPS files**, transforming them into **“track”
objects**. These objects can be used to **calculate and plot key
parameters** such as **turning angles, trackway distances, step lengths,
sinuosity, and straightness**. Additional functions enable the
**quantification of velocity** and the testing of scenarios involving
**acceleration, deceleration, or stable velocity**, as well as comparing
velocity differences between tracks.

The tool also supports **visualizing movement patterns**, plotting
tracks according to velocity. Furthermore, **QuAnTeTrack facilitates the
testing of hypotheses related to gregarious behavior or predatory
events**. For interactions involving multiple tetrapods, such as **group
movement or predation**, trackways are expected to covary more
significantly and present fewer intersections than those generated by
independent events. QuAnTeTrack estimates this by computing **Frechet
and Dynamic Time Warping (DTW) trajectory similarity metrics** and
**quantifying track intersections**. These metrics are evaluated through
**random simulation procedures** to determine how actual trackway
similarities differ from randomly generated ones under various
scenarios, including **geographical constraints and resource
directionality**.

## Installation

You can install the development version of QuAnTeTrack from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MacroFunUV/QuAnTeTrack")
```

## Example

``` r
# Load package QuAnTeTrack
=======

``` load
library(QuAnTeTrack)
```

# Load data (Replace with your file paths)

PaluxyRiver \<- tps_to_track(“path/to/PaluxyRiver.tps”) MountTom \<-
tps_to_track(“path/to/MountTom.tps”)

# Extract track parameters

track_params_paluxy \<- track_param(PaluxyRiver) track_params_mt \<-
track_param(MountTom)

# Calculate velocities

H_paluxyriver \<- c(3.472, 2.200) H_mounttom \<- c(1.380, 1.404, 1.320,
1.736, 1.364, 1.432, 1.508, 1.768, 1.600, 1.848, 1.532, 1.532, 0.760,
1.532, 1.688, 1.620, 0.636, 1.784, 1.676, 1.872, 1.648, 1.760, 1.612)

vel_paluxy \<- velocity_track(PaluxyRiver, H = H_paluxyriver) vel_mt \<-
velocity_track(MountTom, H = H_mounttom)

# Plot track and direction

plot_track(PaluxyRiver) plot_direction(PaluxyRiver)

# Test directionality (ANOVA)

test_dir_paluxy \<- test_direction(PaluxyRiver, analysis = “ANOVA”)
test_dir_mt \<- test_direction(MountTom, analysis = “ANOVA”)

# Plot velocities

plot_velocity(PaluxyRiver, vel_paluxy) plot_velocity(MountTom, vel_mt)

# Test velocity differences (Kruskal-Wallis)

test_vel_paluxy \<- test_velocity(PaluxyRiver, vel_paluxy, analysis =
“Kruskal-Wallis”, plot = TRUE) test_vel_mt \<- test_velocity(MountTom,
vel_mt, analysis = “Kruskal-Wallis”, plot = TRUE)

# Simulate tracks using the “Directed” model

sim_paluxy \<- simulate_track(PaluxyRiver, nsim = 1000, model =
“Directed”) plot_sim(sim_paluxy, PaluxyRiver)

# Compare track similarity using DTW

dtw_similarity \<- simil_DTW_metric(PaluxyRiver, PaluxyRiver)
frechet_similarity \<- simil_Frechet_metric(PaluxyRiver, PaluxyRiver)

# Calculate intersections

intersection_metrics \<- track_intersection(PaluxyRiver, test = TRUE,
sim = sim_paluxy, origin.permutation = “None”)

# Combine probabilities

combined_prob \<- combined_prob(dtw_similarity, intersection_metrics)

# Cluster tracks

clusters \<- cluster_track(dtw_similarity)
