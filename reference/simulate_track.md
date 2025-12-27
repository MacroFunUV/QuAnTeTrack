# Simulate tracks using different models

`simulate_track()` simulates movement trajectories based on an original
set of tracks. Three movement models are available for simulation, each
reflecting different levels of constraint in movement patterns. These
models can represent biological or environmental constraints, such as
movement along coastlines, rivers, or towards resources like water or
food.

## Usage

``` r
simulate_track(data, nsim = NULL, model = NULL)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- nsim:

  The number of simulations to run. Defaults to `1000` if not specified.

- model:

  The type of movement model to use. Options are `"Directed"`,
  `"Constrained"`, or `"Unconstrained"`. Defaults to `"Unconstrained"`
  if not provided.

## Value

A `track simulation` R object consisting of a list of simulated
trajectories stored as `track` R objects.

## Details

This function simulates movement trajectories based on the following
models:

- **Directed**: This model simulates movement that follows a specific
  direction navigating with a compass (i.e., a directed walk/allothetic
  directed walk/oriented path) (Cheung et al., 2007, 2008). The
  trajectory is constrained by both the angular and linear properties of
  the original track, with minor deviations allowed to reflect natural
  variability.

  - **Angular constraints**: The trajectory closely follows a specific
    direction, maintaining the overall angular orientation of the
    original track. Deviations of consecutive steps are governed by the
    angular standard deviation calculated from the original track using
    `TrajAngles()`.

  - **Linear constraints**: Step lengths are constrained to the mean of
    the original track's step lengths, with variability allowed
    according to the standard deviation of step lengths computed with
    `TrajStepLengths()`.

  - **Starting direction**: Fixed to the original direction (overall
    angular orientation) of the track.

This model is ideal for simulating movement directed toward a specific
resource or constrained by natural barriers, with a relatively fixed
direction and minor deviations.

- **Constrained**: This model simulates movement that correspond to a
  correllated random walk/idiothetic directed walk (Kareiva & Shigesada,
  1983), corresponding to an animal navigating without a compass (Cheung
  et al., 2008), while still maintaining certain angular and linear
  characteristics of the original track. It provides more flexibility
  than the Directed model but is not entirely random like the
  Unconstrained model.

  - **Angular constraints**: The trajectory does not follow a specific
    direction. Deviations of consecutive steps are governed by the
    angular standard deviation calculated from the original track using
    `TrajAngles()`.

  - **Linear constraints**: Step lengths are constrained to the mean of
    the original track's step lengths, with variability allowed
    according to the standard deviation of step lengths computed with
    `TrajStepLengths()`.

  - **Starting direction**: Fixed to the original direction (overall
    angular orientation) of the track.

This model is suitable for scenarios where movement is influenced by
external constraints but allows for some degree of random exploration.

- **Unconstrained**: This model simulates movement that correspond to a
  correllated random walk/idiothetic directed walk (Kareiva & Shigesada,
  1983), corresponding to an animal navigating without a compass (Cheung
  et al., 2008), while still maintaining certain angular and linear
  characteristics of the original track.

  - **Angular constraints**: The trajectory does not follow a specific
    direction. Deviations of consecutive steps are governed by the
    angular standard deviation calculated from the original track using
    `TrajAngles()`.

  - **Linear constraints**: Step lengths are constrained to the mean of
    the original track's step lengths, with variability allowed
    according to the standard deviation of step lengths computed with
    `TrajStepLengths()`.

  - **Starting direction**: Randomly determined.

This model is suitable for simulating exploratory or dispersal behavior
in open environments, where movement is random and not influenced by
specific constraints.

Note: Simulations cannot be applied to trajectories with fewer than four
steps as the standard deviations of angles and step lengths cannot be
computed for such short trajectories. Consider using the
[`subset_track()`](https://macrofunuv.github.io/QuAnTeTrack/reference/subset_track.md)
function to filter tracks with four or more steps.

The function utilizes the trajr package for key calculations:

- `TrajGenerate()`: Generates a new trajectory based on random or
  directed movement models, constrained by specified parameters.

- `TrajStepLengths()`: Calculates the step lengths (distances between
  consecutive points) of the original trajectory.

- `TrajAngles()`: Computes the angles between consecutive segments of
  the trajectory, used to maintain directional movement in constrained
  models.

- `TrajRotate()`: Rotates the trajectory by a specified angle to match
  the original direction or a random angle for unconstrained models.

- `TrajTranslate()`: Translates the simulated trajectory to start at the
  same geographic location as the original.

The `NISTdegTOradian()` function from the NISTunits package is used to
convert angles from degrees to radians.

## Logo

![](figures/Logo.png)

## References

Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007). Animal
navigation: the difficulty of moving in a straight line. Biological
cybernetics, 97, 47-61.

Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2008). Animal
navigation: general properties of directed walks. Biological
cybernetics, 99, 197-217.

## See also

[`tps_to_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md),
[`plot_sim`](https://macrofunuv.github.io/QuAnTeTrack/reference/plot_sim.md),
[`subset_track`](https://macrofunuv.github.io/QuAnTeTrack/reference/subset_track.md),
[`TrajGenerate`](https://rdrr.io/pkg/trajr/man/TrajGenerate.html),
[`TrajStepLengths`](https://rdrr.io/pkg/trajr/man/TrajStepLengths.html),
[`TrajAngles`](https://rdrr.io/pkg/trajr/man/TrajAngles.html),
[`TrajRotate`](https://rdrr.io/pkg/trajr/man/TrajRotate.html),
[`TrajTranslate`](https://rdrr.io/pkg/trajr/man/TrajTranslate.html),
[`NISTdegTOradian`](https://rdrr.io/pkg/NISTunits/man/NISTdegTOradian.html)

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
# Example 1: Simulate tracks using data from the Paluxy River
# Default model (Unconstrained movement)
simulated_tracks <- simulate_track(PaluxyRiver, nsim = 3)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.

# Example 2: Simulate tracks using the "Directed" model, representing movement
# toward a resource (e.g., water source)
simulated_tracks_directed <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")

# Example 3: Simulate tracks using the "Constrained" model, representing movement
# along a geographic feature (e.g., coastline)
simulated_tracks_constrained <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")

# Example 4: Simulate tracks using the "Unconstrained" model (random exploratory movement)
simulated_tracks_unconstrained <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")

# Subsetting trajectories with four or more steps in the MountTom dataset
sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))

# Example 5: Simulate tracks using data from Mount Tom
# Default model (Unconstrained movement)
simulated_tracks_mt <- simulate_track(sbMountTom, nsim = 3)
#> Warning: `model` is NULL. Defaulting to 'Unconstrained'.

# Example 6: Simulate tracks using the "Directed" model for Mount Tom, representing
# directed movement
simulated_tracks_mt_directed <- simulate_track(sbMountTom, nsim = 3, model = "Directed")

# Example 7: Simulate tracks using the "Constrained" model for Mount Tom, representing
# constrained movement
simulated_tracks_mt_constrained <- simulate_track(sbMountTom, nsim = 3, model = "Constrained")

# Example 8: Simulate tracks using the "Unconstrained" model for Mount Tom, representing
# random exploratory movement
simulated_tracks_mt_unconstrained <- simulate_track(sbMountTom, nsim = 3, model = "Unconstrained")
```
