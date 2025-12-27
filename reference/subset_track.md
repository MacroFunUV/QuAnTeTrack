# Subset tracks

`subset_track()` is a function that subsets tracks from a list of track
data based on the specified indices.

## Usage

``` r
subset_track(data, tracks = NULL)
```

## Arguments

- data:

  A `track` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of interpolated trajectories, where each
    trajectory is a series of midpoints between consecutive footprints.

  - **`Footprints`**: A list of data frames containing footprint
    coordinates, metadata (e.g., image reference, ID), and a marker
    indicating whether the footprint is actual or inferred.

- tracks:

  A numeric vector specifying the indices of tracks to subset. The
  default is to include all tracks.

## Value

A `track` R object that contains only the specified subset of tracks.
The structure of the returned object mirrors the input structure but
includes only the selected tracks.

## Details

This function subsets both the **`Trajectories`** and **`Footprints`**
elements of the input data based on the provided vector of indices. It
allows users to focus on a specific subset of tracks for further
analysis or visualization, particularly when working with large datasets
containing numerous tracks.

## Logo

![](figures/Logo.png)

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
# Example 1: Subset the first three tracks of MountTom dataset.
subset_data <- subset_track(MountTom, tracks = c(1:3))

# Example 2:  Subset the tracks at indices 5, 7, and 10.
subset_data <- subset_track(MountTom, tracks = c(5, 7, 10))
```
