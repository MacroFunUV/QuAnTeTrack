# Subset trackways

`subset_track()` is a function that subsets trackways from a `trackway`
R object based on the specified indices.

## Usage

``` r
subset_track(data, tracks = NULL)
```

## Arguments

- data:

  A `trackway` R object, which is a list consisting of two elements:

  - **`Trajectories`**: A list of trajectories representing trackway
    midlines, interpolated by connecting the midpoints of successive
    left–right footprint pairs (i.e., footprints linked by pace lines).
    Includes columns `X`, `Y`, `IMAGE`, `ID`, and `Side` (set to
    `"Medial"`).

  - **`Footprints`**: A list of data frames containing footprint
    coordinates and associated metadata, with a `Side` column (`"R"` or
    `"L"`) and a `missing` marker (`"Actual"` or `"Inferred"`).

- tracks:

  A numeric vector specifying the indices of trackways to subset. The
  default is to include all trackways.

## Value

A `trackway` R object that contains only the specified subset of
trackways. The structure of the returned object mirrors the input
structure but includes only the selected trackways.

## Details

This function subsets both the **`Trajectories`** and **`Footprints`**
elements of the input data based on the provided vector of indices. It
allows users to focus on a specific subset of trackways for further
analysis or visualization, particularly when working with large datasets
containing numerous trackways.

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
# Example 1: Subset the first three trackways of MountTom dataset.
subset_data <- subset_track(MountTom, tracks = c(1:3))

# Example 2:  Subset the trackways at indices 5, 7, and 10.
subset_data <- subset_track(MountTom, tracks = c(5, 7, 10))
```
