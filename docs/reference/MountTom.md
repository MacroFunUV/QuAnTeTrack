# MountTom Dinosaur Track Dataset

A 'track' R object representing dinosaur tracks from the Mount Tom site.

## Usage

``` r
MountTom
```

## Format

A list consisting of two elements:

- **Trajectories**: A list of interpolated trajectories, where each
  trajectory is a series of midpoints between consecutive footprints.

- **Footprints**: A list of data frames containing footprint
  coordinates, metadata (e.g., image reference, ID), and a marker
  indicating whether the footprint is actual or inferred.

## Source

Ostrom, J. H. (1972). Were some dinosaurs gregarious?. Palaeogeography,
Palaeoclimatology, Palaeoecology, 11(4), 287-301.
