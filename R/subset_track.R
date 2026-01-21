#' Subset trackways
#'
#' \code{subset_track()} is a function that subsets trackways from a \code{trackway} R object based on the specified indices.
#'
#' @param data A \code{trackway} R object, which is a list consisting of two elements:
#'
#'    * \strong{\code{Trajectories}}: A list of trajectories representing trackway midlines,
#'      interpolated by connecting the midpoints of successive left–right footprint pairs
#'      (i.e., footprints linked by pace lines). Includes columns \code{X}, \code{Y},
#'      \code{IMAGE}, \code{ID}, and \code{Side} (set to \code{"Medial"}).
#'
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates
#'      and associated metadata, with a \code{Side} column (\code{"R"} or \code{"L"})
#'      and a \code{missing} marker (\code{"Actual"} or \code{"Inferred"}).
#' @param tracks A numeric vector specifying the indices of trackways to subset. The default is to include all trackways.
#'
#' @details This function subsets both the \strong{\code{Trajectories}} and \strong{\code{Footprints}} elements of the input data based on the provided vector of indices.
#' It allows users to focus on a specific subset of trackways for further analysis or visualization, particularly when working with large datasets containing numerous trackways.
#'
#' @return A \code{trackway} R object that contains only the specified subset of trackways. The structure of the returned object mirrors the input structure but includes only the selected trackways.
#'
#' @section Logo:
#' \if{html}{\figure{Logo.png}{options: width=120}}
#'
#' @author Humberto G. Ferrón
#' @author humberto.ferron@uv.es
#' @author Macroevolution and Functional Morphology Research Group (www.macrofun.es)
#' @author Cavanilles Institute of Biodiversity and Evolutionary Biology
#' @author Calle Catedrático José Beltrán Martínez, nº 2
#' @author 46980 Paterna - Valencia - Spain
#' @author Phone: +34 (9635) 44477
#'
#' @examples
#' # Example 1: Subset the first three trackways of MountTom dataset.
#' subset_data <- subset_track(MountTom, tracks = c(1:3))
#'
#' # Example 2:  Subset the trackways at indices 5, 7, and 10.
#' subset_data <- subset_track(MountTom, tracks = c(5, 7, 10))
#'
#' @seealso \code{\link{tps_to_track}}
#'
#' @export

subset_track <- function(data, tracks = NULL) {
  ## Errors and Warnings----

  # Check if 'data' is a list with at least two elements
  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'trackway' R object, which is a list consisting of two elements.")
  }

  # Check if the two elements of 'data' are lists
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("The two elements of 'data' must be lists.")
  }

  # Check if 'tracks' is provided and is not NULL, it must be a numeric vector
  if (!is.null(tracks) && (!is.numeric(tracks) || any(tracks <= 0))) {
    stop("The 'tracks' argument must be a numeric vector of positive indices.")
  }

  # Warning if 'tracks' contains indices that are out of bounds
  if (!is.null(tracks) && any(tracks > length(data[[1]]))) {
    warning("The 'tracks' argument contains indices that exceed the length of the data. These indices will be ignored.")
    tracks <- tracks[tracks <= length(data[[1]])]
  }

  ## Code----

  # Set default values if arguments are NULL
  if (is.null(tracks)) {
    tracks <- c(1:length(data[[1]]))
  } # If 'tracks' is NULL, select all elements in the first element of 'data'

  # Subset the first element of 'data' using 'tracks'
  data[[1]] <- data[[1]][tracks]

  # Subset the second element of 'data' using 'tracks'
  data[[2]] <- data[[2]][tracks]

  # Return the subsetted 'data'
  return(data)
}
