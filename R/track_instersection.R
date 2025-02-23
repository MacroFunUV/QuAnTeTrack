#' Calculate intersection metrics in tracks
#'
#' \code{track_intersection()} calculates the number of unique intersections between trajectories.
#' The function also supports testing with simulations and different permutation procedures for the coordinates
#' of the simulated trajectories' origins to compute p-values. This allows for a robust assessment of the intersection metrics,
#' enabling users to evaluate the significance of the observed intersections in relation to simulated trajectories.
#'
#' @param data A 'track' R object, which is a list consisting of two elements:
#'    * \strong{Trajectories}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{Footprints}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param test Logical; if \code{TRUE}, the function compares the observed DTW distances against. Default is \code{FALSE}.
#' @param sim A 'track simulation' R object consisting of a list of simulated trajectories to use for comparison when \code{test = TRUE}.
#' @param origin.permutation A character string specifying the method for permutation of the coordinates of the simulated trajectories' origins.
#' Options include "None", "Min.Box", "Conv.Hull", or "Custom". Default is "None".
#' @param custom.coord A matrix of custom coordinates that define the vertices of an area for permutation of the coordinates of the simulated
#' trajectories' origins.
#'
#' @details
#' The \code{track_instersection()} function calculates the number of unique intersections between trajectories.
#' The function also provides options for conducting hypothesis testing through simulated data with permutations of simulated
#' trajectory origins, enabling the calculation of p-values to evaluate the significance of observed intersections.
#'
#' The \code{origin.permutation} parameter determines whether any permutation will occur for the simulated trajectories
#' and specifies the method used to permute the starting coordinates of these trajectories. The available options
#' include \code{"None"}, which indicates that no permutation will be applied, allowing the function to compute
#' intersections based solely on the simulated trajectories starting at the same coordinates as the original data.
#' The \code{"Min.Box"} option conducts permutation within the minimum bounding box surrounding the original
#' coordinates of origin. Alternatively, the \code{"Conv.Hull"} option utilizes the convex hull around the original
#' coordinates of origin to perform permutations. Lastly, the \code{"Custom"} option allows for permutation based
#' on user-defined coordinates that delimit a specific area of interest, as specified in the \code{custom.coord} parameter.
#' The \code{custom.coord} parameter must be a matrix specifying custom coordinates representing the vertices of the desired area,
#' effectively constraining the random placement of trajectory origins to this defined region.
#'
#' Intersections between trajectories can provide valuable insights into (palaeo)ethological behaviors,
#' such as group movement or hunting dynamics. In cases of side-by-side movement, like a group walking
#' in parallel or a hunting pack, fewer intersections in the actual data  would be expected compared to random
#' scenarios, as coordinated movements reduce the likelihood of trajectories crossing. Similarly,
#' in a chasing or hunting scene or queuing, where one trackmaker is moving ahead of the other, the paths may intersect more frequently than
#' random movement would predict, leading to an increase in intersections. In both cases, significant similarity metrics may also emerge,
#' as these structured, non-random movement patterns deviate considerably from those in random simulations. Therefore, combining results
#' from intersection counts and similarity metrics from \code{simil_DTW_metric} and \code{simil_Frechet_metric} functions might be advisable
#' when testing specific behavioral hypotheses.
#'
#' The application of permutation to the origins of simulated trajectories, and the method selected, depends on the hypothesis being tested.
#' If the objective is to determine how many intersections should be expected when trackmakers originate from specific points,
#' and compare this to the actual data, the origin points should be constrained by the original coordinates, making the \code{"None"}
#' option the most appropriate. Conversely, if the hypothesis involves testing how many intersections would occur if trackmakers
#' passed through a broader or defined area, the \code{"Min.Box"}, \code{"Conv.Hull"}, or \code{"Custom"} options are more suitable.
#'
#' The choice of permutation method also relies on available information about the spatial region where tracks might have been recorded.
#' For example, the \code{"Min.Box"} or \code{"Conv.Hull"} methods may be employed when there is a defined region that reflects
#' the likely constraints of where trackmakers left their traces. The \code{"Custom"} option is particularly useful when there is specific
#' knowledge of the area (e.g., terrain features or environmental conditions) and permutations need to be restricted or excluded from certain zones.
#'
#' Ultimately, the selection of an appropriate permutation method (or none at all) should align with the behavioral assumptions
#' underlying the hypothesis—whether testing for gregarism, parallel movement, hunting strategies, or random movement—and
#' how these behaviors are expected to manifest in the number of intersections observed in both the actual and simulated datasets.
#'
#' @return A 'track intersection' R object consisting of a list containing the following elements:
#' \item{Intersection_metric}{A matrix of unique intersection counts between trajectories. Each entry
#' represents the number of unique intersection points between the corresponding pair of trajectories.}
#' \item{Intersection_metric_p_values}{(If `test` is TRUE) A matrix of p-values associated with
#' the intersection metrics, calculated through permutations of simulated trajectory origins. Each entry
#' reflects the probability of observing an intersection count as extreme as the observed one,
#' given the null hypothesis of no difference.}
#' \item{Intersection_metric_p_values_combined}{(If `test` is TRUE) A numeric value representing
#' the combined p-value for all intersections, indicating the overall significance of the
#' intersection metrics across all pairs of trajectories.}
#' \item{Intersection_metric_simulations}{(If `test` is TRUE) A list containing matrices of
#' intersection counts for each simulation iteration, allowing for further inspection of the
#' distribution of intersections across multiple randomized scenarios.}
#'
#'@section Logo:
#'\if{html}{\figure{Logo.png}{options: width=30\%}}
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
#' # Example 1: Simulating tracks and comparing intersection metrics in the PaluxyRiver dataset.
#' # No origin permutation is applied ("None").
#' s1 <- simulate_track(PaluxyRiver, nsim = 1000, model = "Directed")
#' int1 <- track_instersection(PaluxyRiver, test = TRUE, sim = s1, origin.permutation = "None")
#' print(int1)
#'
#' # Example 2: Simulating tracks and comparing intersection metrics in the PaluxyRiver dataset.
#' # The origin permutation is applied using the minimum bounding box ("Min.Box").
#' s2 <- simulate_track(PaluxyRiver, nsim = 1000, model = "Constrained")
#' int2 <- track_instersection(PaluxyRiver, test = TRUE, sim = s2, origin.permutation = "Min.Box")
#' print(int2)
#'
#' # Example 3: Simulating tracks and comparing intersection metrics in the PaluxyRiver dataset.
#' # The origin permutation is applied using the convex hull ("Conv.Hull").
#' s3 <- simulate_track(PaluxyRiver, nsim = 1000, model = "Unconstrained")
#' int3 <- track_instersection(PaluxyRiver, test = TRUE, sim = s3, origin.permutation = "Conv.Hull")
#' print(int3)
#'
#' # Example 4: Simulating tracks and comparing intersection metrics in a subsample of the MoutTom dataset.
#' # The "Min.Box" origin permutation is applied.
#' sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
#' s4 <- simulate_track(sbMountTom, nsim = 1000)
#' int4 <- track_instersection(sbMountTom, test = TRUE, sim = s4, origin.permutation = "Min.Box")
#' print(int4)
#'
#' # Example 5: Simulating tracks and comparing intersection metrics in a subsample of the MoutTom dataset.
#' # The origin permutation is customized ("Custom") using predefined coordinates.
#' sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
#' s5 <- simulate_track(sbMountTom, nsim = 1000)
#' area_origin <- matrix(c(50, 5,
#'                         10, 5,
#'                         10, 20,
#'                         50, 20),
#'                       ncol = 2, byrow = TRUE)
#' int5 <- track_instersection(sbMountTom, test = TRUE, sim = s5, origin.permutation = "Custom", custom.coord = area_origin)
#' print(int5)
#'
#' @importFrom grDevices chull
#' @importFrom shotGroups getMinBBox
#' @importFrom splancs csr
#' @importFrom trajr TrajTranslate
#'
#' @seealso \code{\link[tps_to_track]{tps_to_track}}, \code{\link[simulate_track]{simulate_track}}, \code{\link[simil_DTW_metric]{simil_DTW_metric}}, \code{\link[simil_Frechet_metric]{simil_Frechet_metric}}
#'
#' @export


track_instersection <- function(data, test = NULL, sim = NULL, origin.permutation = NULL, custom.coord = NULL) {

  ## Set default values if arguments are NULL----
  if (is.null(test)) test <- FALSE  # Set default if 'test' is NULL
  if (is.null(origin.permutation)) origin.permutation <- "None"  # Set default if 'origin.permutation' is NULL

  ## Errors and Warnings----

  # Check if 'data' is a list with at least two elements
  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'track' R object, which is a list consisting of two elements.")
  }

  # Check if the two elements of 'data' are lists
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("The two elements of 'data' must be lists.")
  }

  # Warn if the 'test' argument is not a boolean
  if (!is.logical(test)) {
    stop("'test' argument should be TRUE or FALSE.")
  }

  # Check if 'sim' is provided when test is TRUE
  if (test == TRUE && is.null(sim)) {
    stop("A 'sim' argument must be provided when 'test' is TRUE.")
  }

  # If 'sim' is provided, ensure it is a list and has the same structure as 'data'
  if (!is.null(sim)) {
    if (!is.list(sim)) {
      stop("The 'sim' argument must be a list.")
    }

    # Check that 'sim' contains the same number of tracks as 'data'
    if (length(sim[[1]]) != length(data[[1]])) {
      stop("The 'sim' list must have the same number of trajectories as 'data'.")
    }
  }

  # Check if 'origin.permutation' is valid
  valid_permutations <- c("None", "Min.Box", "Conv.Hull", "Custom")
  if (!origin.permutation %in% valid_permutations) {
    stop(paste("Invalid 'origin.permutation'. Valid options are:", paste(valid_permutations, collapse = ", ")))
  }

  # If 'origin.permutation' is "Custom", check if 'custom.coord' is provided
  if (origin.permutation == "Custom" && is.null(custom.coord)) {
    stop("If 'origin.permutation' is set to 'Custom', the 'custom.coord' must be provided.")
  }

  # Check if 'custom.coord' is a matrix or a data frame with two columns
  if (!is.null(custom.coord) && !is.matrix(custom.coord) && !is.data.frame(custom.coord)) {
    stop("The 'custom.coord' must be a matrix or a data frame.")
  }

  if (!is.null(custom.coord) && ncol(as.matrix(custom.coord)) != 2) {
    stop("The 'custom.coord' must have exactly two columns.")
  }


  ## Code----
  data <- data[[1]]

  # Function to determine if two line segments intersect
  find_intersection <- function(p1, p2, p3, p4) {
    x1 <- p1[1]; y1 <- p1[2]
    x2 <- p2[1]; y2 <- p2[2]
    x3 <- p3[1]; y3 <- p3[2]
    x4 <- p4[1]; y4 <- p4[2]

    # Compute determinants
    denom <- (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

    # Check if the lines are parallel
    if (denom == 0) {
      return(c(NA, NA))  # Lines are parallel or coincident
    }

    # Compute intersection point
    intersect_x <- ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom
    intersect_y <- ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom

    # Check if the intersection point is within both line segments
    if (min(x1, x2) <= intersect_x && intersect_x <= max(x1, x2) &&
        min(y1, y2) <= intersect_y && intersect_y <= max(y1, y2) &&
        min(x3, x4) <= intersect_x && intersect_x <= max(x3, x4) &&
        min(y3, y4) <= intersect_y && intersect_y <= max(y3, y4)) {
      return(c(intersect_x, intersect_y))
    } else {
      return(c(NA, NA))
    }
  }

  # Function to count unique intersections between two trajectories
  intersect <- function(traj1, traj2) {
    intersections <- list()

    for (i in 1:(length(traj1$x) - 1)) {
      for (j in 1:(length(traj2$x) - 1)) {
        int <- find_intersection(
          c(traj1$x[i], traj1$y[i]), c(traj1$x[i+1], traj1$y[i+1]),
          c(traj2$x[j], traj2$y[j]), c(traj2$x[j+1], traj2$y[j+1])
        )

        if (!is.na(int[1])) {
          intersections[[length(intersections) + 1]] <- int
        }
      }
    }

    # Remove duplicate intersections
    unique_intersections <- unique(do.call(rbind, intersections))

    # Return the number of unique intersections
    if (is.null(unique_intersections)){
      return(0)} else {
        return(nrow(unique_intersections))}
  }


  # Calculate actual intersection metrics
  Matrixsim<-data.frame(matrix(nrow=length(data),ncol = length(data)))
  colnames(Matrixsim)<-names(data)
  rownames(Matrixsim)<-names(data)

  Intersect<-Matrixsim
  for (i in 1:length(data)){
    for (j in 1:length(data)){
      if(i<=j) next
      Intersect[j,i]<-intersect(data[[i]],data[[j]])
    }
  }


  # Calculate simulation intersection metrics
  # Permutation of coordinates at origin
  if (test == TRUE) {

    nsim <- length(sim)

    if (origin.permutation != "None") {

      mat <- matrix(ncol = 2, nrow = length(data))

      for (i in 1:length(data)) {
        mat[, 1][i] <- data[[i]]$x[1]
        mat[, 2][i] <- data[[i]]$y[1]
      }

      if (length(data) == 2) {
        mat <- rbind(mat, mat * 1.000001)
      }

      if (origin.permutation == "Min.Box") {
        mat <- getMinBBox(mat)
        mat <- mat$pts
      }

      if (origin.permutation == "Conv.Hull") {
        mat <- mat[chull(mat), ]
      }

      if (origin.permutation == "Custom") {
        mat <- custom.coord
      }

      for (i in 1:length(sim)) {
        for (j in 1:length(sim[[1]])) {
          sim[[i]][[j]] <- TrajTranslate(sim[[i]][[j]], csr(mat, 1)[1], csr(mat, 1)[2])
        }

        writeLines(paste(Sys.time(), paste("Permutation", i)))
        writeLines(" ")
        writeLines(paste("Permutation of coordinates at origin using", origin.permutation))
        writeLines("------------------------------------")
        if (i == nsim) {
          writeLines("PERMUTATION COMPLETED")
          writeLines("------------------------------------")
          writeLines(" ")
        }
      }
    }

    # Calculate metrics
    Intersectsim <- Matrixsim

    listIntersect <- list()
    listnegIntersect <- c()

    for (i in 1:nsim) {
      for (c in 1:length(data)) {
        for (r in 1:length(data)) {
          if (c <= r) next
          Intersectsim[r, c] <- intersect(sim[[i]][[r]], sim[[i]][[c]])
        }
      }
      listIntersect[[i]] <- Intersectsim

      positive <- c(as.matrix(Intersect - listIntersect[[i]]))
      positive <- positive[!is.na(positive)]
      listnegIntersect[i] <- all(is.real.positive(positive))

      writeLines(paste(Sys.time(), paste("Iteration", i)))
      writeLines(" ")
      writeLines("Intersect metric")
      print(Intersectsim)
      writeLines("------------------------------------")
      if (i == nsim) {
        writeLines("ANALYSIS COMPLETED")
        writeLines("------------------------------------")
        writeLines(" ")
      }
    }

    # Calculate p-values
    Intersectsim_pval <- Matrixsim

    vector <- c()
    for (c in 1:length(data)) {
      for (r in 1:length(data)) {
        if (c <= r) next
        for (i in 1:nsim) {
          vector[i] <- listIntersect[[i]][r, c]
          Intersectsim_pval[r, c] <- length(which(vector <= Intersect[r, c])) / nsim
        }
      }
    }

    vector <- c()

    Intersect_together_pval <- length(which(listnegIntersect == TRUE)) / nsim
  }


  if (test == TRUE) {
    list <- list()
    list[[1]] <- Intersect
    list[[2]] <- Intersectsim_pval
    list[[3]] <- Intersect_together_pval
    list[[4]] <- listIntersect

    names(list) <- c("Intersection_metric", "Intersection_metric_p_values", "Intersection_metric_p_values_combined", "Intersection_metric_simulations")
    print(list[1])
    return(list)
  } else {
    return(Intersect)
  }
}
