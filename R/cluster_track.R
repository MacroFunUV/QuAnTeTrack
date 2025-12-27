#' Cluster tracks based on movement parameters
#'
#' \code{cluster_track()} clusters trajectories based on various movement and velocity parameters calculated for each track.
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param veltrack A \code{track velocity} R object consisting of a list of lists, where each sublist contains the computed parameters for a corresponding track.
#' @param variables A character vector specifying the movement parameters to be used in the clustering analysis. Valid parameter names include:
#'   \code{"TurnAng"}, \code{"sdTurnAng"}, \code{"Distance"}, \code{"Length"}, \code{"StLength"}, \code{"sdStLength"},
#'   \code{"Sinuosity"}, \code{"Straightness"}, \code{"Velocity"}, \code{"sdVelocity"}, \code{"MaxVelocity"}, \code{"MinVelocity"},
#'   \code{"TrackWidth"}, \code{"PaceAng"}.
#' @param analysis Character; clustering backend. One of \code{"hclust"} (default; distance-based hierarchical clustering)
#'   or \code{"mclust"} (model-based Gaussian mixtures).
#' @param k Integer or \code{NULL}; number of clusters to cut the dendrogram when \code{analysis = "hclust"}.
#'   Default \code{2}. If \code{NULL}, no cut is made and \code{classification} is returned as \code{NULL}.
#' @param dist_method Character; distance used when \code{analysis = "hclust"}.
#'   One of \code{"euclidean"} (default), \code{"maximum"}, \code{"manhattan"}, \code{"canberra"},
#'   \code{"binary"}, or \code{"minkowski"}.
#' @param hclust_method Character; agglomeration method when \code{analysis = "hclust"}.
#'   One of  \code{"complete"} (default), \code{"ward.D"}, \code{"ward.D2"}, \code{"single"},
#'   \code{"average"} (= UPGMA), \code{"mcquitty"} (= WPGMA), \code{"median"} (= WPGMC), or \code{"centroid"} (= UPGMC).
#' @param transform Logical; if \code{TRUE} (default), applies monotone transformations suited to variable support
#'   (e.g., log/logit) before clustering.
#' @param scale Logical; if \code{TRUE} (default), z-transforms (centers and scales) the variables used in clustering.
#' @param max_clusters Integer. Maximum number of clusters (mixture components) to consider when \code{analysis = "mclust"}.
#'   If \code{NULL} (default), uses \code{min(9, n)}, where \code{n} is the number of valid tracks (observations) entering the model.
#'   The minimum number of clusters considered by default is \code{1}.
#'
#' @details
#' The \code{cluster_track()} function can perform distance-based hierarchical clustering via the \code{hclust()} function from the \pkg{stats} package (\code{analysis = "hclust"}, default) or model-based clustering via the \code{Mclust()} function from the \pkg{mclust} package
#' (\code{analysis = "mclust"}).
#'
#' The function first filters out tracks with fewer than four steps, as these tracks may not provide reliable movement data.
#' It then calculates various movement parameters for each remaining track, including turning angles, distances, lengths, sinuosity, straightness, velocities,
#' and additional ichnological descriptors (trackway width, pace angulation).
#' Finally, the selected movement parameters are used as input for clustering the tracks.
#'
#' @details
#' When \code{analysis = "hclust"}: the function performs agglomerative hierarchical clustering
#' on a pairwise dissimilarity matrix computed from the working feature matrix (which already reflects the
#' internal expansion of circular variables to sine/cosine components, and any requested transformations and
#' z-scaling).
#'
#' The algorithm starts with each track as its own cluster and iteratively merges the two closest clusters
#' according to the chosen \code{hclust_method} (linkage criterion) until a single cluster remains. The full
#' merge history is represented by a dendrogram that can be inspected and cut at any desired number of clusters.
#'
#' Advantages include being distribution-free, since the method does not assume Gaussianity or any parametric mixture form.
#' It is also flexible, supporting multiple distance metrics, and highly interpretable, as the dendrogram reveals
#' multiscale clustering structure. In addition, the number of clusters does not need to be fixed in advance,
#' since it can be chosen a posteriori by cutting the tree. Disadvantages are that
#' the greedy nature of the algorithm makes merges irreversible, so no reassignment
#' is possible once clusters are combined. Results can also be sensitive to variable scaling, as well as the
#' choice of distance metric and linkage method. Certain linkage strategies, such as \code{"median"} or
#' \code{"centroid"}, may produce dendrogram inversions, and the \code{"ward.*"} methods require Euclidean
#' geometry to preserve their variance-minimization interpretation.

#'
#' If \code{transform = TRUE} and/or \code{scale = TRUE}, the distance is computed on the
#' transformed/scaled space. This is generally recommended so that variables on different scales do not dominate the
#' distances. Missing values in the selected variables will cause \code{stats::dist()} to fail; ensure inputs are
#' complete (or impute externally) before clustering.
#'
#' The argument \code{k} controls the number of clusters to cut the dendrogram. By default it is set to \code{2}.
#' If \code{NULL}, no cut is made and \code{classification} is returned as \code{NULL}. Setting \code{k = 1}
#' yields a single cluster containing all tracks.
#'
#' The argument \code{dist_method} specifies the distance metric used to build the dissimilarity matrix,
#' passed to \code{stats::dist}. For a detailed explanation of distance metrics see \code{?dist}.
#'
#' The argument \code{hclust_method} specifies the agglomeration (linkage) method used for hierarchical clustering,
#' passed to \code{stats::hclust}. For a detailed explanation of linkage procedures see \code{?hclust}.
#'
#' When \code{analysis = "mclust"}, clustering is performed using finite Gaussian mixture models
#' via the \code{Mclust()} function from the \pkg{mclust} package.
#' This approach assumes that the data are generated from a mixture of multivariate Gaussian distributions,
#' and selects both the number of mixture components and the covariance parameterization by optimizing
#' the Bayesian Information Criterion (BIC).
#'
#' Advantages include statistical rigor, automatic model selection, and the ability to capture ellipsoidal
#' cluster shapes. Limitations are the reliance on distributional assumptions (normality, homoscedasticity
#' depending on the model), potential sensitivity to outliers, and increased computational cost for
#' high-dimensional data. These assumptions are partially mitigated by the variable transformations
#' applied before clustering (e.g., log or logit transforms), which help approximate Gaussian-like
#' distributions.
#'
#' If only one parameter is selected, the clustering is performed using equal variance ("E") and variable
#' variance ("V") Gaussian models. If more than one parameter is selected, all covariance models available
#' in \code{mclust.options("emModelNames")} are considered.
#'
#' The argument \code{max_clusters} controls the maximum number of clusters (mixture components) to be explored.
#' If \code{NULL} (default), the maximum is set to \code{min(9, n)}, where \code{n} is the number of valid tracks
#' (observations) entering the model. The minimum number of clusters considered by default is \code{1}.
#'
#' The following movement parameters can be included in the clustering:
#'
#'   * \code{"TurnAng"}: Turning angles for the track, measured in degrees. Internally expanded to
#'     sine and cosine components (circular treatment) before analysis. This measures how much the
#'     direction of movement changes at each step.
#'
#'   * \code{"sdTurnAng"}: The circular standard deviation of the turning angles, in degrees.
#'     Internally expanded to sine and cosine components before analysis to respect angular geometry.
#'
#'   * \code{"Distance"}: The total distance covered by the track, calculated as the sum of the
#'     straight-line distances between consecutive points (meters).
#'
#'   * \code{"Length"}: The overall length of the track, a straight-line distance between the
#'     starting and ending points (meters).
#'
#'   * \code{"StLength"}: Step lengths for each step of the track, representing how far the object
#'     moved between two consecutive points (meters).
#'
#'   * \code{"sdStLength"}: The standard deviation of the step lengths, showing how consistent the
#'     steps are in length (meters).
#'
#'   * \code{"Sinuosity"}: A measure of the track's winding nature, calculated as the ratio of the
#'     actual track length to the straight-line distance (dimensionless).
#'
#'   * \code{"Straightness"}: The straightness of the track, calculated as the straight-line
#'     distance divided by the total path length (dimensionless).
#'
#'   * \code{"Velocity"}: The average velocity of the track, calculated as the total distance divided
#'     by the time elapsed between the first and last footprint (meters per second).
#'
#'   * \code{"sdVelocity"}: The standard deviation of the velocity, indicating how much the velocity
#'     fluctuates throughout the track (meters per second).
#'
#'   * \code{"MaxVelocity"}: The maximum velocity achieved during the track, identifying the fastest
#'     point (meters per second).
#'
#'   * \code{"MinVelocity"}: The minimum velocity during the track, identifying the slowest point
#'     (meters per second).
#'
#'   * \code{"TrackWidth"}: Trackway width, i.e., the lateral spacing between left and right
#'     footprint series (units consistent with input coordinates, typically meters).
#'
#'   * \code{"PaceAng"}: Pace angulation (degrees), a classical ichnological descriptor summarizing
#'     the angular relation of successive steps within a trackway.
#'
#' The \code{cluster_track()} function has biological relevance in identifying groups of tracks with
#' similar movement parameters, providing insights into ecological and behavioral patterns. By
#' clustering tracks based on characteristics such as sinuosity, velocity, and turning angles (treated
#' circularly via sine/cosine), it helps detect movement patterns associated with specific behaviors.
#' This can reveal tracks potentially made by individuals moving together, useful for investigating
#' hypotheses on gregarious behavior, predation strategies, or coordinated movement. Additionally,
#' clustering serves as a preliminary step before similarity tests and simulations, refining track
#' selection and improving hypothesis testing in movement ecology studies.
#'
#' @return
#' A track clustering R object consisting of a list containing the following elements:
#'
#'   * \code{matrix}: A data frame containing the movement parameters calculated for each track
#'     (original scale/units as produced by \code{track_param()} and \code{velocity_track()}).
#'
#'   * \code{clust}:
#'       - When \code{analysis = "hclust"}: a list with:
#'           - \code{hclust}: The \code{stats::hclust} object (dendrogram and merge history).
#'           - \code{dist}: The \code{stats::dist} object used to fit \code{hclust}.
#'           - \code{k}: The requested number of clusters used to cut the tree (if provided).
#'           - \code{classification}: Cluster labels obtained by \code{stats::cutree(hclust, k)}; \code{NULL} if \code{k = NULL}.
#'
#'       - When \code{analysis = "mclust"}: an \code{Mclust} object containing the results of the
#'         model-based clustering analysis (optimal model selected by BIC). Main components include:
#'           - \code{call}: The matched call.
#'           - \code{data}: The input data matrix used by \code{Mclust}.
#'           - \code{modelName}: The covariance model at which the optimal BIC occurs.
#'           - \code{n}: Number of observations.
#'           - \code{d}: Data dimensionality.
#'           - \code{G}: Optimal number of mixture components.
#'           - \code{BIC}: Matrix of BIC values for models and components considered.
#'           - \code{loglik}: Log-likelihood for the selected model.
#'           - \code{df}: Number of estimated parameters.
#'           - \code{bic}: BIC value of the selected model.
#'           - \code{icl}: ICL value of the selected model.
#'           - \code{hypvol}: Hypervolume parameter for the noise component (if applicable), otherwise \code{NULL}.
#'           - \code{parameters}: List of fitted parameters, including:
#'               * \code{pro}: Mixing proportions (length \code{G}).
#'               * \code{mean}: Component means (vector or matrix with columns as components).
#'               * \code{variance}: List of variance parameters (structure depends on \code{modelName}).
#'           - \code{z}: Matrix of posterior probabilities; entry \code{(i, k)} is P(observation \code{i} in class \code{k}).
#'           - \code{classification}: Hard labels (MAP from \code{z}).
#'           - \code{uncertainty}: Classification uncertainty for each observation.
#'
#' @section Logo:
#' \if{html}{\figure{Logo.png}{options: width=30\%}}
#'
#' @author Humberto G. Ferrón
#' @author humberto.ferron@uv.es
#' @author Macroevolution and Functional Morphology Research Group (www.macrofun.es)
#' @author Cavanilles Institute of Biodiversity and Evolutionary Biology
#' @author Calle Catedrático José Beltrán Martínez, nº 2
#' @author 46980 Paterna - Valencia - Spain
#' @author Phone: +34 (9635) 44477
#'
#' @references
#' Alexander, R. M. (1976). Estimates of speeds of dinosaurs. Nature, 261(5556), 129-130.
#'
#' Ruiz, J., & Torices, A. (2013). Humans running at stadiums and beaches and the accuracy of speed estimations from fossil trackways. Ichnos, 20(1), 31-35.
#'
#'  Scrucca L., Fop M., Murphy T. B., & Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models. The R Journal, 8(1), 289-317.
#'
#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
#' @importFrom mclust mclust.options
#' @importFrom circular circular
#' @importFrom stats dist hclust cutree sd
#'
#' @seealso \code{\link{track_param}}, \code{\link{velocity_track}}, \code{\link[mclust]{Mclust}}, \code{\link[stats]{hclust}}
#'
#' @export
cluster_track <- function(
    data, veltrack, variables,
    analysis = c("hclust","mclust"),
    k = 2,
    dist_method = "euclidean",
    hclust_method = "complete",
    transform = TRUE, scale = TRUE,
    max_clusters = NULL
) {

  analysis <- match.arg(analysis)

  ## Errors and Warnings----

  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'track' R object, which is a list consisting of two elements.")
  }
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("The two elements of 'data' must be lists.")
  }
  if (length(variables) == 0) {
    stop("Error: No movement parameters specified for clustering. Please provide a valid vector of variables.")
  }
  if (length(data[[1]]) == 0) {
    stop("Error: No tracks available for clustering. Please check the input data.")
  }

  # Identify tracks with fewer than 4 steps
  short <- sapply(data[[1]], function(track) nrow(track) < 4)
  if (any(short)) {
    warning(paste(
      sum(short),
      "tracks were discarded for having fewer than 4 footprints. Discarded track indices:",
      paste(which(short), collapse = ", ")
    ))
  }

  ## Code----

  # Subset to valid tracks
  valid_indices <- which(!short)
  data <- subset_track(data, valid_indices)
  veltrack <- veltrack[valid_indices]

  # Compute parameters
  data_param <- track_param(data)

  # Create results matrix (original-scale summary, incl. new vars)
  matrix <- data.frame(matrix(nrow = length(data[[1]]), ncol = 14))
  colnames(matrix) <- c(
    "TurnAng", "sdTurnAng", "Distance", "Length", "StLength", "sdStLength",
    "Sinuosity", "Straightness", "Velocity", "sdVelocity",
    "MaxVelocity", "MinVelocity", "TrackWidth", "PaceAng"
  )
  rownames(matrix) <- paste0("Track ", valid_indices)

  for (i in 1:length(data[[1]])) {
    matrix[i, "TurnAng"]        <- data_param[[i]]$Mean_turning_angle               # degrees (circular mean)
    matrix[i, "sdTurnAng"]      <- data_param[[i]]$Standard_deviation_turning_angle # degrees (circular sd)
    matrix[i, "Distance"]       <- data_param[[i]]$Distance
    matrix[i, "Length"]         <- data_param[[i]]$Length
    matrix[i, "StLength"]       <- data_param[[i]]$Mean_step_length
    matrix[i, "sdStLength"]     <- data_param[[i]]$Standard_deviation_step_length
    matrix[i, "Sinuosity"]      <- data_param[[i]]$Sinuosity
    matrix[i, "Straightness"]   <- data_param[[i]]$Straightness
    matrix[i, "Velocity"]       <- veltrack[[i]]$Mean_velocity
    matrix[i, "sdVelocity"]     <- veltrack[[i]]$Standard_deviation_velocity
    matrix[i, "MaxVelocity"]    <- veltrack[[i]]$Maximum_velocity
    matrix[i, "MinVelocity"]    <- veltrack[[i]]$Minimum_velocity
    tw <- data_param[[i]]$TrackWidth
    pa <- data_param[[i]]$PaceAng
    matrix[i, "TrackWidth"] <- if (is.null(tw)) NA_real_ else tw
    matrix[i, "PaceAng"]    <- if (is.null(pa)) NA_real_ else pa
  }

  ## Build working matrix for clustering (expand circular vars to sin/cos)
  work <- data.frame(row.names = rownames(matrix))
  add_var <- function(vn) {
    if (vn == "TurnAng") {
      ang <- circular::circular(matrix$TurnAng, units = "degrees", modulo = "2pi")
      work$TurnAng_sin <<- sin(ang)
      work$TurnAng_cos <<- cos(ang)
    } else if (vn == "sdTurnAng") {
      sdang <- circular::circular(matrix$sdTurnAng, units = "degrees", modulo = "2pi")
      work$sdTurnAng_sin <<- sin(sdang)
      work$sdTurnAng_cos <<- cos(sdang)
    } else if (vn %in% colnames(matrix)) {
      work[[vn]] <<- matrix[[vn]]
    } else {
      stop(sprintf("Variable '%s' not recognized.", vn))
    }
  }
  invisible(lapply(variables, add_var))

  ## Transformations
  if (transform) {
    eps <- 1e-8
    logit <- function(x) log(x/(1 - x))
    transform_one <- function(x, name) {
      if (grepl("(_sin|_cos)$", name)) return(x)
      if (name %in% c("Straightness")) {
        x <- pmin(pmax(x, eps), 1 - eps); return(logit(x))
      }
      if (name %in% c("Sinuosity")) {
        return(log(pmax(x, 1 + eps)))
      }
      if (name %in% c("Distance","Length","StLength","sdStLength",
                      "Velocity","sdVelocity","MaxVelocity","MinVelocity",
                      "TrackWidth","PaceAng")) {
        shift <- if (any(x <= 0, na.rm = TRUE)) eps - min(x, na.rm = TRUE) + eps else 0
        return(log(x + shift))
      }
      x
    }
    for (nm in names(work)) work[[nm]] <- transform_one(work[[nm]], nm)
  }

  ## Z-scale
  if (scale) {
    for (nm in names(work)) {
      v <- work[[nm]]; m <- mean(v, na.rm = TRUE); s <- stats::sd(v, na.rm = TRUE)
      work[[nm]] <- if (is.na(s) || s == 0) v*0 else (v - m)/s
    }
  }

  ## Backend selection
  xmat  <- as.matrix(work)

  if (analysis == "mclust") {
    # Mclust with minimum clusters = 1 (as in mclust) and user-defined maximum
    n_obs <- nrow(xmat)
    if (is.null(max_clusters)) {
      max_clusters <- min(9L, n_obs)
    } else {
      max_clusters <- as.integer(max_clusters)
      if (is.na(max_clusters) || max_clusters < 1L) stop("'max_clusters' must be >= 1.")
      max_clusters <- min(max_clusters, n_obs)
    }
    Gmin <- 1L
    if (Gmin > max_clusters) Gmin <- max_clusters
    Gseq <- seq.int(from = Gmin, to = max_clusters, by = 1L)

    if (ncol(xmat) == 1) {
      clust <- Mclust(xmat, G = Gseq, modelNames = c("E", "V"))
    } else {
      clust <- Mclust(xmat, G = Gseq, modelNames = mclust.options("emModelNames"))
    }

  } else { # analysis == "hclust"
    if (nrow(xmat) < 2) {
      warning("Fewer than 2 valid tracks after filtering; returning without clustering.")
      clust <- list(hclust = NULL, dist = NULL, k = NULL, classification = NULL)
      return(list(matrix = matrix, clust = clust))
    }
    D  <- stats::dist(xmat, method = dist_method)
    hc <- stats::hclust(D, method = hclust_method)

    if (is.null(k)) {
      groups <- NULL
    } else {
      k <- as.integer(k)
      if (is.na(k) || k < 1) stop("'k' must be >= 1 or NULL.")
      k <- min(k, nrow(xmat))
      groups <- stats::cutree(hc, k = k)
    }

    clust <- list(hclust = hc, dist = D, k = if (is.null(k)) NULL else k, classification = groups)
  }

  return(list(matrix = matrix, clust = clust))
}
