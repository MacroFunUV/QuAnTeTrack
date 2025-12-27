#' Similarity metric using Dynamic Time Warping (DTW)
#'
#' \code{simil_DTW_metric()}  computes similarity metrics between two or more trajectories using
#' Dynamic Time Warping (DTW). It allows for different superposition methods
#' to align trajectories before calculating the DTW metric. The function also supports
#' testing with simulations to calculate *p*-values for the DTW distance metrics.
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param test Logical; if \code{TRUE}, the function compares the observed DTW distances against
#' simulated trajectories and calculates *p*-values. Default is \code{FALSE}.
#' @param sim A \code{track simulation} R object consisting of a list of simulated trajectories to use for comparison when \code{test = TRUE}.
#' @param superposition A character string indicating the method used to align trajectories.
#' Options are \code{"None"}, \code{"Centroid"}, or \code{"Origin"}. Default is \code{"None"}.
#'
#' @details
#' The \code{simil_DTW_metric()} function calculates the similarity between trajectories using
#' the Dynamic Time Warping (DTW) algorithm from the \pkg{dtw} package. The \code{dtw()} function
#' is used with the \code{dist.method} argument set to \code{"Euclidean"} for computing the local distances
#' between points in the trajectories.
#'
#' DTW aligns two time series by minimizing the cumulative distance between their points, creating an optimal alignment despite
#' variations in length or temporal distortions. The algorithm constructs a distance matrix where each element represents the cost of aligning
#' points between the two series and finds a warping path through this matrix that minimizes the total distance. The warping path is contiguous and monotonic,
#' starting from the bottom-left corner and ending at the top-right corner (Cleasby et al., 2019).
#'
#' DTW measures are non-negative and unbounded, with larger values indicating greater dissimilarity between the time series. This method has been used in various contexts, including ecological studies to analyze and cluster trajectory data (Cleasby et al., 2019).
#'
#' Potential limitations and biases of DTW include sensitivity to noise and outliers, computational complexity, and the need for appropriate distance metrics.
#' Additionally, DTW may not always account for all structural differences between trajectories and can be biased by the chosen alignment constraints.
#' While DTW can handle trajectories of different lengths due to its elastic nature, having trajectories of similar lengths can
#' improve the accuracy and interpretability of the similarity measure. Similar lengths result in a more meaningful alignment and
#' can make the computation more efficient. When trajectories differ significantly in length, preprocessing or normalization might
#' be necessary, and careful analysis is required to understand the alignment path. The function’s flexibility in handling different
#' lengths allows it to be applied in various contexts. However, large differences in trajectory lengths might introduce potential biases that should be
#' considered when interpreting the results.
#'
#' The function offers three different superposition methods to align the trajectories
#' before \code{DTW()} is applied:
#' \itemize{
#'   \item \code{"None"}: No superposition is applied.
#'   \item \code{"Centroid"}: Trajectories are shifted to align based on their centroids.
#'   \item \code{"Origin"}: Trajectories are shifted to align based on their starting point.
#' }
#'
#' If \code{test = TRUE}, the function can compute *p*-values by comparing the observed DTW
#' distances with those generated from a set of simulated trajectories. The *p*-values
#' are calculated for both individual trajectory pairs and for the entire set of trajectories.
#' Pairwise *p*-values are computed with a Monte Carlo tail test using the (+1) correction to avoid
#' zero-values (see Phipson & Smyth, 2010):
#' \deqn{p = (1 + \#\{\text{extreme}\}) / (nsim + 1)}
#'
#' These raw *p*-values are then adjusted for multiple comparisons across the set
#' of unique pairwise tests using the Benjamini–Hochberg (BH) procedure for
#' false discovery rate (FDR) control (Benjamini & Hochberg, 1995). Both the raw and the BH-adjusted *p*-value
#' matrices are returned in the output object, allowing users to inspect either
#' uncorrected or corrected results. In addition, a global combined *p*-value is
#' provided, summarizing the overall deviation from the null across all pairs.
#'
#' @return A \code{track similarity} R object consisting of a list containing the following elements:
#'
#' \item{DTW_distance_metric}{A numeric matrix of pairwise Dynamic Time Warping (DTW) distances
#' between trajectories. Each entry represents the DTW distance between the corresponding pair of trajectories.}
#'
#' \item{DTW_distance_metric_p_values}{(If \code{test = TRUE}) A numeric matrix of raw pairwise *p*-values,
#' computed by Monte Carlo tail tests with the (+1) correction
#' (Phipson & Smyth, 2010):
#' \deqn{p = (1 + \#\{\text{extreme}\}) / (nsim + 1)}.
#' Each entry reflects the probability of observing a DTW distance as extreme as the observed one,
#' given the null hypothesis of no difference.}
#'
#' \item{DTW_distance_metric_p_values_BH}{(If \code{test = TRUE}) A numeric matrix of Benjamini–Hochberg (BH)
#' adjusted *p*-values controlling the false discovery rate (FDR), applied across the set of unique pairwise tests
#' (Benjamini & Hochberg, 1995).}
#'
#' \item{DTW_metric_p_values_combined}{(If \code{test = TRUE}) A single numeric value representing
#' the combined *p*-value across all DTW distances (based on the global statistic: the sum of pairwise distances).
#' This indicates the overall significance of the observed DTW distances relative to simulations.}
#'
#' \item{DTW_distance_metric_simulations}{(If \code{test = TRUE}) A list containing matrices of DTW distances
#' for each simulation iteration, allowing for inspection of the distribution of DTW distances across randomized scenarios.}
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
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing.
#' Journal of the Royal statistical society: series B (Methodological), 57(1), 289-300.
#'
#' Cleasby, I. R., Wakefield, E. D., Morrissey, B. J., Bodey, T. W., Votier, S. C., Bearhop, S., & Hamer, K. C. (2019). Using time-series similarity measures to compare animal movement trajectories in ecology. Behavioral Ecology and Sociobiology, 73, 1-19.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn.
#' Statistical applications in genetics and molecular biology, 9(1).
#'
#' @examples
#' # Example 1: Simulating tracks using the "Directed" model and comparing DTW distance
#' # in the PaluxyRiver dataset
#' s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
#' simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#'
#' # Example 2: Simulating tracks using the "Constrained" model and comparing DTW distance
#' # in the PaluxyRiver dataset
#' s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
#' simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2, superposition = "None")
#'
#' # Example 3: Simulating tracks using the "Unconstrained" model and comparing DTW distance
#' # in the PaluxyRiver dataset
#' s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
#' simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3, superposition = "None")
#'
#' # Example 4: Simulating and comparing DTW distance in the MountTom dataset using the
#' # "Centroid" superposition method
#' sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
#' s4 <- simulate_track(sbMountTom, nsim = 3)
#' simil_DTW_metric(sbMountTom, test = TRUE, sim = s4, superposition = "Centroid")
#'
#' # Example 5: Simulating and comparing DTW distance in the MountTom dataset using the
#' # "Origin" superposition method
#' sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
#' s5 <- simulate_track(sbMountTom, nsim = 3)
#' simil_DTW_metric(sbMountTom, test = TRUE, sim = s5, superposition = "Origin")
#'
#' @importFrom dtw dtw
#' @importFrom dplyr bind_rows
#' @importFrom stats na.omit
#' @importFrom schoolmath is.real.positive
#'
#' @seealso \code{\link{tps_to_track}}, \code{\link{simulate_track}}, \code{\link[dtw]{dtw}}
#'
#' @export


simil_DTW_metric <- function(data, test = FALSE, sim = NULL, superposition = "None") {

  ## ---- Defaults & checks ----
  if (is.null(test)) test <- FALSE
  if (is.null(superposition)) superposition <- "None"

  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'track' R object, which is a list consisting of two elements.")
  }
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("The two elements of 'data' must be lists.")
  }
  if (!is.logical(test)) stop("'test' argument should be TRUE or FALSE.")
  if (test && is.null(sim)) stop("A 'sim' argument must be provided when 'test' is TRUE.")

  valid_superpositions <- c("None","Centroid","Origin")
  if (!superposition %in% valid_superpositions) {
    stop("Invalid 'superposition' argument. One of 'None', 'Centroid', or 'Origin' must be chosen.")
  }

  if (!is.null(sim)) {
    if (!is.list(sim)) stop("The 'sim' argument must be a list.")
    if (length(sim[[1]]) != length(data[[1]])) {
      stop("The 'sim' list must have the same number of trajectories as 'data'.")
    }
  }

  ## ---- Extract trajectories ----
  trks <- data[[1]]

  ## ---- Global per-trajectory superposition ----
  if (superposition == "Centroid") {
    for (i in seq_along(trks)) {
      trks[[i]][,1] <- trks[[i]][,1] - mean(trks[[i]][,1], na.rm = TRUE)
      trks[[i]][,2] <- trks[[i]][,2] - mean(trks[[i]][,2], na.rm = TRUE)
    }
  } else if (superposition == "Origin") {
    for (i in seq_along(trks)) {
      trks[[i]][,1] <- trks[[i]][,1] - trks[[i]][1,1]
      trks[[i]][,2] <- trks[[i]][,2] - trks[[i]][1,2]
    }
  }

  ## ---- Output matrix ----
  Matrixsim <- data.frame(matrix(nrow = length(trks), ncol = length(trks)))
  colnames(Matrixsim) <- names(trks)
  rownames(Matrixsim) <- names(trks)
  DTW <- Matrixsim

  ## ---- Observed pairwise DTW ----
  for (c in seq_along(trks)) {
    for (r in seq_along(trks)) {
      if (c <= r) next
      val <- dtw::dtw(
        as.matrix(trks[[r]][,1:2, drop = FALSE]),
        as.matrix(trks[[c]][,1:2, drop = FALSE]),
        distance.only = TRUE
      )$distance
      DTW[r, c] <- val
      DTW[c, r] <- val  # make symmetric (minimal change)
    }
  }
  diag(DTW) <- NA  # explicit NA diagonal (minimal change)

  if (test) {
    ## ---- Apply global superposition to simulations ----
    if (superposition %in% c("Centroid","Origin")) {
      for (i in seq_along(sim)) {
        for (j in seq_along(sim[[i]])) {
          if (superposition == "Centroid") {
            sim[[i]][[j]][,1] <- sim[[i]][[j]][,1] - mean(sim[[i]][[j]][,1], na.rm = TRUE)
            sim[[i]][[j]][,2] <- sim[[i]][[j]][,2] - mean(sim[[i]][[j]][,2], na.rm = TRUE)
          } else {
            sim[[i]][[j]][,1] <- sim[[i]][[j]][,1] - sim[[i]][[j]][1,1]
            sim[[i]][[j]][,2] <- sim[[i]][[j]][,2] - sim[[i]][[j]][1,2]
          }
        }
      }
    }

    ## ---- Collapse simulations ----
    listSIM <- vector("list", length(sim))
    for (i in seq_along(sim)) listSIM[[i]] <- dplyr::bind_rows(sim[[i]])
    sim <- listSIM

    ## ---- Simulated metrics ----
    nsim <- length(sim)
    DTWsim <- Matrixsim
    listDTW <- vector("list", nsim)
    listnegDTW <- logical(nsim)

    for (i in seq_len(nsim)) {
      sim[[i]]$Trajectory <- as.factor(sim[[i]]$Trajectory)
      levs <- levels(sim[[i]]$Trajectory)

      for (c in seq_along(trks)) {
        for (r in seq_along(trks)) {
          if (c <= r) next
          A <- sim[[i]][sim[[i]]$Trajectory == levs[r], 1:2, drop = FALSE]
          B <- sim[[i]][sim[[i]]$Trajectory == levs[c], 1:2, drop = FALSE]
          DTWsim[r, c] <- suppressWarnings(
            dtw::dtw(as.matrix(A), as.matrix(B), distance.only = TRUE)$distance
          )
        }
      }

      listDTW[[i]] <- DTWsim
      diff_vec <- c(as.matrix(DTW - DTWsim))
      diff_vec <- diff_vec[!is.na(diff_vec)]
      listnegDTW[i] <- if (length(diff_vec)) all(schoolmath::is.real.positive(diff_vec)) else FALSE

      ## ---- Progress messages ----
      message(paste(Sys.time(), paste("Iteration", i)))
      message(" ")
      message("DTW metric")
      print(DTWsim)
      message("------------------------------------")
      if (i == nsim) {
        message("ANALYSIS COMPLETED")
        message("------------------------------------")
        message(" ")
      }
    }

    ## ---- p-values (with +1 correction) ----
    DTWsim_pval <- Matrixsim
    for (c in seq_along(trks)) {
      for (r in seq_along(trks)) {
        if (c <= r) next
        vec <- numeric(nsim)
        for (i in seq_len(nsim)) vec[i] <- listDTW[[i]][r, c]
        # Monte Carlo lower-tail with (+1) correction (minimal change)
        DTWsim_pval[r, c] <- (1 + sum(vec <= DTW[r, c], na.rm = TRUE)) / (nsim + 1)
      }
    }
    diag(DTWsim_pval) <- NA  # NA diagonal for p-values (consistency)

    # Combined (kept as in original logic)
    DTW_together_pval <- sum(listnegDTW) / nsim

    ## ---- BH correction (minimal, robust to data.frame) ----
    tmp <- as.matrix(DTWsim_pval)
    vals <- as.vector(tmp)
    keep <- !is.na(vals)
    vals_adj <- vals
    vals_adj[keep] <- p.adjust(vals[keep], method = "BH")
    DTWsim_pval_BH <- as.data.frame(matrix(vals_adj,
                                           nrow = nrow(tmp),
                                           ncol = ncol(tmp),
                                           byrow = FALSE))
    colnames(DTWsim_pval_BH) <- colnames(DTWsim_pval)
    rownames(DTWsim_pval_BH) <- rownames(DTWsim_pval)
    diag(DTWsim_pval_BH) <- NA

    return(list(
      DTW_distance_metric = DTW,
      DTW_distance_metric_p_values = DTWsim_pval,
      DTW_distance_metric_p_values_BH = DTWsim_pval_BH,
      DTW_metric_p_values_combined = DTW_together_pval,
      DTW_distance_metric_simulations = listDTW
    ))
  }

  ## ---- Return ----
  list(DTW_distance_metric = DTW)
}
