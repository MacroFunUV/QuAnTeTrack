#' Similarity metric using Fréchet distance
#'
#' \code{simil_Frechet_metric()} computes similarity metrics between two or more trajectories using
#' the Fréchet distance. It allows for different superposition methods
#' to align trajectories before calculating the Fréchet distance  metrics. The function also supports
#' testing with simulations to calculate *p*-values for the Fréchet distance metrics.
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param test Logical; if \code{TRUE}, the function compares the observed Fréchet distances against
#' simulated trajectories and calculates *p*-values. Default is \code{FALSE}.
#' @param sim A \code{track simulation} R object consisting of a list of simulated trajectories to use for comparison when \code{test = TRUE}.
#' @param superposition A character string indicating the method used to align trajectories.
#' Options are \code{"None"}, \code{"Centroid"}, or \code{"Origin"}. Default is \code{"None"}.
#'
#' @details
#' The \code{simil_Frechet_metric()} function calculates the similarity between trajectories using
#' the \code{Frechet()} function from the \pkg{SimilarityMeasures} package.
#'
#' The Fréchet distance is a measure of similarity between two curves or continuous trajectories, which takes into account
#' both the order and location of points within the trajectories (Besse et al. 2015). The distance can be described by the
#' analogy of a person walking a dog on an extendable leash (Aronov et al. 2006). Both the person and the dog move along
#' their respective trajectories, with each able to adjust their speed but not retrace their steps. The Fréchet distance is
#' the minimum leash length required to keep the dog connected to the person throughout the walk (Cleasby et al., 2019).
#'
#' Unlike other trajectory comparison techniques, such as Dynamic Time Warping, the Fréchet distance focuses on the overall
#' shape of the trajectories rather than matching specific points. As a result, it is sensitive to noise because all points
#' of the compared trajectories are considered in its calculation. However, it can still be a powerful tool for trajectory
#' clustering and comparison, particularly when shape is the primary concern (Cleasby et al., 2019).
#'
#' Note that when comparing real trajectories that are very disparate or those simulated under an unconstrained method,
#' the resulting trajectories may not be suitable for Fréchet distance calculations. In such cases, the Fréchet distance is
#' returned as -1 to indicate an invalid measurement.
#'
#' The function offers three different superposition methods to align the trajectories
#' before \code{Frechet()} is applied:
#' \itemize{
#'   \item \code{"None"}: No superposition is applied.
#'   \item \code{"Centroid"}: Trajectories are shifted to align based on their centroids.
#'   \item \code{"Origin"}: Trajectories are shifted to align based on their starting point.
#' }
#'
#' If \code{test = TRUE}, the function can compute *p*-values by comparing the observed Fréchet
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
#' \item{Frechet_distance_metric}{A numeric matrix of pairwise Fréchet distances
#' between trajectories. Each entry represents the Fréchet distance between the corresponding pair of trajectories.}
#'
#' \item{Frechet_distance_metric_p_values}{(If \code{test = TRUE}) A numeric matrix of raw pairwise *p*-values,
#' computed by Monte Carlo tail tests with the (+1) correction
#' (Phipson & Smyth, 2010):
#' \deqn{p = (1 + \#\{\text{extreme}\}) / (nsim + 1)}.
#' Each entry reflects the probability of observing a Fréchet distance as extreme as the observed one,
#' given the null hypothesis of no difference.}
#'
#' \item{Frechet_distance_metric_p_values_BH}{(If \code{test = TRUE}) A numeric matrix of Benjamini–Hochberg (BH)
#' adjusted *p*-values controlling the false discovery rate (FDR), applied across the set of unique pairwise tests
#' (Benjamini & Hochberg, 1995).}
#'
#' \item{Frechet_metric_p_values_combined}{(If \code{test = TRUE}) A single numeric value representing
#' the combined *p*-value across all Fréchet distances (based on a global dominance criterion, evaluating in how many simulations the observed distances are smaller than the simulated ones across all trajectory pairs simultaneously).
#' This indicates the overall significance of the observed Fréchet distances relative to simulations.}
#'
#' \item{Frechet_distance_metric_simulations}{(If \code{test = TRUE}) A list containing matrices of Fréchet distances
#' for each simulation iteration, allowing for inspection of the distribution of Fréchet distances across randomized scenarios.}
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
#' # Example 1: Simulating tracks using the "Directed" model and comparing Frechet distance
#' # in the PaluxyRiver dataset
#' s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
#' simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#'
#' # Example 2: Simulating tracks using the "Constrained" model and comparing Frechet distance
#' # in the PaluxyRiver dataset  using the "Centroid" superposition method
#' s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
#' simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2, superposition = "Centroid")
#'
#' # Example 3: Simulating tracks using the "Unconstrained" model and comparing Frechet distance
#' # in the PaluxyRiver dataset using the "Origin" superposition method
#' s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
#' simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3, superposition = "Origin")
#'
#'
#' @importFrom SimilarityMeasures Frechet
#' @importFrom dplyr bind_rows
#' @importFrom stats na.omit
#' @importFrom schoolmath is.real.positive
#'
#' @seealso \code{\link{tps_to_track}}, \code{\link{simulate_track}}, \code{\link[SimilarityMeasures]{Frechet}}
#'
#' @export

simil_Frechet_metric <- function(data, test = FALSE, sim = NULL, superposition = "None") {

  ## ---- Defaults & checks ----
  if (is.null(test)) test <- FALSE
  if (is.null(superposition)) superposition <- "None"

  # data structure checks
  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'track' R object, which is a list consisting of two elements.")
  }
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("The two elements of 'data' must be lists.")
  }
  if (!is.logical(test)) {
    stop("'test' argument should be TRUE or FALSE.")
  }
  if (test && is.null(sim)) {
    stop("A 'sim' argument must be provided when 'test' is TRUE.")
  }

  # superposition option checks
  valid_superpositions <- c("None","Centroid","Origin")
  if (!superposition %in% valid_superpositions) {
    stop("Invalid 'superposition' argument. One of 'None', 'Centroid', or 'Origin' must be chosen.")
  }

  # simulation structure checks
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
      trks[[i]][, 1] <- trks[[i]][, 1] - mean(trks[[i]][, 1], na.rm = TRUE)
      trks[[i]][, 2] <- trks[[i]][, 2] - mean(trks[[i]][, 2], na.rm = TRUE)
    }
  } else if (superposition == "Origin") {
    for (i in seq_along(trks)) {
      trks[[i]][, 1] <- trks[[i]][, 1] - trks[[i]][1, 1]
      trks[[i]][, 2] <- trks[[i]][, 2] - trks[[i]][1, 2]
    }
  }
  # "None": no translation

  ## ---- Output matrix scaffold ----
  Matrixsim <- data.frame(matrix(nrow = length(trks), ncol = length(trks)))
  colnames(Matrixsim) <- names(trks)
  rownames(Matrixsim) <- names(trks)
  Frechet_metric <- Matrixsim  # avoid masking the function name

  ## ---- Observed pairwise Fréchet ----
  for (c in seq_along(trks)) {
    for (r in seq_along(trks)) {
      if (c <= r) next
      val <- suppressWarnings(
        SimilarityMeasures::Frechet(
          as.matrix(trks[[r]][, 1:2, drop = FALSE]),
          as.matrix(trks[[c]][, 1:2, drop = FALSE])
        )
      )
      Frechet_metric[r, c] <- val
      Frechet_metric[c, r] <- val  # make symmetric (minimal change)
    }
  }
  diag(Frechet_metric) <- NA  # explicit NA diagonal (minimal change)

  if (test) {
    ## ---- Apply global superposition to simulations ----
    if (superposition %in% c("Centroid","Origin")) {
      for (i in seq_along(sim)) {
        for (j in seq_along(sim[[i]])) {
          if (superposition == "Centroid") {
            sim[[i]][[j]][, 1] <- sim[[i]][[j]][, 1] - mean(sim[[i]][[j]][, 1], na.rm = TRUE)
            sim[[i]][[j]][, 2] <- sim[[i]][[j]][, 2] - mean(sim[[i]][[j]][, 2], na.rm = TRUE)
          } else { # "Origin"
            sim[[i]][[j]][, 1] <- sim[[i]][[j]][, 1] - sim[[i]][[j]][1, 1]
            sim[[i]][[j]][, 2] <- sim[[i]][[j]][, 2] - sim[[i]][[j]][1, 2]
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
    listFrechet <- vector("list", nsim)
    listnegFrechet <- logical(nsim)

    for (i in seq_len(nsim)) {

      ## CHANGE: reset Frechet_sim each iteration (robustness)
      Frechet_sim <- Matrixsim

      sim[[i]]$Trajectory <- as.factor(sim[[i]]$Trajectory)
      levs <- levels(sim[[i]]$Trajectory)

      for (c in seq_along(trks)) {
        for (r in seq_along(trks)) {
          if (c <= r) next
          A <- sim[[i]][sim[[i]]$Trajectory == levs[r], 1:2, drop = FALSE]
          B <- sim[[i]][sim[[i]]$Trajectory == levs[c], 1:2, drop = FALSE]
          Frechet_sim[r, c] <- suppressWarnings(
            SimilarityMeasures::Frechet(as.matrix(A), as.matrix(B))
          )

          ## (optional but consistent): make simulated matrix symmetric
          Frechet_sim[c, r] <- Frechet_sim[r, c]
        }
      }
      diag(Frechet_sim) <- NA

      listFrechet[[i]] <- Frechet_sim
      diff_vec <- c(as.matrix(Frechet_metric - Frechet_sim))
      diff_vec <- diff_vec[!is.na(diff_vec)]
      listnegFrechet[i] <- if (length(diff_vec)) all(schoolmath::is.real.positive(diff_vec)) else FALSE

      ## ---- Progress messages (console) ----
      message(paste(Sys.time(), paste("Iteration", i)))
      message(" ")
      message("Frechet metric")
      print(Frechet_sim)
      message("------------------------------------")
      if (i == nsim) {
        message("ANALYSIS COMPLETED")
        message("------------------------------------")
        message(" ")
      }
    }

    ## ---- p-values (+1 correction to avoid zeros) ----
    Frechet_pval <- Matrixsim
    for (c in seq_along(trks)) {
      for (r in seq_along(trks)) {
        if (c <= r) next
        vec <- numeric(nsim)
        for (i in seq_len(nsim)) vec[i] <- as.matrix(listFrechet[[i]])[r, c]
        p_rc <- (1 + sum(vec <= Frechet_metric[r, c], na.rm = TRUE)) / (nsim + 1)
        Frechet_pval[r, c] <- p_rc
        Frechet_pval[c, r] <- p_rc
      }
    }
    diag(Frechet_pval) <- NA

    # Combined p-value kept as original logic (no change requested)
    Frechet_together_pval <- sum(listnegFrechet) / nsim

    ## ---- Benjamini–Hochberg (BH) correction ----
    tmp <- as.matrix(Frechet_pval)
    vals <- as.vector(tmp)
    keep <- !is.na(vals)
    vals_adj <- vals
    vals_adj[keep] <- p.adjust(vals[keep], method = "BH")
    Frechet_pval_BH <- as.data.frame(matrix(vals_adj,
                                            nrow = nrow(tmp),
                                            ncol = ncol(tmp),
                                            byrow = FALSE))
    colnames(Frechet_pval_BH) <- colnames(Frechet_pval)
    rownames(Frechet_pval_BH) <- rownames(Frechet_pval)
    diag(Frechet_pval_BH) <- NA

    return(list(
      Frechet_distance_metric = Frechet_metric,
      Frechet_distance_metric_p_values = Frechet_pval,
      Frechet_distance_metric_p_values_BH = Frechet_pval_BH,
      Frechet_metric_p_values_combined = Frechet_together_pval,
      Frechet_distance_metric_simulations = listFrechet
    ))
  }

  ## ---- Return ----
  list(Frechet_distance_metric = Frechet_metric)
}
