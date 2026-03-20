#' Calculate combined probabilities of similarity or intersection metrics of tracks
#'
#' \code{combined_prob()} calculates the combined probabilities of similarity and intersection metrics
#' derived from different models. The function uses simulation data to extract *p*-values, providing insight into
#' the significance of combined metrics across various similarity assessments.
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param metrics A list of \code{track similarity} and/or \code{track intersection} R objects derived from different tests. All tests must be
#' based on the same number of simulations.
#' @param H1 Alternative hypothesis for intersection metrics. One of \code{"Higher"} or \code{"Lower"}.
#'
#' @details
#' The \code{combined_prob()} function combines evidence from multiple trajectory-based
#' metrics, including distance-based similarity measures (e.g., DTW and Fréchet)
#' and intersection metrics, using simulation-based hypothesis testing.
#'
#' For distance-based similarity metrics (such as DTW and Fréchet), smaller values
#' indicate greater similarity between trajectories. Accordingly, pairwise *p*-values
#' are computed using a lower-tail Monte Carlo test with the (+1) correction
#' (Phipson & Smyth, 2010):
#' \deqn{p = (1 + \#\{D_{sim} \le D_{obs}\}) / (nsim + 1)}
#' This tests whether the observed trajectories are at least as similar as expected
#' under the null model of random movement.
#'
#' For intersection metrics, the direction of the test depends on the alternative
#' hypothesis specified by \code{H1}. When \code{H1 = "Lower"}, the test evaluates
#' whether the observed number of intersections is significantly lower than expected
#' under random simulations (e.g., coordinated or gregarious movement). When
#' \code{H1 = "Higher"}, the test evaluates whether the observed number of intersections
#' is significantly higher than expected (e.g., chasing or predatory interactions).
#'
#' Combined pairwise *p*-values are obtained by jointly evaluating all metrics across
#' simulations, identifying those simulation replicates in which all metrics are
#' simultaneously as extreme as the observed values, given their respective directions.
#' A Monte Carlo test with the (+1) correction is used, and the resulting pairwise
#' *p*-values are adjusted for multiple comparisons using the Benjamini–Hochberg (BH)
#' procedure.
#'
#' In addition, a global combined *p*-value is computed based on a single global
#' statistic summarizing all pairwise comparisons, providing an overall assessment
#' of whether the observed set of trajectories departs from random expectations
#' across all metrics considered.
#'
#' @return A list containing:
#' \item{P_values}{A matrix of raw *p*-values for the combined metrics across all trajectories (pairwise).}
#' \item{P_values_BH}{A matrix of BH-adjusted *p*-values for the combined metrics (pairwise).}
#' \item{P_values_global}{A single numeric value: the overall probability of observing the combined metrics across all pairs of trajectories.}
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
#' # Example 1: "Directed" model and similarity metrics.
#' s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
#' DTW1 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#' Frechet1 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#' int1 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1,
#'   origin.permutation = "None")
#' combined_prob(PaluxyRiver, metrics = list(DTW1, Frechet1, int1), H1 = "Lower")
#'
#' # Example 2: "Constrained" model and similarity metrics.
#' s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
#' DTW2 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2,
#'   superposition = "None")
#' Frechet2 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2,
#'   superposition = "None")
#' int2 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2,
#'   origin.permutation = "Min.Box")
#' combined_prob(PaluxyRiver, metrics = list(DTW2, Frechet2, int2), H1 = "Lower")
#'
#' # Example 3: "Unconstrained" model and similarity metrics.
#' s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
#' DTW3 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3,
#'   superposition = "None")
#' Frechet3 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3,
#'   superposition = "None")
#' int3 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3,
#'   origin.permutation = "Conv.Hull")
#' combined_prob(PaluxyRiver, metrics = list(DTW3, Frechet3, int3), H1 = "Lower")
#'
#' @seealso \code{\link{tps_to_track}}, \code{\link{simulate_track}}, \code{\link{track_intersection}}, \code{\link{simil_DTW_metric}}, \code{\link{simil_Frechet_metric}}
#'
#' @export


combined_prob <- function(data, metrics = NULL, H1 = c("Higher", "Lower")) {

  ## ---- Checks ----
  if (!is.list(data) || length(data) < 2)
    stop("The 'data' argument must be a 'track' R object (list of two elements).")
  if (!is.list(data[[1]]) || !is.list(data[[2]]))
    stop("The two elements of 'data' must be lists.")
  if (is.null(metrics) || !is.list(metrics))
    stop("'metrics' must be a list of track similarity or intersection metric objects.")

  H1 <- match.arg(H1)
  is_lower <- identical(H1, "Lower")

  # detect which metrics are intersection metrics
  is_intersection_metric <- sapply(metrics, function(m) {
    nms <- names(m)
    any(c("Intersection_metric",
          "Intersection_metric_simulations",
          "Intersection_metric_p_values",
          "Intersection_metric_p_values_BH") %in% nms)
  })

  # Expect simulations in slot [[5]] for each metric object
  nsim_lengths <- sapply(metrics, function(m) length(m[[5]]))
  if (length(unique(nsim_lengths)) != 1)
    stop("All metric objects must have the same number of simulations (looked in [[5]]).")
  nsim <- nsim_lengths[1]

  trks   <- data[[1]]
  n      <- length(trks)
  pnames <- names(trks)

  make_blank <- function() matrix(NA_real_, n, n, dimnames = list(pnames, pnames))
  upper_idx  <- upper.tri(matrix(NA_real_, n, n), diag = FALSE)

  ## ---- Precompute diffs: diff = obs - sim (upper triangle only) ----
  # listc[[j]][[i]] stores a matrix of diffs for metric i, simulation j
  listc <- vector("list", nsim)
  for (j in seq_len(nsim)) {
    per_metric <- vector("list", length(metrics))
    for (i in seq_along(metrics)) {
      obs  <- as.matrix(metrics[[i]][[1]])
      simj <- as.matrix(metrics[[i]][[5]][[j]])

      diffs <- make_blank()
      diffs[upper_idx] <- obs[upper_idx] - simj[upper_idx]   # diff = obs - sim
      per_metric[[i]] <- diffs
    }
    listc[[j]] <- per_metric
  }

  ## ---- Helper: is simulation "as extreme as observed" for metric i? ----
  # Distances (DTW/Frechet): lower-tail => sim <= obs  <=>  obs - sim >= 0
  # Intersections:
  #   H1="Lower"  => sim <= obs  <=>  obs - sim >= 0
  #   H1="Higher" => sim >= obs  <=>  obs - sim <= 0
  is_extreme <- function(diff_value, is_intersection, is_lower) {
    if (!is.finite(diff_value)) return(FALSE)

    if (!is_intersection) {
      return(diff_value >= 0)             # distances: sim <= obs
    } else {
      if (is_lower) return(diff_value >= 0)  # intersections Lower: sim <= obs
      return(diff_value <= 0)                # intersections Higher: sim >= obs
    }
  }

  ## ---- Pairwise combined p-values (raw, +1 correction) ----
  P_pair <- make_blank()

  for (c in seq_len(n)) {
    for (r in seq_len(n)) {
      if (c <= r) next

      hits <- logical(nsim)
      for (j in seq_len(nsim)) {
        ok_all_metrics <- TRUE

        for (i in seq_along(metrics)) {
          d <- listc[[j]][[i]][r, c]
          if (!is_extreme(d, is_intersection_metric[i], is_lower)) {
            ok_all_metrics <- FALSE
            break
          }
        }

        hits[j] <- ok_all_metrics
      }

      p_rc <- (1 + sum(hits)) / (nsim + 1)
      P_pair[r, c] <- p_rc
      P_pair[c, r] <- p_rc
    }
  }
  diag(P_pair) <- NA

  ## ---- BH correction (over unique pairwise tests only) ----
  vals <- P_pair[upper_idx]
  keep <- is.finite(vals)

  adj_vals <- vals
  adj_vals[keep] <- p.adjust(vals[keep], method = "BH")

  P_pair_BH <- make_blank()
  P_pair_BH[upper_idx] <- adj_vals
  P_pair_BH[t(upper_idx)] <- adj_vals
  diag(P_pair_BH) <- NA

  ## ---- Global combined p-value (scalar, +1 correction) ----
  # "global hit" if, for each metric, ALL unique pairwise diffs satisfy the appropriate extreme condition
  hits_global <- logical(nsim)

  for (j in seq_len(nsim)) {
    ok_all_metrics <- TRUE

    for (i in seq_along(metrics)) {
      M  <- as.matrix(listc[[j]][[i]])
      vv <- as.numeric(M[upper_idx])  # unique pairs only

      vv <- vv[is.finite(vv)]
      if (length(vv) == 0) { ok_all_metrics <- FALSE; break }

      if (!is_intersection_metric[i]) {
        # distances: sim <= obs  <=>  diff >= 0
        if (!all(vv >= 0)) { ok_all_metrics <- FALSE; break }
      } else {
        if (is_lower) {
          # intersections Lower: sim <= obs  <=>  diff >= 0
          if (!all(vv >= 0)) { ok_all_metrics <- FALSE; break }
        } else {
          # intersections Higher: sim >= obs  <=>  diff <= 0
          if (!all(vv <= 0)) { ok_all_metrics <- FALSE; break }
        }
      }
    }

    hits_global[j] <- ok_all_metrics
  }

  P_global <- (1 + sum(hits_global)) / (nsim + 1)

  ## ---- Return ----
  list(
    P_values        = P_pair,
    P_values_BH     = P_pair_BH,
    P_values_global = P_global
  )
}

