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
#'
#' @details
#' The \code{combined_prob()} function combines *p*-values derived from multiple similarity metric tests and intersection tests.
#' It calculates the combined *p*-values by assessing the probability of observing the combined metrics across simulated datasets.
#' Pairwise *p*-values use a Monte Carlo tail test with the (+1) correction and are FDR-adjusted with Benjamini–Hochberg (BH)
#' across unique trajectory pairs.
#'
#' @return A list containing:
#' \item{P_values}{A matrix of raw *p*-values for the combined metrics across all trajectories (pairwise).}
#' \item{P_values_BH}{A matrix of BH-adjusted *p*-values for the combined metrics (pairwise).}
#' \item{P_values_global}{A single numeric value: the overall probability of observing the combined metrics across all pairs of trajectories.}
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
#' @examples
#' # Example 1: "Directed" model and similarity metrics.
#' s1 <- simulate_track(PaluxyRiver, nsim = 3, model = "Directed")
#' DTW1 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#' Frechet1 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s1, superposition = "None")
#' int1 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1,
#'   origin.permutation = "None")
#' combined_prob(PaluxyRiver, metrics = list(DTW1, Frechet1, int1))
#'
#' # Example 2: "Constrained" model and similarity metrics.
#' s2 <- simulate_track(PaluxyRiver, nsim = 3, model = "Constrained")
#' DTW2 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s2,
#'   superposition = "None")
#' Frechet2 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s2,
#'   superposition = "None")
#' int2 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2,
#'   origin.permutation = "Min.Box")
#' combined_prob(PaluxyRiver, metrics = list(DTW2, Frechet2, int2))
#'
#' # Example 3: "Unconstrained" model and similarity metrics.
#' s3 <- simulate_track(PaluxyRiver, nsim = 3, model = "Unconstrained")
#' DTW3 <- simil_DTW_metric(PaluxyRiver, test = TRUE, sim = s3,
#'   superposition = "None")
#' Frechet3 <- simil_Frechet_metric(PaluxyRiver, test = TRUE, sim = s3,
#'   superposition = "None")
#' int3 <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3,
#'   origin.permutation = "Conv.Hull")
#' combined_prob(PaluxyRiver, metrics = list(DTW3, Frechet3, int3))
#'
#' @seealso \code{\link{tps_to_track}}, \code{\link{simulate_track}}, \code{\link{track_intersection}}, \code{\link{simil_DTW_metric}}, \code{\link{simil_Frechet_metric}}
#'
#' @export


combined_prob <- function(data, metrics = NULL) {
  ## ---- Checks ----
  if (!is.list(data) || length(data) < 2)
    stop("The 'data' argument must be a 'track' R object (list of two elements).")
  if (!is.list(data[[1]]) || !is.list(data[[2]]))
    stop("The two elements of 'data' must be lists.")
  if (is.null(metrics) || !is.list(metrics))
    stop("'metrics' must be a list of track similarity or intersection metric objects.")

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

  # listc[[j]][[i]] will store (Observed - Sim_j) for metric i, simulation j
  listc <- vector("list", nsim)
  for (j in seq_len(nsim)) {
    per_metric <- vector("list", length(metrics))
    for (i in seq_along(metrics)) {
      obs  <- as.matrix(metrics[[i]][[1]])       # observed matrix
      simj <- as.matrix(metrics[[i]][[5]][[j]])  # j-th simulated matrix
      diffs <- make_blank()
      diffs[upper_idx] <- obs[upper_idx] - simj[upper_idx]
      per_metric[[i]] <- diffs
    }
    listc[[j]] <- per_metric
  }

  ## ---- Pairwise combined p-values (raw, +1 correction) ----
  P_pair <- make_blank()
  for (c in seq_len(n)) for (r in seq_len(n)) {
    if (c <= r) next
    hits <- logical(nsim)
    for (j in seq_len(nsim)) {
      v <- vapply(seq_along(metrics), function(i) listc[[j]][[i]][r, c], numeric(1))
      hits[j] <- all(is.finite(v) & v > 0)
    }
    P_pair[r, c] <- (1 + sum(hits)) / (nsim + 1)
  }

  ## ---- BH correction ----
  vals      <- P_pair[upper_idx]; keep <- is.finite(vals)
  adj_vals  <- vals; adj_vals[keep] <- p.adjust(vals[keep], method = "BH")
  P_pair_BH <- make_blank(); P_pair_BH[upper_idx] <- adj_vals

  ## ---- Global combined p-value (scalar, +1 correction) ----
  hits_global <- logical(nsim)
  for (j in seq_len(nsim)) {
    mats <- do.call(rbind, listc[[j]])
    vv   <- as.numeric(mats); vv <- vv[is.finite(vv)]
    hits_global[j] <- length(vv) > 0 && all(vv > 0)
  }
  P_global <- (1 + sum(hits_global)) / (nsim + 1)

  ## ---- Return ----
  list(
    P_values        = P_pair,
    P_values_BH     = P_pair_BH,
    P_values_global = P_global
  )
}
