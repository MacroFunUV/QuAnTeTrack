#' Test for differences in direction means with pairwise comparisons
#'
#' \code{test_direction()} assesses differences in mean direction across tracks using a selected circular statistical test. It provides two options: the parametric Watson-Williams test (default), which assumes that directions follow a von Mises distribution and that dispersion is homogeneous across tracks—these assumptions are automatically checked before testing—and the non-parametric Watson-Wheeler test, which does not require these assumptions. For datasets with more than two tracks, the function also performs pairwise comparisons to identify which tracks differ significantly in mean direction.
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param analysis A character string specifying the type of circular analysis: \code{"Watson-Williams"} (default; parametric, assumes von Mises and similar concentration) or \code{"Watson-Wheeler"} (non-parametric).
#' @param permutation Logical or \code{NULL}. For \code{"Watson-Wheeler"} only, whether to compute a permutation (Monte-Carlo) *p*-value. If \code{NULL} (default), the function automatically switches to permutation when some groups have \eqn{n < 10} or ties are detected; otherwise honors the user choice.
#' @param B Integer. Number of permutations for the Watson-Wheeler permutation *p*-value (default \code{1000}).
#' @param seed Optional integer seed for reproducibility when using permutation (default \code{NULL}).
#'
#' @details
#' The \code{test_direction()} function performs circular data analyses using functions primarily from the \pkg{circular} package. Estimation of the concentration parameter (\eqn{\kappa}) relies on the \code{est.kappa()} function of the \pkg{CircStats} package. It includes:
#'
#' - **Condition Testing:**
#'   - **Non-uniformity within tracks:** Rayleigh test on step directions within each track.
#'   - **Similarity of concentrations:** Per-track concentration \eqn{(\kappa)} is estimated via \code{est.kappa}. Large heterogeneity (e.g., \eqn{\max(\kappa)/\min(\kappa) > 2}) or unstable estimates will trigger a warning, since the Watson-Williams test assumes approximately similar concentration across tracks.
#'
#' - **Statistical Analysis:**
#'   - **Watson-Williams (default):** Parametric comparison of mean directions across tracks (assumes von Mises and similar concentration).
#'   - **Watson-Wheeler:** Non-parametric k-sample test for differences in direction. When some groups have \eqn{n < 10} or when there are ties (identical angles), the asymptotic chi-squared *p*-value can be unreliable. To address this, the function can compute a permutation (Monte-Carlo) *p*-value for the same Watson-Wheeler statistic (W). If \code{permutation = NULL} (default), this permutation *p*-value is used automatically under small samples or ties; otherwise, the user's choice is respected. A tiny internal jitter (in radians) is applied to angles only when ties are detected to stabilize the test; this is not exposed as an argument.
#' - **Permutation (Monte-Carlo) p-values:** For the Watson-Wheeler test, when group sizes are small (\eqn{n < 10}) or when ties are present, the default chi-squared *p*-value can be unreliable. In these cases, the function can compute a permutation (Monte-Carlo) *p*-value. This is done by repeatedly shuffling group labels, recalculating the Watson-Wheeler statistic each time, and comparing the observed value to the resulting null distribution. The *p*-value is then the proportion of permuted statistics at least as extreme as the observed one. This approach provides a more accurate significance level for small or tied samples.
#'   - **Pairwise comparisons:** For datasets with more than two tracks, pairwise tests are performed using the same family as the selected global test. *P*-values are adjusted using Holm's method. For Watson-Wheeler, permutation *p*-values are used pairwise whenever applicable by the same rule (small n or ties, or when \code{permutation = TRUE}).
#'
#' - **Direction Measurement:**
#'   - The direction is measured in degrees. The reference direction, or 0 degrees, is considered to be along the positive x-axis. Angles are measured counterclockwise from the positive x-axis, with 0 degrees pointing directly along this axis.
#'
#' @return A list with the results of the statistical analysis and diagnostic tests:
#'   - \code{assumption_results}: A list with per-track Rayleigh test results (statistic and *p*-value), estimated \eqn{\kappa} by track, and summary diagnostics (\code{kappa_range}, \code{kappa_ratio}).
#'   - \code{global_test}: Result of the selected k-sample circular test (Watson-Williams or Watson-Wheeler). If permutation was used for Watson-Wheeler, the object is an \code{"htest"} with the permutation *p*-value and the method labeled as \code{"Watson-Wheeler (permutation, B=...)"}.
#'   - \code{pairwise}: Data frame of pairwise comparisons (test statistic, raw *p*-value, and Holm-adjusted *p*-value) when more than two tracks are present.
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
#' # Example 1: Parametric Circular Test (Watson-Williams) for Differences in Direction Means with Pairwise Comparisons in MountTom dataset
#' test_direction(MountTom, analysis = "Watson-Williams")
#'
#' # Example 2: Non-Parametric Circular Test (Watson-Wheeler) with automatic permutation under small n/ties
#' test_direction(MountTom, analysis = "Watson-Wheeler",
#'                permutation = TRUE, B = 10, seed = 42)
#'
#' # Example 3: Force permutation with more replicates and a seed
#' test_direction(PaluxyRiver, analysis = "Watson-Wheeler",
#'                permutation = TRUE, B = 100, seed = 42)
#'
#' # Example 4: Keep asymptotic p-values (even if small n or ties), but warn
#' test_direction(MountTom, analysis = "Watson-Wheeler", permutation = FALSE)
#'
#' @importFrom circular circular watson.williams.test watson.wheeler.test rayleigh.test rho.circular
#' @importFrom CircStats est.kappa
#' @importFrom stats p.adjust runif
#'
#' @seealso \code{\link{tps_to_track}}, \code{\link{plot_direction}}
#'
#' @export
test_direction <- function(data, analysis = NULL,
                           permutation = NULL, B = 1000, seed = NULL) {

  ## --- Internal helpers (ties + permutation) ----
  .ww_stat <- function(circ_list) {
    suppressWarnings(circular::watson.wheeler.test(circ_list)$statistic)
  }
  .has_ties <- function(circ_list) {
    any(vapply(circ_list, function(z) {
      zz <- as.numeric(z)  # radians
      any(duplicated(signif(zz, 12)))
    }, logical(1)))
  }
  .apply_jitter_internal <- function(circ_list) {
    # tiny jitter (radians) used ONLY if ties are present
    eps <- 1e-6
    lapply(circ_list, function(z) {
      zr <- as.numeric(z)
      zr <- zr + stats::runif(length(zr), min = -eps, max = eps)
      circular::circular(zr, units = "radians", modulo = "2pi", zero = 0, rotation = "counter")
    })
  }
  .ww_permutation <- function(circ_list, B = 1000, seed = NULL) {
    angles <- unlist(lapply(circ_list, as.numeric))          # radians
    sizes  <- vapply(circ_list, length, integer(1))
    groups <- names(circ_list)
    W_obs  <- .ww_stat(circ_list)
    if (!is.null(seed)) set.seed(seed)
    W_perm <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(length(angles), replace = FALSE)
      split_idx <- split(idx, rep.int(seq_along(sizes), sizes))
      perm_list <- lapply(split_idx, function(i) circular::circular(angles[i],
                                                                    units = "radians",
                                                                    modulo = "2pi",
                                                                    zero = 0, rotation = "counter"))
      names(perm_list) <- groups
      W_perm[b] <- .ww_stat(perm_list)
    }
    p_emp <- (sum(W_perm >= W_obs) + 1) / (B + 1)
    list(statistic = W_obs, p_value = p_emp, method = sprintf("Watson-Wheeler (permutation, B=%d)", B))
  }

  ## Set default values if arguments are NULL----
  if (is.null(analysis)) analysis <- "Watson-Williams" # Default to parametric circular test

  ## Errors and Warnings----

  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'track' R object, which is a list consisting of two elements: 'Trajectories' and 'Footprints'.")
  }
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("Both elements of 'data' must be lists. Ensure that 'Trajectories' and 'Footprints' are provided.")
  }
  if (any(sapply(data, function(x) length(x) == 0))) {
    stop("The elements within 'data' must not be empty. Ensure that both 'Trajectories' and 'Footprints' contain data.")
  }
  if (!analysis %in% c("Watson-Williams", "Watson-Wheeler")) {
    stop("Invalid 'analysis' argument. Choose 'Watson-Williams' (default) or 'Watson-Wheeler'.")
  }

  ## Code----

  track_param <- track_param(data)
  data <- data[[1]]

  n <- c(track_param[[1]][[1]], track_param[[1]][[1]][[length(track_param[[1]][[1]])]])
  if (length(data) > 1) {
    for (i in 2:length(data)) {
      n <- c(n, c(track_param[[i]][[1]], track_param[[i]][[1]][[length(track_param[[i]][[1]])]]))
    }
  }

  M <- data.frame(dir = n, track = rep(names(data), sapply(data, nrow)))
  M$track <- gsub("_", " ", M$track)
  M_original <- M

  track_counts <- table(M$track)
  valid_tracks <- names(track_counts[track_counts > 3])
  invalid_tracks <- names(track_counts[track_counts <= 3])
  if (length(invalid_tracks) > 0) {
    warning("The following tracks were removed from the analysis due to having 3 or fewer footprints: ", paste(invalid_tracks, collapse = ", "), ".")
  }
  if (length(valid_tracks) < 2) {
    stop("Not enough tracks with more than 3 footprints for meaningful analysis. Ensure at least two tracks have more than 3 footprints.")
  }
  M_analysis <- subset(M, track %in% valid_tracks)

  ## Circular assumption checks & data prep ----
  split_dirs <- split(M_analysis$dir, M_analysis$track)
  circ_list <- lapply(
    split_dirs,
    function(x) circular::circular(as.numeric(x),
                                   units = "degrees",
                                   modulo = "2pi",
                                   zero = 0,
                                   rotation = "counter")
  )

  # Rayleigh per track
  rayleigh_list <- lapply(circ_list, function(z) suppressWarnings(circular::rayleigh.test(z)))
  rayleigh_results <- t(vapply(rayleigh_list, function(obj) {
    c(statistic = unname(obj$statistic), p_value = obj$p.value)
  }, numeric(2)))
  rownames(rayleigh_results) <- names(circ_list)

  # Kappa per track
  kappas <- vapply(circ_list, function(z) {
    kz <- suppressWarnings(try(CircStats::est.kappa(as.numeric(z)), silent = TRUE))
    if (inherits(kz, "try-error") || !is.finite(kz)) NA_real_ else as.numeric(kz)
  }, numeric(1))

  kappa_range <- if (all(is.na(kappas))) NA_real_ else (max(kappas, na.rm = TRUE) - min(kappas, na.rm = TRUE))
  kappa_ratio <- if (all(is.na(kappas)) || min(kappas, na.rm = TRUE) <= 0) NA_real_ else max(kappas, na.rm = TRUE) / min(kappas[kappas > 0], na.rm = TRUE)

  assumption_results <- list(
    rayleigh = rayleigh_results,
    kappa = kappas,
    kappa_range = kappa_range,
    kappa_ratio = kappa_ratio
  )

  if (any(!is.na(assumption_results$rayleigh[, "p_value"]) & assumption_results$rayleigh[, "p_value"] > 0.05)) {
    warning("One or more tracks show near-uniform directions (Rayleigh p > 0.05). Parametric assumptions (von Mises) may be questionable; consider 'Watson-Wheeler'.")
  }
  if (any(is.na(kappas))) {
    warning("Some tracks have unstable concentration (kappa) estimates; Watson-Williams assumptions may not hold. Consider using 'Watson-Wheeler'.")
  } else if (is.finite(kappa_ratio) && kappa_ratio > 2) {
    warning("Estimated concentration (kappa) appears heterogeneous across tracks (ratio > 2). Watson-Williams assumes similar concentration; consider 'Watson-Wheeler'.")
  }

  results <- list(assumption_results = assumption_results)

  ## Decide about small n and ties (for WW decisions) ----
  n_per_group <- vapply(circ_list, length, integer(1))
  small_group <- any(n_per_group < 10)
  has_ties    <- .has_ties(circ_list)

  ## Statistical analysis ----
  if (analysis == "Watson-Williams") {
    all_rad <- unlist(lapply(circ_list, as.numeric))
    all_circ <- circular::circular(all_rad, units = "radians",
                                   modulo = "2pi", zero = 0, rotation = "counter")
    groups <- factor(rep(names(circ_list), times = vapply(circ_list, length, integer(1))))
    results$global_test <- suppressWarnings(circular::watson.williams.test(all_circ, groups))

    pairs <- utils::combn(names(circ_list), 2, simplify = FALSE)
    pw <- lapply(pairs, function(pp) {
      a <- circ_list[[pp[1]]]; b <- circ_list[[pp[2]]]
      xx <- c(as.numeric(a), as.numeric(b))
      grp <- factor(c(rep(pp[1], length(a)), rep(pp[2], length(b))))
      ww <- suppressWarnings(circular::watson.williams.test(
        circular::circular(xx, units = "radians", modulo = "2pi", zero = 0, rotation = "counter"),
        grp
      ))
      data.frame(track1 = pp[1], track2 = pp[2],
                 statistic = unname(ww$statistic),
                 p_value = ww$p.value,
                 method = "Watson-Williams",
                 stringsAsFactors = FALSE)
    })
    pairwise_df <- do.call(rbind, pw)
    if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
      pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "holm")
    }
    results$pairwise <- pairwise_df

  } else if (analysis == "Watson-Wheeler") {

    if (!is.null(seed)) set.seed(seed)

    # apply internal jitter ONLY if ties present
    circ_list_j <- if (has_ties) .apply_jitter_internal(circ_list) else circ_list

    # decide permutation use:
    use_perm <- if (is.null(permutation)) (small_group || has_ties) else isTRUE(permutation)

    if (use_perm) {
      if (small_group) warning("Using permutation p-value for Watson-Wheeler because some groups have n < 10.")
      if (has_ties && is.null(permutation))
        warning("Ties detected. Using permutation p-value for Watson-Wheeler.")
      ww_perm <- .ww_permutation(circ_list_j, B = B, seed = seed)
      results$global_test <- structure(list(statistic = ww_perm$statistic,
                                            p.value   = ww_perm$p_value,
                                            method    = ww_perm$method),
                                       class = "htest")
    } else {
      results$global_test <- suppressWarnings(circular::watson.wheeler.test(circ_list_j))
    }

    pairs <- utils::combn(names(circ_list), 2, simplify = FALSE)
    pw <- lapply(pairs, function(pp) {
      a <- circ_list[[pp[1]]]; b <- circ_list[[pp[2]]]
      ab_list <- list(a, b)
      ties_pair <- .has_ties(ab_list)
      ab_list_j <- if (ties_pair) .apply_jitter_internal(ab_list) else ab_list
      small_pair <- any(vapply(ab_list, length, integer(1)) < 10)

      use_perm_pair <- if (is.null(permutation)) (small_pair || ties_pair) else isTRUE(permutation)

      if (use_perm_pair) {
        ww <- .ww_permutation(ab_list_j, B = B, seed = seed)
        data.frame(track1 = pp[1], track2 = pp[2],
                   statistic = unname(ww$statistic),
                   p_value = ww$p_value,
                   method = sprintf("Watson-Wheeler (permutation, B=%d)", B),
                   stringsAsFactors = FALSE)
      } else {
        tt <- suppressWarnings(circular::watson.wheeler.test(ab_list_j))
        data.frame(track1 = pp[1], track2 = pp[2],
                   statistic = unname(tt$statistic),
                   p_value = tt$p.value,
                   method = "Watson-Wheeler",
                   stringsAsFactors = FALSE)
      }
    })
    pairwise_df <- do.call(rbind, pw)
    if (!is.null(pairwise_df) && nrow(pairwise_df) > 0) {
      pairwise_df$p_adj <- p.adjust(pairwise_df$p_value, method = "holm")
    }
    results$pairwise <- pairwise_df
  }

  return(results)
}
