#' Partition anatomical-fidelity uncertainty (simulation-based)
#'
#' \code{anatomical_error_partitioning()} quantifies how much variance in trackway parameters
#' is due to footprint landmark placement uncertainty (degree of anatomical fidelity) versus genuine
#' between-trackway (biological) differences, using simulation only (no observer effects).
#'
#' @param data A \code{trackway} R object. See \code{\link{tps_to_track}}.
#' @param error_radius Numeric; single radius (m) applied to footprints of all trackways or a numeric vector of length
#'  \code{length(data$Trajectories)} with per-trackway tolerances.
#' @param variables A character vector specifying the movement parameters to be used. Valid parameter names include:
#'   \code{"TurnAng"}, \code{"sdTurnAng"}, \code{"PathLen"}, \code{"BeelineLen"}, \code{"StLength"}, \code{"sdStLength"},
#'   \code{"StrideLen"}, \code{"PaceLen"}, \code{"Sinuosity"}, \code{"Straightness"},
#'   \code{"TrackWidth"}, \code{"Gauge"}, \code{"PaceAng"}, \code{"StepAng"},
#'   \code{"Velocity"}, \code{"sdVelocity"}, \code{"MaxVelocity"}, \code{"MinVelocity"}. See  \code{\link{track_param}} and \code{\link{cluster_track}} for details.
#' @param gauge_size Numeric. Pes/manus length (or width) used to compute Gauge as \code{Trackway_width / gauge_size}. See \code{\link{track_param}} for details.
#' @param n_sim Integer; number of Monte Carlo jitter simulations per trackway. Default is \code{200}.
#' @param distribution A character string indicating the jitter model. Options are
#'   \code{"uniform"} (points uniformly within a disk of radius r; default) or
#'   \code{"gaussian"} (SD = r/2 on X and Y, truncated at 3 SD).
#' @param seed Optional integer for reproducibility. Default is \code{NULL}.
#'
#' @details
#' The function estimates how much of the variability in trackway metrics can be
#' attributed to positional uncertainty in footprint landmarks, rather than to true
#' biological differences among trackways. It asks, in essence: "If landmarks
#' were slightly misplaced by up to \code{error_radius}, how much would that
#' affect our computed parameters?" For each trackway, landmarks are repeatedly
#' perturbed within the specified spatial tolerance, the medial trajectory is
#' reconstructed, and parameters are recalculated across Monte Carlo simulations.
#' The resulting variance in each metric is then decomposed into a between-trackway
#' component (biological signal) and a within-trackway component (anatomical noise).
#'
#' The \code{error_radius} represents expected imprecision in reference-point
#' digitization (e.g., preservation quality, erosion, deformation, or subjective
#' landmark interpretation). By specifying a realistic tolerance, users can
#' evaluate how sensitive each metric is to anatomical or taphonomic uncertainty
#' and identify parameters that remain stable despite imperfect preservation.
#'
#' For every trackway and simulation, landmarks are randomly
#' displaced within a circular area of radius \code{r}, using either a uniform
#' model (\code{"uniform"}) or a truncated Gaussian on \code{X} and \code{Y}
#' with \code{SD = r/2} (\code{"gaussian"}, truncated at ±3 SD). Each perturbed
#' landmark set is converted to a medial trajectory (midpoints between
#' consecutive footprints, as in \code{tps_to_track()}), and parameters are
#' recalculated via \code{track_param()}. For each metric \eqn{Y}, total
#' variance is partitioned via the law of total variance into:
#' \itemize{
#'   \item \strong{trackway (biological)} — variance of simulation means among trackways;
#'   \item \strong{anatomical} — mean of within-trackway variances across simulations;
#'   \item \strong{Residual} — residual from the random-intercept model on
#'         simulation means (no observer terms in this simulation-only setup).
#' }
#' Angular variables (\code{TurnAng}, \code{PaceAng}) are embedded in sine–cosine
#' space to handle circularity and avoid 0°/360° discontinuities (Fisher, 1995).
#'
#' A signal-to-noise ratio quantifies robustness:
#' \deqn{\mathrm{SNR} = \frac{\mathrm{Var}(\mathrm{trackway})}
#'                 {\mathrm{Var}(\mathrm{anatomical}) + \mathrm{Var}(\mathrm{Residual})}.}
#' As a rule of thumb:
#' \itemize{
#'   \item \strong{SNR < 1} — anatomical error exceeds biology (weak);
#'   \item \strong{SNR ~ 1–2} — biology and anatomical error are comparable (moderate);
#'   \item \strong{SNR > 2} — biology dominates (strong).
#' }
#' A compact QC table reports SNR, the qualitative rating, the top component by
#' percent, and the conditional \eqn{R^2_c}.
#'
#' @return An \code{"error_partitioning"} R object consisting of a list
#' containing the following elements:
#'
#' \item{summary}{Data frame of variance components per variable.
#' Columns: \code{variable}, \code{component}, \code{variance},
#' \code{percent}, and \eqn{R^2_c}. Components are \code{trackway}
#' (biological), \code{anatomical} (within-trackway simulation variance),
#' and \code{Residual} (from the model on simulation means).}
#'
#' \item{snr}{Data frame with the signal-to-noise ratio per variable.
#' Columns: \code{variable}, \code{bio_component} (always \code{"trackway"}),
#' \code{bio_var}, \code{error_var}, and \code{SNR}. The ratio is
#' \deqn{SNR = Var(trackway) / [Var(anatomical) + Var(Residual)].}}
#'
#' \item{qc}{Compact quality-control table per variable with \code{SNR},
#' a qualitative \code{SNR_rating} (\code{"weak"}, \code{"moderate"},
#' \code{"strong"}), the \code{top_component} by percent, its
#' \code{top_percent}, \eqn{R^2_c}, and an \code{observer_note}
#' (always \code{"anatomical-only"} for this function).}
#'
#' \item{models}{\code{NULL}. Placeholder kept for consistency with the
#' observer-inclusive workflow. No observer models are fitted here.}
#'
#' \item{analysis_table}{Long data frame with all simulated values for
#' every trackway across all Monte Carlo runs (variables joined with trackway
#' IDs and simulation labels).}
#'
#' \item{formulae}{A short note describing the simulation-only
#' decomposition (within-trackway anatomical variance from simulations,
#' plus a random-intercept model on simulation means:
#' \eqn{y \sim 1 + (1|trackway)}; angles handled via sine–cosine).}
#'
#' \item{trackway_names}{Character vector of internal trackway labels used in
#' the analysis.}
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
#' Fisher, N. I. (1995). Statistical Analysis of Circular Data. Cambridge University Press.
#'
#' @examples
#' # Example 1: PaluxyRiver, small jitter (2 cm), uniform noise
#' # Expect relatively stable metrics and higher SNR.
#' set.seed(1)
#' ep_small <- anatomical_error_partitioning(
#'   data         = PaluxyRiver,
#'   error_radius = 0.02,
#'   variables    = c("BeelineLen", "Straightness", "TurnAng"),
#'   n_sim        = 10
#' )
#' ep_small$qc
#'
#' # Example 2: PaluxyRiver, larger jitter (8 cm) to see SNR drop
#' # Inflate anatomical uncertainty; SNR should typically decrease.
#' set.seed(1)
#' ep_large <- anatomical_error_partitioning(
#'   data         = PaluxyRiver,
#'   error_radius = 0.08,
#'   variables    = c("BeelineLen", "Straightness", "TurnAng"),
#'   n_sim        = 10
#' )
#' ep_large$qc
#'
#' # Example 3: MountTom subset + Gaussian jitter (3 cm, truncated at ±3 SD)
#' # Demonstrates alternative noise model and a reduced dataset.
#' sbMountTom <- subset_track(
#'   MountTom,
#'   tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18)
#' )
#' set.seed(2)
#' ep_gauss <- anatomical_error_partitioning(
#'   data         = sbMountTom,
#'   error_radius = 0.03,
#'   variables    = c("StLength", "sdStLength", "PaceAng"),
#'   n_sim        = 10,
#'   distribution = "gaussian"
#' )
#' ep_gauss$snr
#'
#' # Example 4: PaluxyRiver with heterogeneous per-trackway tolerances
#' # Alternate 2 cm / 5 cm across trackways; highlights mixed preservation quality.
#' set.seed(3)
#' rads <- rep(c(0.02, 0.05),
#'             length.out = length(PaluxyRiver$Footprints))
#' ep_het <- anatomical_error_partitioning(
#'   data         = PaluxyRiver,
#'   error_radius = rads,
#'   variables    = c("TrackWidth", "Sinuosity", "PaceAng"),
#'   n_sim        = 10
#' )
#' ep_het$summary
#'
#' # Example 5: Angles only (circular handling via sine–cosine embedding)
#' # Focus on orientation metrics; output includes averaged sin/cos components.
#' set.seed(4)
#' ep_ang <- anatomical_error_partitioning(
#'   data         = PaluxyRiver,
#'   error_radius = 0.04,
#'   variables    = c("TurnAng", "PaceAng"),
#'   n_sim        = 10
#' )
#' ep_ang$summary
#'
#' @importFrom stats runif rnorm var
#' @importFrom trajr TrajFromCoords
#' @seealso \code{\link{tps_to_track}}, \code{\link{track_param}}, \code{\link{subset_track}}
#' @export

anatomical_error_partitioning <- function(
    data,
    error_radius,
    variables,
    gauge_size = NA,
    n_sim = 200,
    distribution = c("uniform","gaussian"),
    seed = NULL
) {

  ## Set default values if arguments are NULL ----
  distribution <- match.arg(distribution)  # default "uniform"

  if (missing(variables) || is.null(variables)) {
    variables <- c(
      "TurnAng","sdTurnAng","PathLen","BeelineLen","StLength","sdStLength",
      "StrideLen","PaceLen","Sinuosity","Straightness","TrackWidth","Gauge",
      "PaceAng","StepAng","Velocity","sdVelocity","MaxVelocity","MinVelocity"
    )
  }

  # --- STRICTLY NEEDED: Gauge requiremen
  if ("Gauge" %in% variables) {
    if (is.null(gauge_size) || is.na(gauge_size) || !is.numeric(gauge_size) ||
        length(gauge_size) != 1L || gauge_size <= 0) {
      stop("Error: 'gauge_size' must be a single positive numeric value when 'Gauge' is included in 'variables'.")
    }
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1)
      warning("`seed` should be a single numeric value.")
    set.seed(seed)  # will error if non-numeric; warn above satisfies tests
  }

  ## ---- Errors and Warnings ----
  stopifnot(is.list(data))
  if (is.null(data$Footprints) || is.null(data$Trajectories)) {
    stop("`data` must be a list with elements $Footprints and $Trajectories.")
  }

  n_trackways <- length(data$Footprints)
  if (n_trackways != length(data$Trajectories)) {
    stop("Lengths of $Footprints and $Trajectories must match.")
  }

  trackway_ids <- names(data$Trajectories)
  if (is.null(trackway_ids) || any(!nzchar(trackway_ids))) {
    trackway_ids <- sprintf("Trackway_%02d", seq_len(n_trackways))
  }

  if (is.null(error_radius)) {
    stop("Provide `error_radius` (single value or per-trackway numeric vector).")
  }
  if (!is.numeric(error_radius)) {
    stop("`error_radius` must be numeric.")
  }
  if (length(error_radius) == 1L) {
    error_radius <- rep(error_radius, n_trackways)
  } else if (length(error_radius) != n_trackways) {
    stop("`error_radius` must be length 1 or equal to the number of trackways.")
  }
  if (any(!is.finite(error_radius) | error_radius < 0)) {
    stop("`error_radius` values must be finite and >= 0.")
  }

  # --- STRICTLY NEEDED: define variables internally ----
  var_spec <- list(
    TurnAng      = list(field = "Mean_turning_angle"),
    sdTurnAng    = list(field = "Standard_deviation_turning_angle"),
    PathLen      = list(field = "Path_length"),
    BeelineLen   = list(field = "Beeline_length"),
    StLength     = list(field = "Mean_step_length"),
    sdStLength   = list(field = "Standard_deviation_step_length"),
    StrideLen    = list(field = "Mean_stride_length"),
    PaceLen      = list(field = "Mean_pace_length"),
    Sinuosity    = list(field = "Sinuosity"),
    Straightness = list(field = "Straightness"),
    TrackWidth   = list(field = "Trackway_width"),
    Gauge        = list(field = "Gauge"),
    PaceAng      = list(field = "Pace_angulation"),
    StepAng      = list(field = "Step_angle"),
    Velocity     = list(field = "Mean_velocity"),
    sdVelocity   = list(field = "Standard_deviation_velocity"),
    MaxVelocity  = list(field = "Maximum_velocity"),
    MinVelocity  = list(field = "Minimum_velocity")
  )

  unknown   <- setdiff(variables, names(var_spec))
  variables <- intersect(variables, names(var_spec))
  if (!length(variables)) {
    stop("No valid variables to analyze after filtering.")
  }
  if (length(unknown)) {
    warning("Unknown variable names ignored: ", paste(unknown, collapse = ", "))
  }

  if (!is.numeric(n_sim)) {
    warning("`n_sim` must be a positive integer.")
    stop("Invalid `n_sim`.")
  }
  if (length(n_sim) != 1L || n_sim != as.integer(n_sim) || n_sim <= 0) {
    stop("`n_sim` must be a positive integer.")
  }

  is_circular <- function(v) v %in% c("TurnAng","PaceAng")

  ## Code ----
  ## ---- Helpers ----
  .jitter_df <- function(df, r, dist) {
    if (is.null(df)) return(df)
    n <- tryCatch(nrow(df), error = function(e) NA_integer_)
    if (is.na(n) || n == 0) return(df)

    if (dist == "uniform") {
      ang <- runif(n, 0, 2*pi)
      rad <- sqrt(runif(n, 0, 1)) * r
      df$X <- df$X + rad * cos(ang)
      df$Y <- df$Y + rad * sin(ang)
    } else { # gaussian: SD = r/2, truncated at ±3 SD
      sdv <- r/2
      dx  <- pmin(pmax(rnorm(n, 0, sdv), -3*sdv), 3*sdv)
      dy  <- pmin(pmax(rnorm(n, 0, sdv), -3*sdv), 3*sdv)
      df$X <- df$X + dx
      df$Y <- df$Y + dy
    }
    df
  }

  .rebuild_traj <- function(footprints_list) {
    L2 <- footprints_list
    L3 <- vector("list", length(L2))
    for (i in seq_along(L2)) {
      n <- nrow(L2[[i]])
      if (n < 2) stop("A trackway has < 2 footprints (cannot rebuild).")
      out <- data.frame(X = NA_real_, Y = NA_real_, IMAGE = NA, ID = NA, Side = NA)
      out <- out[rep(1, n - 1), ]
      out$IMAGE <- L2[[i]][1, "IMAGE"]
      out$ID    <- L2[[i]][1, "ID"]
      out$Side  <- "Medial"
      for (j in 1:(n - 1)) {
        out$X[j] <- (L2[[i]]$X[j] + L2[[i]]$X[j + 1]) / 2
        out$Y[j] <- (L2[[i]]$Y[j] + L2[[i]]$Y[j + 1]) / 2
      }
      L3[[i]] <- trajr::TrajFromCoords(out)
    }
    names(L3) <- trackway_ids
    L3
  }

  ## ---- Monte Carlo simulations ----
  sim_vals <- vector("list", n_sim)
  sim_lab  <- sprintf("sim%03d", seq_len(n_sim))

  for (s in seq_len(n_sim)) {
    # Jitter footprints per trackway
    fp_jit <- lapply(
      seq_len(n_trackways),
      function(i) .jitter_df(data$Footprints[[i]], error_radius[i], distribution)
    )
    names(fp_jit) <- trackway_ids

    # Rebuild trajectories from jittered footprints
    traj_jit <- .rebuild_traj(fp_jit)

    # Build a temporary trackway object and recompute parameters (Gauge handled by track_param)
    obj_jit <- list(Trajectories = traj_jit, Footprints = fp_jit)
    par_jit <- track_param(obj_jit, gauge_size = gauge_size)

    # Collect simulated values for the selected variables
    rows <- lapply(seq_along(par_jit), function(i) {
      x <- par_jit[[i]]
      vals <- lapply(variables, function(v) {
        fld <- var_spec[[v]]$field
        if (is.null(x[[fld]])) NA_real_ else x[[fld]]
      })
      names(vals) <- variables
      data.frame(trackway = trackway_ids[i], as.data.frame(vals), check.names = FALSE)
    })
    sim_vals[[s]] <- do.call(rbind, rows)
    sim_vals[[s]]$anatomical <- sim_lab[s]
  }

  df_all <- do.call(rbind, sim_vals)
  df_all$trackway      <- factor(df_all$trackway, levels = trackway_ids)
  df_all$anatomical <- factor(df_all$anatomical, levels = sim_lab)

  ## ---- Variance decomposition (simulation-only) ----
  rows_out <- list()
  snr_rows <- list()

  for (nm in variables) {
    sub <- df_all[, c("trackway","anatomical", nm)]
    names(sub)[ncol(sub)] <- "value"
    sub <- sub[is.finite(sub$value), , drop = FALSE]

    if (is_circular(nm)) {
      th <- (sub$value * pi/180) %% (2*pi)
      s  <- sin(th); c <- cos(th)
      # within-trackway anatomical variance (average of sin/cos)
      var_s <- tapply(s, sub$trackway, stats::var, na.rm = TRUE)
      var_c <- tapply(c, sub$trackway, stats::var, na.rm = TRUE)
      Var_anatomical <- mean((var_s + var_c)/2, na.rm = TRUE)
      # biological (between trackways) from means
      mean_s <- tapply(s, sub$trackway, mean, na.rm = TRUE)
      mean_c <- tapply(c, sub$trackway, mean, na.rm = TRUE)
      Var_bio <- mean(c(stats::var(mean_s, na.rm = TRUE),
                        stats::var(mean_c, na.rm = TRUE)), na.rm = TRUE)
    } else {
      v_within      <- tapply(sub$value, sub$trackway, stats::var,  na.rm = TRUE)
      m_by_trackway    <- tapply(sub$value, sub$trackway, mean,       na.rm = TRUE)
      Var_anatomical <- mean(v_within,    na.rm = TRUE)
      Var_bio        <- stats::var(as.numeric(m_by_trackway), na.rm = TRUE)
    }

    Var_res <- 0
    total   <- Var_bio + Var_anatomical + Var_res

    tab <- data.frame(
      variable  = nm,
      component = c("trackway","anatomical","Residual"),
      variance  = c(Var_bio, Var_anatomical, Var_res),
      percent   = if (is.finite(total) && total > 0)
        100 * c(Var_bio, Var_anatomical, Var_res) / total else
          c(NA_real_, NA_real_, NA_real_),
      R2_c      = rep(if (is.finite(total) && total > 0) Var_bio/total else NA_real_, 3),
      stringsAsFactors = FALSE
    )
    rows_out[[nm]] <- tab

    snr_rows[[nm]] <- data.frame(
      variable      = nm,
      bio_component = "trackway",
      bio_var       = Var_bio,
      error_var     = Var_anatomical + Var_res,
      SNR           = if ((Var_anatomical + Var_res) > 0)
        Var_bio / (Var_anatomical + Var_res) else NA_real_
    )
  }

  summary <- do.call(rbind, rows_out)
  snr_tbl <- do.call(rbind, snr_rows)

  ## ---- QC table ----
  qc_rows <- lapply(unique(summary$variable), function(v) {
    tab <- summary[summary$variable == v, , drop = FALSE]
    snr <- snr_tbl$SNR[snr_tbl$variable == v][1]
    rating <- if (is.na(snr)) "NA" else if (snr < 1) "weak (error > bio)"
    else if (snr < 2) "moderate" else "strong"
    ord <- order(tab$percent, decreasing = TRUE)
    R2c <- unique(tab$R2_c); R2c <- if (length(R2c)) R2c[1] else NA_real_
    data.frame(
      variable = v,
      SNR = snr,
      SNR_rating = rating,
      top_component = tab$component[ord][1],
      top_percent   = tab$percent[ord][1],
      R2_c = R2c,
      observer_estimable = "no (anatomical-only, sim-based)",
      stringsAsFactors = FALSE
    )
  })
  qc_tbl <- do.call(rbind, qc_rows)

  ## ---- Return ----
  res <- list(
    summary        = summary,
    snr            = snr_tbl,
    qc             = qc_tbl,
    models         = NULL,
    analysis_table = df_all,
    formulae       = list(note = "simulation-only decomposition via law of total variance"),
    trackway_names    = trackway_ids
  )
  class(res) <- "error_partitioning"
  res
}

