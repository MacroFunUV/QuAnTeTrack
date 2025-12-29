#' Partition intra- and inter-observer uncertainty (replicate-based)
#'
#' \code{observer_error_partitioning()} quantifies how much variance in trackway parameters
#' is due to observer and replicate effects (inter-/intra-observer) versus genuine
#' between-track (biological) differences, using replicated digitizations and mixed-effects models
#' (no anatomical simulation).
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param metadata A data frame with one row per track (same order as
#'   \code{track_param(data)}). Required columns:
#'   \itemize{
#'     \item \code{replica} — replicate index within (\code{track}, \code{observer}).
#'     \item \code{track} — biological track ID.
#'     \item \code{observer} — observer ID.
#'   }
#' @param variables A character vector of parameters to analyze. Default is
#'   \code{c("TurnAng","sdTurnAng","Distance","Length","StLength","sdStLength",
#'           "Sinuosity","Straightness","TrackWidth","PaceAng")}.
#'
#' @details
#' This function partitions the variability of each selected metric into
#' components attributable to \strong{biology} (consistent differences among
#' tracks) and to the \strong{human sampling process}: differences between
#' observers (\strong{inter-observer sampling error}) and differences between
#' repeated digitizations by the same observer (\strong{intra-observer sampling error}).
#' Use it to audit reproducibility across people and sessions, identify metrics
#' that remain stable despite who digitizes them, and prioritize where protocol
#' changes (clearer landmark definitions, training, calibration) will most reduce
#' measurement noise. Metrics whose variability is dominated by biological signal
#' are more interpretable downstream (testing, clustering, classification),
#' whereas those dominated by observer/replica effects warrant caution or
#' methodological refinement.
#'
#' For each metric the model fits random intercepts for
#' \code{track}, \code{observer}, and their \code{observer:track} interaction.
#' Replicate-to-replicate scatter within the same (\code{observer}, \code{track})
#' cell contributes to the residual (this is the intra-observer sampling error).
#' Random terms that appear with a single level are dropped automatically.
#'
#' \deqn{y \sim 1 + (1|track) + (1|observer) + (1|observer:track).}
#'
#' Angular variables (\code{TurnAng}, \code{PaceAng}) are handled in sine–cosine
#' space to respect circularity and avoid 0°/360° discontinuities: the same
#' random-effects structure is fit separately to \eqn{\sin(\theta)} and \eqn{\cos(\theta)},
#' and the two variance decompositions are averaged (Fisher, 1995).
#'
#' Robustness is summarized by the signal-to-noise ratio
#' \deqn{SNR = \frac{Var(track)}{Var(observer) + Var(observer:track) + Var(Residual)}.}
#' As a rule of thumb:
#' \itemize{
#'   \item \strong{SNR < 1} — observer/replicate error exceeds biology; the metric is \emph{weak};
#'   \item \strong{SNR ~ 1–2} — biology and error are comparable; \emph{moderate} robustness;
#'   \item \strong{SNR > 2} — biology dominates; \emph{strongly robust}.
#' }
#' A compact QC table reports SNR, the qualitative rating, the top variance
#' component by percent, and the conditional \eqn{R^2_c}; it also flags whether
#' the inter-observer term was estimable (>1 level).
#'
#' Practical designs:
#' \itemize{
#'   \item \emph{Full inter + intra}: provide multiple observers \emph{and} multiple
#'         replicates per observer. Both \code{observer} and \code{observer:track} are estimable; replicate scatter is residual.
#'   \item \emph{Intra-only}: set \code{metadata$observer} to the same value for all rows
#'         and provide \code{replica} with >1 levels. Only \code{track} remains; replicate scatter is residual.
#'   \item \emph{Inter-only}: set \code{metadata$replica} to a single level and supply multiple observers.
#'         \code{observer} is estimable; the residual still represents replicate-level noise (if present).
#' }
#'
#' Note that estimates can be sensitive when designs are highly
#' unbalanced or sample sizes are small; singular fits (random variances at/near
#' zero) may occur (REML) and should be interpreted cautiously.
#' This function does not simulate anatomical landmark jitter; to quantify
#' sensitivity to landmark placement itself (e.g., preservation, landmark
#' ambiguity), use \code{anatomical_error_partitioning()}.
#'
#' All \strong{replicated digitizations} (across observers and/or replicate
#' attempts) must be bundled together in this single \code{data}
#' object (i.e., in the same file), with each replicate represented
#' as its own track entry. Do not split replicates across multiple
#' files/objects, as \code{track_param(data)} and \code{metadata}
#' must align one row per replicate in the same order.
#'
#' @return An \code{"error_partitioning"} R object consisting of a list
#' containing the following elements:
#'
#' \item{summary}{Data frame of variance components per variable.
#' Columns: \code{variable}, \code{component}, \code{variance},
#' \code{percent}, and \eqn{R^2_c}. Components reflect what was
#' estimable (e.g., \code{track}, \code{observer},
#' \code{observer:track}, \code{Residual}).}
#'
#' \item{snr}{Data frame with the signal-to-noise ratio per variable.
#' Columns: \code{variable}, \code{bio_component} (always \code{"track"}),
#' \code{bio_var}, \code{error_var}, and \code{SNR}.
#' The SNR is
#' \deqn{SNR = Var(track) / [Var(observer) + Var(observer:track) + Var(Residual)].}}
#'
#' \item{qc}{Compact quality-control table per variable with \code{SNR},
#' a qualitative \code{SNR_rating} (\code{"weak"}, \code{"moderate"},
#' \code{"strong"}), the \code{top_component} by percent, its
#' \code{top_percent}, \eqn{R^2_c}, and whether the observer term was
#' estimable (\code{observer_estimable}).}
#'
#' \item{analysis_table}{Data frame used in the analysis: selected
#' variables joined with \code{metadata} (columns typically include
#' \code{replica}, \code{track}, and \code{observer}).}
#'
#' \item{models}{Fitted mixed-effects models used to extract variance
#' components. For linear variables, each entry stores the \code{lmer}
#' fit and its variance components. For angular variables, each entry is
#' a list with separate \code{sin} and \code{cos} fits.}
#'
#' \item{formulae}{Model formula(e) per variable. For angular variables,
#' a list containing the \code{sin} and \code{cos} formulas.}
#'
#' \item{track_names}{Character vector of internal track labels used in
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
#' Fisher, N. I. (1995). Statistical Analysis of Circular Data. Cambridge University Press
#'
#' @examples
#' # Example 1: Full partition (inter + intra via residual)
#' # Model: (1|track) + (1|observer) + (1|observer:track).
#' # Var(track)=biology; Var(observer)=inter-observer;
#' # Var(observer:track)=observer×track; Residual=intra-observer.
#' tps <- system.file("extdata", "PaluxyRiverObsRep.tps",
#'                    package = "QuAnTeTrack")
#' RL <- rep(c("R","L"), length.out = 34)
#' trks <- tps_to_track(file = tps, scale = 0.004341493,
#'                      R.L.side = RL, missing = FALSE, NAs = NULL)
#'
#' n  <- length(trks$Trajectories)
#' md <- data.frame(
#'   track    = rep(c("T01","T02"), length.out = n),
#'   observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
#'   stringsAsFactors = FALSE
#' )
#' md$replica <- ave(seq_len(n),
#'                   interaction(md$observer, md$track, drop = TRUE),
#'                   FUN = seq_along)
#' md <- md[, c("replica","track","observer")]
#'
#' vars <- c("Distance","Straightness","TurnAng","PaceAng")
#' res_full <- observer_error_partitioning(data = trks,
#'                                         metadata = md,
#'                                         variables = vars)
#' res_full$qc
#'
#' # Example 2: Intra-only, single observer; residual = intra-observer
#' md_intra <- md
#' md_intra$observer <- "obs1"
#' res_intra <- observer_error_partitioning(data = trks,
#'                                          metadata = md_intra,
#'                                          variables = vars)
#' res_intra$qc
#'
#' # Example 3: Inter-only, single replicate; estimates inter-observer
#' md_inter <- md
#' md_inter$replica <- 1L
#' res_inter <- observer_error_partitioning(data = trks,
#'                                          metadata = md_inter,
#'                                          variables = vars)
#' res_inter$qc
#'
#' # Example 4: Circular metrics only (angles via sine–cosine embedding)
#' res_ang <- observer_error_partitioning(data = trks,
#'                                        metadata = md,
#'                                        variables = c("TurnAng","PaceAng"))
#' res_ang$summary
#'
#' @importFrom lme4 lmer VarCorr lmerControl
#' @importFrom MuMIn r.squaredGLMM
#' @importFrom stats as.formula sigma
#' @seealso \code{\link{tps_to_track}}, \code{\link{track_param}}
#' @export
observer_error_partitioning <- function(data,
                                        metadata,
                                        variables = c("TurnAng","sdTurnAng","Distance","Length",
                                                      "StLength","sdStLength","Sinuosity","Straightness",
                                                      "TrackWidth","PaceAng")) {

  ## Set default values if arguments are NULL ----
  # (no optional defaults beyond 'variables')

  ## Errors and Warnings----
  # compute parameter list and basic checks
  param_obj   <- track_param(data)
  track_names <- names(param_obj)

  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  req <- c("replica","track","observer")
  if (!all(req %in% names(metadata)))
    stop("`metadata` must include columns: 'replica', 'track', 'observer'.")
  if (nrow(metadata) != length(param_obj))
    stop("`metadata` must have one row per track entry, in the same order as track_param(data).")

  # map short -> long names
  name_map <- c(
    TurnAng      = "Mean_turning_angle",
    sdTurnAng    = "Standard_deviation_turning_angle",
    Distance     = "Distance",
    Length       = "Length",
    StLength     = "Mean_step_length",
    sdStLength   = "Standard_deviation_step_length",
    Sinuosity    = "Sinuosity",
    Straightness = "Straightness",
    TrackWidth   = "Trackway_width",
    PaceAng      = "Pace_angulation"
  )
  unknown   <- setdiff(variables, names(name_map))
  variables <- intersect(variables, names(name_map))

  # All unknown -> error (NO warning)
  if (!length(variables))
    stop("No valid variables to analyze after filtering.")

  # Mixed known/unknown -> warn, then proceed
  if (length(unknown))
    warning("Unknown variable names ignored: ", paste(unknown, collapse = ", "))
  is_circular_short <- function(v) v %in% c("TurnAng","PaceAng")

  ## Code----
  # build analysis table
  rows_list <- lapply(seq_along(param_obj), function(i) {
    x <- param_obj[[i]]
    vals <- lapply(variables, function(short) {
      nm <- name_map[[short]]
      if (is.null(x[[nm]])) NA_real_ else x[[nm]]
    })
    names(vals) <- variables
    data.frame(.row = i, as.data.frame(vals), check.names = FALSE)
  })
  vals_df <- do.call(rbind, rows_list)
  df <- cbind(vals_df, metadata)

  # cast factors
  df$track    <- factor(df$track)
  df$observer <- factor(df$observer)
  df$replica  <- factor(df$replica)

  # ---- RANDOM-EFFECTS SELECTION (Option A, robust) ---------------------------
  n_track <- nlevels(df$track)
  n_obs   <- nlevels(df$observer)

  # interaction groups and replicates per cell
  int_grp <- interaction(df$observer, df$track, drop = TRUE)
  rep_per_cell <- as.integer(table(int_grp))

  include_track <- n_track > 1
  include_obs   <- n_obs   > 1
  # include observer:track only if BOTH observer and track have >1 levels
  # AND at least one cell has >=2 replicates; otherwise it is redundant or not estimable
  include_int   <- (include_obs && include_track &&
                      nlevels(int_grp) > 1 &&
                      any(rep_per_cell >= 2))

  rand_terms <- character(0)
  if (include_track) rand_terms <- c(rand_terms, "(1|track)")
  if (include_obs)   rand_terms <- c(rand_terms, "(1|observer)")
  if (include_int)   rand_terms <- c(rand_terms, "(1|observer:track)")

  if (!length(rand_terms))
    stop("No random effect has > 1 level; cannot partition variance.")
  random_part   <- paste(rand_terms, collapse = " + ")
  full_formula  <- function(resp) as.formula(sprintf("%s ~ 1 + %s", resp, random_part))

  # helpers
  .fit_lmm <- function(fm, dat) {
    lme4::lmer(
      fm, data = dat, REML = TRUE,
      control = lme4::lmerControl(
        check.conv.singular = "ignore",
        optimizer = "bobyqa",
        optCtrl   = list(maxfun = 1e5)
      )
    )
  }
  .varcomps <- function(fit) {
    vc <- lme4::VarCorr(fit)
    comps <- sapply(vc, attr, "stddev")^2
    names(comps) <- names(vc)
    c(comps, Residual = sigma(fit)^2)
  }
  .to_table <- function(v) {
    tot <- sum(v, na.rm = TRUE)
    data.frame(component = names(v),
               variance  = as.numeric(v),
               percent   = as.numeric(100 * v / tot),
               stringsAsFactors = FALSE)
  }
  .fit_circular <- function(dat, short_name) {
    th <- (dat[[short_name]] * pi/180) %% (2*pi)
    tmp <- dat
    tmp$.__sin <- sin(th); tmp$.__cos <- cos(th)
    fm_sin <- full_formula(".__sin")
    fm_cos <- full_formula(".__cos")
    fit_sin <- .fit_lmm(fm_sin, tmp)
    fit_cos <- .fit_lmm(fm_cos, tmp)
    vc <- (.varcomps(fit_sin) + .varcomps(fit_cos)) / 2
    list(varcomps = vc, fits = list(sin = fit_sin, cos = fit_cos),
         formulas = list(sin = fm_sin, cos = fm_cos))
  }

  # fit per variable
  models <- list(); rows_out <- list(); forms <- list(); snr_rows <- list()

  for (short in variables) {
    if (is_circular_short(short)) {
      circ <- .fit_circular(df, short)
      vc   <- circ$varcomps
      tab  <- .to_table(vc)
      r2   <- MuMIn::r.squaredGLMM(circ$fits$sin)
      tab$variable <- short
      tab$R2_c <- as.numeric(r2[2])
      rows_out[[short]] <- tab
      models[[short]]   <- circ
      forms[[short]]    <- circ$formulas
    } else {
      fm  <- full_formula(short)
      fit <- .fit_lmm(fm, df)
      vc  <- .varcomps(fit)
      tab <- .to_table(vc)
      r2  <- MuMIn::r.squaredGLMM(fit)
      tab$variable <- short
      tab$R2_c <- as.numeric(r2[2])
      rows_out[[short]] <- tab
      models[[short]]   <- list(fit = fit, varcomps = vc, formulas = fm)
      forms[[short]]    <- fm
    }

    # SNR: track vs (observer + observer:track + residual)
    st <- rows_out[[short]]
    bio_var   <- sum(st$variance[st$component == "track"], na.rm = TRUE)
    error_var <- sum(st$variance[st$component != "track"], na.rm = TRUE)
    snr_rows[[short]] <- data.frame(
      variable = short,
      bio_component = "track",
      bio_var   = bio_var,
      error_var = error_var,
      SNR       = if (is.finite(error_var) && error_var > 0) bio_var / error_var else NA_real_
    )
  }

  summary <- do.call(rbind, rows_out)
  if (!is.null(summary) && nrow(summary))
    summary <- summary[, c("variable","component","variance","percent","R2_c")]
  snr_tbl <- do.call(rbind, snr_rows)

  # QC table
  qc_rows <- lapply(unique(summary$variable), function(v) {
    tab <- summary[summary$variable == v, , drop = FALSE]
    snr <- snr_tbl$SNR[snr_tbl$variable == v][1]
    rating <- if (is.na(snr)) "NA" else if (snr < 1) "weak (error > bio)" else if (snr < 2) "moderate" else "strong"
    ord <- order(tab$percent, decreasing = TRUE)
    R2c <- unique(tab$R2_c); R2c <- if (length(R2c)) R2c[1] else NA_real_
    obs_est <- if (any(tab$component == "observer")) "yes" else "no (1 level)"
    data.frame(variable = v,
               SNR = snr,
               SNR_rating = rating,
               top_component = tab$component[ord][1],
               top_percent   = tab$percent[ord][1],
               R2_c = R2c,
               observer_estimable = obs_est,
               stringsAsFactors = FALSE)
  })
  qc_tbl <- do.call(rbind, qc_rows)

  res <- list(
    summary        = summary,
    snr            = snr_tbl,
    qc             = qc_tbl,
    models         = models,
    analysis_table = df,
    formulae       = forms,
    track_names    = track_names
  )
  class(res) <- "error_partitioning"
  res
}
