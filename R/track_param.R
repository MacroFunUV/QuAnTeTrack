#' Print track parameters
#'
#' \code{track_param()} is a function to compute and print various parameters of trackways.
#'
#' @param data A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#' @param gauge_size Numeric. Pes/manus length (or width) used to compute Gauge as \code{Trackway_width / gauge_size}.
#' If \code{NULL} or \code{NA}, Gauge is returned as \code{NA}.
#'
#' @details This function calculates various movement parameters for each track in the provided data.
#'
#' It uses the following helper functions:
#'
#' From the \pkg{trajr} package:
#'    * \code{TrajAngles()}: Calculates the turning angles of the track.
#'    * \code{TrajDistance()}: Calculates the total distance covered by the track.
#'    * \code{TrajLength()}: Calculates the length of the track.
#'    * \code{TrajStepLengths()}: Calculates the step lengths of the track.
#'    * \code{TrajSinuosity2()}: Calculates the sinuosity of the track.
#'    * \code{TrajStraightness()}: Calculates the straightness of the track.
#'
#' From the \pkg{circular} package:
#'    * \code{circular()}: Converts raw angles (in radians) into a circular data type.
#'    * \code{mean.circular()}: Computes the circular mean of the turning angles.
#'    * \code{sd.circular()}: Computes the circular standard deviation of the turning angles.
#'
#' The circular mean and circular standard deviation are returned in degrees in this function.
#'
#' The reference direction, or 0 degrees, is considered to be along the positive x-axis. This means that angles are measured counterclockwise from the positive x-axis, with 0 degrees pointing directly along this axis. For a detailed explanation and appropriate methods for analyzing circular data, refer to Batschelet (1981).
#'
#' Circular mean of turning angles is computed as:
#'
#' \deqn{\overline{\theta} = atan2\left(\frac{1}{n}\sum_{i=1}^n \sin \theta_i, \frac{1}{n}\sum_{i=1}^n \cos \theta_i\right)}
#'
#' where:
#' \describe{
#'   \item{\eqn{\theta_i}}{ is the \eqn{i^{th}} turning angle in radians.}
#'   \item{\eqn{n}}{ is the total number of turning angles.}
#'   \item{\eqn{\sin \theta_i}, \eqn{\cos \theta_i}}{ are the sine and cosine components of each turning angle.}
#' \item{\eqn{\mathrm{atan2}(y,x)}}{ is the two-argument arctangent that returns the angle in the correct quadrant.}
#' }
#'
#' Circular standard deviation of turning angles is computed as:
#'
#' \deqn{s_c = \sqrt{-2 \ln(\overline{R})}, \quad \overline{R} = \sqrt{\left(\frac{1}{n}\sum_{i=1}^n \cos \theta_i\right)^2 + \left(\frac{1}{n}\sum_{i=1}^n \sin \theta_i\right)^2}}
#'
#' where:
#' \describe{
#'   \item{\eqn{\overline{R}}}{ is the mean resultant length, measuring concentration of angles around the mean direction.}
#'   \item{\eqn{n}}{ is the total number of turning angles.}
#'   \item{\eqn{\cos \theta_i}, \eqn{\sin \theta_i}}{ are the trigonometric components of each angle.}
#'   \item{\eqn{s_c}}{ is the circular standard deviation in radians (converted to degrees in this function).}
#' }
#'
#' Sinuosity is calculated according to Benhamou (2004), as defined in equation 8.
#' The formula used here is a refined version of the sinuosity index presented by Bovet & Benhamou (1988),
#' which is applicable to a broader range of turning angle distributions and does not require a constant step length.
#'
#' The sinuosity is computed using the formula:
#' \deqn{S = 2 \left[ p \left( \frac{1 + c}{1 - c} + b^2 \right) \right]^{-0.5}}
#' where:
#' \describe{
#' \item{\eqn{p}}{ is the mean step length (in meters),}
#' \item{\eqn{c}}{ is the mean cosine of turning angles (in radians), and}
#' \item{\eqn{b}}{ is the coefficient of variation of the step length (in meters).}
#' }
#'
#' The straightness index is defined as the ratio D/L, where:
#' \describe{
#' \item{\eqn{D}}{ is the beeline distance between the first and last points in the trajectory (in meters), and}
#' \item{\eqn{L}}{ is the total path length traveled (in meters).}
#' }
#'
#' Straightness index is based on the method described by Batschelet (1981). According to Benhamou (2004),
#' the straightness index serves as a reliable measure of the efficiency of a directed walk. However, it is not suitable
#' for random trajectories, as the index for a random walk tends towards zero with increasing steps. Thus, it is recommended
#' to use this measure to compare the tortuosity of random walks only if they consist of a similar number of steps.
#'
#' @return A list of lists, where each sublist contains the computed parameters for a corresponding track.
#' The parameters included are:
#' \describe{
#'   \item{\code{Turning_angles}}{A vector of turning angles for the track (in degrees).}
#'   \item{\code{Mean_turning_angle}}{The mean of the turning angles (in degrees).}
#'   \item{\code{Standard_deviation_turning_angle}}{The standard deviation of the turning angles (in degrees).}
#'
#'   \item{\code{Path_length}}{Total path length (in meters), computed as the sum of distances between consecutive trajectory points.}
#'   \item{\code{Beeline_length}}{Straight-line distance (in meters) between the first and last points of the trajectory.}
#'
#'   \item{\code{Step_lengths}}{A vector of step lengths for the track  (in meters).}
#'   \item{\code{Mean_step_length}}{The mean of the step lengths  (in meters).}
#'   \item{\code{Standard_deviation_step_length}}{The standard deviation of the step lengths  (in meters).}
#'
#'   \item{\code{Stride_length}}{A vector of stride lengths (in meters), computed as distances between consecutive footprints of the same side (L–L and R–R), in footprint order.}
#'   \item{\code{Mean_stride_length}}{Mean stride length (in meters). Computed as \code{Mean_step_length * 2}.}
#'
#'   \item{\code{Pace_length}}{A vector of pace lengths (in meters), computed as distances between consecutive contralateral footprints (L–R or R–L) in footprint order.}
#'   \item{\code{Mean_pace_length}}{Mean pace length (in meters).}
#'
#'   \item{\code{Sinuosity}}{The sinuosity of the track (dimensionless).}
#'   \item{\code{Straightness}}{The straightness of the track (dimensionless).}
#'
#'   \item{\code{Trackway_width}}{Mean lateral separation between left and right footprints (in meters), measured perpendicular to the inferred trackway axis.}
#'   \item{\code{Gauge}}{Trackway gauge (dimensionless), computed as \code{Trackway_width / gauge_size}.}
#'   \item{\code{Pace_angulation}}{Mean interior angle (in degrees) computed from alternating triplets (L–R–L or R–L–R).}
#'   \item{\code{Step_angle}}{Mean step angle (in degrees), computed as the mean angle between each pace segment and the inferred stride/trackway axis.}
#' }
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
#'
#' Batschelet, E. (1981). Circular statistics in biology. Academic press, 111 Fifth Ave., New York, NY 10003, 1981, 388.
#'
#' Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's path:: straightness, sinuosity, or fractal dimension?. Journal of theoretical biology, 229(2), 209-220.
#'
#' Bovet, P., & Benhamou, S. (1988). Spatial analysis of animals' movements using a correlated random walk model. Journal of theoretical biology, 131(4), 419-433.
#'
#' @examples
#' # Example 1:
#' track_param(PaluxyRiver)
#'
#' # Example 2:
#' track_param(MountTom)
#'
#' @importFrom trajr TrajAngles
#' @importFrom trajr TrajDistance
#' @importFrom trajr TrajLength
#' @importFrom trajr TrajStepLengths
#' @importFrom trajr TrajSinuosity2
#' @importFrom trajr TrajStraightness
#' @importFrom circular circular mean.circular sd.circular
#' @importFrom stats prcomp
#'
#' @seealso \code{\link{tps_to_track}}
#'
#' @export
track_param <- function(data, gauge_size = NA) {
  ## Errors and Warnings----

  # Check structure
  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'track' R object.")
  }
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("Both elements of 'data' must be lists. Ensure that 'Trajectories' and 'Footprints' are provided.")
  }

  ## Code----
  trajectories <- data[[1]]
  footprints_list <- data[[2]]

  list <- list()

  # helper: mean circular angle at middle of alternating triplets
  .pace_angulation <- function(fp) {
    if (!all(c("X","Y","Side") %in% names(fp))) return(NA_real_)
    side <- as.character(fp$Side)
    # normalize to "L"/"R"
    side <- ifelse(grepl("^[lL]", side), "L",
                   ifelse(grepl("^[rR]", side), "R", NA))
    ok <- !is.na(fp$X) & !is.na(fp$Y) & !is.na(side)
    if (sum(ok) < 3) return(NA_real_)
    X <- fp$X[ok]; Y <- fp$Y[ok]; s <- side[ok]
    n <- length(s)
    angs <- numeric(0)
    for (k in 2:(n-1)) {
      if (s[k-1] != s[k] && s[k] != s[k+1] && s[k-1] == s[k+1]) {
        p1 <- c(X[k-1], Y[k-1]); p2 <- c(X[k], Y[k]); p3 <- c(X[k+1], Y[k+1])
        u <- p1 - p2; v <- p3 - p2
        nu <- sqrt(sum(u*u)); nv <- sqrt(sum(v*v))
        if (nu > 0 && nv > 0) {
          cosang <- sum(u*v)/(nu*nv)
          cosang <- max(min(cosang, 1), -1)
          angs <- c(angs, acos(cosang) * 180/pi)
        }
      }
    }
    if (length(angs) == 0) NA_real_ else mean(angs)
  }

  # helper: trackway width using trajectory axis + perpendicular distances of footprints
  .trackway_width <- function(traj, fp) {
    if (!all(c("x","y") %in% names(traj))) return(NA_real_)
    if (!all(c("X","Y","Side") %in% names(fp))) return(NA_real_)
    side <- as.character(fp$Side)
    side <- ifelse(grepl("^[lL]", side), "L",
                   ifelse(grepl("^[rR]", side), "R", NA))
    ok_fp <- !is.na(fp$X) & !is.na(fp$Y) & !is.na(side)
    if (sum(ok_fp) < 2 || length(unique(side[ok_fp])) < 2) return(NA_real_)

    coords_traj <- cbind(traj$x, traj$y)
    if (nrow(coords_traj) < 2) return(NA_real_)
    pc <- prcomp(coords_traj, center = TRUE, scale. = FALSE)
    dir <- pc$rotation[, 1]; dir <- dir / sqrt(sum(dir^2))
    ctr <- colMeans(coords_traj)

    # perpendicular unit vector to axis
    nvec <- c(-dir[2], dir[1]); nvec <- nvec / sqrt(sum(nvec^2))

    # signed perpendicular distances of footprints to axis
    coords_fp <- cbind(fp$X[ok_fp], fp$Y[ok_fp])
    d <- as.numeric((coords_fp[,1] - ctr[1]) * nvec[1] + (coords_fp[,2] - ctr[2]) * nvec[2])
    s_ok <- side[ok_fp]
    Ld <- d[s_ok == "L"]; Rd <- d[s_ok == "R"]
    if (length(Ld) == 0 || length(Rd) == 0) return(NA_real_)
    abs(mean(Ld, na.rm = TRUE) - mean(Rd, na.rm = TRUE))
  }

  # helper: stride lengths (same side: L-L and R-R) from footprints
  .stride_length <- function(fp) {
    if (!all(c("X","Y","Side") %in% names(fp))) return(numeric(0))
    side <- as.character(fp$Side)
    side <- ifelse(grepl("^[lL]", side), "L",
                   ifelse(grepl("^[rR]", side), "R", NA))
    ok <- !is.na(fp$X) & !is.na(fp$Y) & !is.na(side)
    if (sum(ok) < 2) return(numeric(0))
    X <- fp$X[ok]; Y <- fp$Y[ok]; s <- side[ok]

    out <- numeric(0)
    for (ss in c("L","R")) {
      idx <- which(s == ss)
      if (length(idx) >= 2) {
        dx <- diff(X[idx]); dy <- diff(Y[idx])
        out <- c(out, sqrt(dx*dx + dy*dy))
      }
    }
    out
  }

  # helper: pace lengths (contralateral consecutive: L-R or R-L) from footprints
  .pace_length <- function(fp) {
    if (!all(c("X","Y","Side") %in% names(fp))) return(numeric(0))
    side <- as.character(fp$Side)
    side <- ifelse(grepl("^[lL]", side), "L",
                   ifelse(grepl("^[rR]", side), "R", NA))
    ok <- !is.na(fp$X) & !is.na(fp$Y) & !is.na(side)
    if (sum(ok) < 2) return(numeric(0))
    X <- fp$X[ok]; Y <- fp$Y[ok]; s <- side[ok]

    out <- numeric(0)
    n <- length(s)
    for (k in 2:n) {
      if (!is.na(s[k]) && !is.na(s[k-1]) && s[k] != s[k-1]) {
        dx <- X[k] - X[k-1]; dy <- Y[k] - Y[k-1]
        out <- c(out, sqrt(dx*dx + dy*dy))
      }
    }
    out
  }

  # helper: step angle (angle between pace segments and inferred stride axis)
  .step_angle <- function(traj, fp) {
    if (!all(c("x","y") %in% names(traj))) return(NA_real_)
    if (!all(c("X","Y","Side") %in% names(fp))) return(NA_real_)

    side <- as.character(fp$Side)
    side <- ifelse(grepl("^[lL]", side), "L",
                   ifelse(grepl("^[rR]", side), "R", NA))
    ok <- !is.na(fp$X) & !is.na(fp$Y) & !is.na(side)
    if (sum(ok) < 2) return(NA_real_)

    # stride axis from trajectory PCA (same approach as trackway width)
    coords_traj <- cbind(traj$x, traj$y)
    if (nrow(coords_traj) < 2) return(NA_real_)
    pc <- prcomp(coords_traj, center = TRUE, scale. = FALSE)
    dir <- pc$rotation[, 1]
    dir <- dir / sqrt(sum(dir^2))

    X <- fp$X[ok]; Y <- fp$Y[ok]; s <- side[ok]
    n <- length(s)

    angs <- numeric(0)
    for (k in 2:n) {
      if (s[k] != s[k-1]) {
        v <- c(X[k] - X[k-1], Y[k] - Y[k-1])
        nv <- sqrt(sum(v*v))
        if (nv > 0) {
          cosang <- sum(v * dir) / nv
          cosang <- max(min(abs(cosang), 1), -1) # ignore direction sign
          angs <- c(angs, acos(cosang) * 180/pi)
        }
      }
    }
    if (length(angs) == 0) NA_real_ else mean(angs)
  }

  for (i in seq_along(trajectories)) {

    # Turning angles (degrees) and circular summaries
    angles_rad <- TrajAngles(trajectories[[i]], compass.direction = 0)
    Turning_angles <- angles_rad * (180 / pi)

    ang_circ <- circular(angles_rad, units = "radians", modulo = "2pi")
    Mean_turning_angle <- as.numeric(mean(ang_circ)) * (180 / pi)
    Standard_deviation_turning_angle <- as.numeric(sd(ang_circ)) * (180 / pi)

    # Distance/length/steps/sinuosity/straightness
    Beeline_length <- TrajDistance(trajectories[[i]])
    Path_length <- TrajLength(trajectories[[i]])
    Step_lengths <- TrajStepLengths(trajectories[[i]])
    Mean_step_length <- mean(Step_lengths)
    Standard_deviation_step_length <- sd(Step_lengths)
    Sinuosity <- TrajSinuosity2(trajectories[[i]])
    Straightness <- TrajStraightness(trajectories[[i]])

    # Derived footprint-based metrics (X,Y, Side)
    fp_i <- footprints_list[[i]]
    Trackway_width <- Pace_angulation <- NA_real_
    Stride_length <- Pace_length <- numeric(0)
    Mean_pace_length <- Gauge <- Step_angle <- NA_real_

    if (is.data.frame(fp_i)) {
      if (!all(c("X","Y") %in% names(fp_i))) {
        warning(sprintf("Track %d: Footprints missing 'X'/'Y'; derived footprint metrics set to NA.", i))
      } else if (!("Side" %in% names(fp_i))) {
        warning(sprintf("Track %d: Footprints missing 'Side'; derived footprint metrics set to NA.", i))
      } else {
        Trackway_width <- .trackway_width(trajectories[[i]], fp_i)
        Pace_angulation <- .pace_angulation(fp_i)

        Stride_length <- .stride_length(fp_i)
        Mean_stride_length <- Mean_step_length * 2

        Pace_length <- .pace_length(fp_i)
        Mean_pace_length <- if (length(Pace_length) == 0) NA_real_ else mean(Pace_length)

        Gauge <- if (is.null(gauge_size) || is.na(gauge_size)) NA_real_ else (Trackway_width / gauge_size)

        Step_angle <- .step_angle(trajectories[[i]], fp_i)
      }
    } else {
      warning(sprintf("Track %d: Footprints element is not a data.frame; derived footprint metrics set to NA.", i))
      Mean_stride_length <- Mean_step_length * 2
    }

    # Coherent, logical output order
    sublist <- list(
      Turning_angles = Turning_angles,
      Mean_turning_angle = Mean_turning_angle,
      Standard_deviation_turning_angle = Standard_deviation_turning_angle,

      Beeline_length = Beeline_length,
      Path_length = Path_length,

      Step_lengths = Step_lengths,
      Mean_step_length = Mean_step_length,
      Standard_deviation_step_length = Standard_deviation_step_length,

      Stride_length = Stride_length,
      Mean_stride_length = Mean_step_length * 2,

      Pace_length = Pace_length,
      Mean_pace_length = Mean_pace_length,

      Sinuosity = Sinuosity,
      Straightness = Straightness,

      Trackway_width = Trackway_width,
      Gauge = Gauge,
      Pace_angulation = Pace_angulation,
      Step_angle = Step_angle
    )

    list[[i]] <- sublist
  }

  return(list)
}
