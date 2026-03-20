#' Transform a *.tps file into a \code{trackway} R object
#'
#' \code{tps_to_track()} reads a *.tps file containing footprint coordinates of one or several trackways and transforms it into a \code{trackway} R object.
#'
#' @param file A *.tps file containing (x,y) coordinates of footprints in trackways.
#' @param scale A numeric value specifying the scale in meters per pixel.
#' @param R.L.side A character vector specifying the laterality of the first footprint of each trackway. The length of the vector must be equal to the total number of trackways in the sample.
#' \itemize{
#'   \item \code{"L"}: first footprint corresponds to the left foot.
#'   \item \code{"R"}: first footprint corresponds to the right foot.
#' }
#' @param missing A logical value indicating whether there are missing footprints in any trackway to be interpolated: \code{TRUE}, or \code{FALSE} (the default).
#' @param NAs A matrix with two columns indicating which missing footprints will be interpolated.
#'           The first column gives the number of the trackway containing missing footprints, and the second column gives the number of the footprint that is missing within this trackway.
#'           The number of rows is equal to the total number of missing footprints in the sample.
#'
#' @details It is highly recommended that the *.tps file is built using the TPS software (Rohlf 2008, 2009).
#'          Trackways with a different number of footprints (i.e., landmarks) are allowed.
#'          This function transforms the coordinates of the footprints of each trackway into a set of trajectory coordinates.
#'          Each point of the trajectory is calculated as: \deqn{Point_i(x,y)= (Footprint_i(x,y) + Footprint_{i+1}(x,y)/2}
#'
#' @details The number of points of the resulting trajectory is \eqn{n_{footprints} - 1}.
#'
#' @details If \code{missing} is set to \code{TRUE}, missing footprints can be interpolated.
#'           The interpolated footprint is then placed at the midpoint between the two adjacent
#'          footprints and shifted laterally along the direction perpendicular to their connecting
#'          segment. The magnitude of this lateral shift is estimated from the mean trackway width
#'          of the corresponding trackway, and its sign is determined by the inferred anatomical
#'          side (left or right) of the missing footprint. Inferred footprints are flagged as such
#'          in the output.
#'
#' @return A \code{trackway} R object, which is a list consisting of two elements:
#'
#'    * \strong{\code{Trajectories}}: A list of trajectories representing trackway midlines,
#'      interpolated by connecting the midpoints of successive left–right footprint pairs
#'      (i.e., footprints linked by pace lines). Includes columns \code{X}, \code{Y},
#'      \code{IMAGE}, \code{ID}, and \code{Side} (set to \code{"Medial"}).
#'
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates
#'      and associated metadata, with a \code{Side} column (\code{"R"} or \code{"L"})
#'      and a \code{missing} marker (\code{"Actual"} or \code{"Inferred"}).
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
#' # Example 1: Trackways without missing footprints.
#' # Based on the Paluxy River dinosaur chase sequence (Farlow et al., 2011).
#'
#' # Load the example TPS file provided in the QuAnTeTrack package.
#' tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")
#'
#' # Convert the TPS data into a trackway object.
#' # The 'scale' argument sets the scaling factor,
#' # 'R.L.side' specifies the starting side for each trackway,
#' # and 'missing = FALSE' indicates no footprints are missing.
#' tps_to_track(
#'   tpsPaluxyRiver,
#'   scale = 0.004341493,
#'   R.L.side = c("R", "L"),
#'   missing = FALSE,
#'   NAs = NULL
#' )
#'
#'
#' # Example 2: Trackways with missing footprints.
#' # Based on dinosaur trackways from Mount Tom (Ostrom, 1972).
#'
#' # Load the example TPS file.
#' tpsMountTom <- system.file("extdata", "MountTom.tps", package = "QuAnTeTrack")
#'
#' # Define a matrix representing the missing footprints.
#' # Here, footprint 7 is missing in trackway 3.
#' NAs <- matrix(c(3, 7), nrow = 1, ncol = 2)
#'
#' # Convert the TPS data into a trackway object.
#' # The 'missing = TRUE' flag activates interpolation for missing footprints,
#' # 'NAs' specifies which footprints are missing,
#' # and 'R.L.side' indicates the starting side for each trackway.
#' tps_to_track(
#'   tpsMountTom,
#'   scale = 0.004411765,
#'   R.L.side = c(
#'     "R", "L", "L", "L", "R", "L", "R", "R", "L", "L", "L",
#'     "L", "L", "R", "R", "L", "R", "R", "L", "R", "R",
#'     "R", "R"
#'   ),
#'   missing = TRUE,
#'   NAs = NAs
#' )
#'
#' @references
#' Farlow, J. O., O’Brien, M., Kuban, G. J., Dattilo, B. F., Bates, K. T., Falkingham, P. L., & Piñuela, L. (2012). Dinosaur Tracksites of the Paluxy River Valley (Glen Rose Formation, Lower Cretaceous), Dinosaur Valley State Park, Somervell County, Texas. In Proceedings of the V International Symposium about Dinosaur Palaeontology and their Environment (pp. 41-69). Burgos: Salas de los Infantes.
#'
#' Ostrom, J. H. (1972). Were some dinosaurs gregarious?. Palaeogeography, Palaeoclimatology, Palaeoecology, 11(4), 287-301.
#'
#' Rohlf, F. J. 2008. TPSUTIL. Version 1.40. Department of Ecology and Evolution, State University of New York.
#'          Available at <https://sbmorphometrics.org/>.
#'
#' Rohlf, F. J. 2009. tpsDig. Version 2.14. Department of Ecology and Evolution, State University of New York.
#'          Available at <https://sbmorphometrics.org/>.
#'
#' @importFrom berryFunctions insertRows
#' @importFrom dplyr bind_rows
#' @importFrom geomorph digit.curves
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 theme_light
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggrepel geom_label_repel
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics legend
#' @importFrom grDevices chull
#' @importFrom grDevices colors
#' @importFrom stats cor
#' @importFrom stringr str_pad
#' @importFrom trajr TrajFromCoords
#' @importFrom utils read.table
#'
#' @export


tps_to_track <- function(file, scale = NULL, R.L.side, missing = FALSE, NAs = NULL) {

  ## ---------- Errors and Warnings ----------
  if (is.null(scale)) {
    stop("The 'scale' argument is missing. Please provide the scale in meters per pixel.")
  }
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("The 'scale' argument must be a single positive numeric value.")
  }
  if (!is.logical(missing) || length(missing) != 1) {
    stop("The 'missing' argument must be a single logical value: TRUE or FALSE.")
  }
  if (!missing && !is.null(NAs)) {
    warning("The 'NAs' argument will be ignored because 'missing' is set to FALSE.")
  }
  if (!is.null(NAs)) {
    if (!is.matrix(NAs) || ncol(NAs) != 2) {
      stop("The 'NAs' argument must be a matrix with two columns.")
    }
    if (any(NAs <= 0)) {
      stop("The 'NAs' matrix must contain positive integers.")
    }
  }
  if (missing(R.L.side) || is.null(R.L.side)) {
    stop("The 'R.L.side' argument is mandatory and must be provided.")
  }
  if (!all(R.L.side %in% c("R", "L"))) {
    stop("The 'R.L.side' vector must contain only 'R' or 'L' values.")
  }

  ## ---------- Read TPS ----------
  a <- readLines(file, warn = FALSE)

  # Indices of key lines
  LM     <- grep("^\\s*LM\\s*=\\s*", a)
  ID.ind <- grep("^\\s*ID\\s*=\\s*", a)

  if (length(LM) == 0 || length(ID.ind) == 0) {
    stop("Could not find required LM= or ID= lines in the TPS file.")
  }

  # Extract IMAGE from the line immediately before each ID line
  # (common TPS convention: ... coordinates ... IMAGE=...  ID=...)
  img_lines <- a[ID.ind - 1]
  images <- trimws(sub("^.*?IMAGE\\s*=\\s*", "", img_lines), which = "both")

  # Extract ID values from the ID lines
  ids <- trimws(sub("^\\s*ID\\s*=\\s*", "", a[ID.ind]), which = "both")

  # Number of rows (landmarks) for each block
  nrows <- as.numeric(sub("^\\s*LM\\s*=\\s*([0-9]+).*", "\\1", a[LM]))
  l <- length(LM)

  # Validate side vector length
  if (length(R.L.side) != l) {
    stop(sprintf("Length of 'R.L.side' (%d) must equal the number of trackways detected (%d).",
                 length(R.L.side), l))
  }

  # ---------- NEW: Drop an error if IMAGE or ID repeats ----------
  dup_img <- unique(images[duplicated(images)])
  dup_id  <- unique(ids[duplicated(ids)])

  if (length(dup_img) > 0) {
    stop(sprintf(
      "Duplicate IMAGE names detected in TPS. IMAGE= must be unique per trackway.\nRepeated: %s",
      paste(dup_img, collapse = ", ")
    ))
  }
  if (length(dup_id) > 0) {
    stop(sprintf(
      "Duplicate ID values detected in TPS. ID= must be unique per trackway.\nRepeated: %s",
      paste(dup_id, collapse = ", ")
    ))
  }

  ## ---------- Build landmarks list ----------
  landmarks <- vector("list", l)
  for (i in seq_len(l)) {
    # read coordinates in block i
    coords <- utils::read.table(file = file, header = FALSE, skip = LM[i],
                                nrows = nrows[i], col.names = c("X", "Y"))
    # attach metadata from parsed IMAGE/ID (single source of truth)
    coords$IMAGE <- images[i]
    coords$ID    <- ids[i]
    landmarks[[i]] <- coords
  }

  data_frame <- landmarks

  ## ---------- Assign Side (mandatory R.L.side) ----------
  for (i in seq_along(data_frame)) {
    nfp <- nrow(data_frame[[i]])
    idx <- seq_len(nfp)
    side_vec <- character(nfp)
    if (R.L.side[i] == "L") {
      side_vec[idx %% 2 != 0] <- "L"
      side_vec[idx %% 2 == 0] <- "R"
    } else {
      side_vec[idx %% 2 != 0] <- "R"
      side_vec[idx %% 2 == 0] <- "L"
    }
    data_frame[[i]]$Side <- side_vec
  }

  ## ---------- Handle missing footprints (optional) ----------
  if (missing) {
    levelsnum <- as.numeric(levels(as.factor(NAs[, 1])))

    # Insert NA rows for missing positions
    for (i in levelsnum) {
      data_frame[[i]] <- berryFunctions::insertRows(
        data_frame[[i]],
        c(NAs[which(NAs[, 1] == i), 2]),
        new = NA
      )
    }

    # Refill IMAGE and ID in inserted rows
    for (i in levelsnum) {
      idx_na <- c(NAs[which(NAs[, 1] == i), 2])
      data_frame[[i]][idx_na, "IMAGE"] <- data_frame[[i]]$IMAGE[1]
      data_frame[[i]][idx_na, "ID"]    <- data_frame[[i]]$ID[1]
    }

    # Recompute Side (keep alternation starting from R.L.side[i])
    for (i in seq_along(data_frame)) {
      nfp <- nrow(data_frame[[i]])
      idx <- seq_len(nfp)
      side_vec <- character(nfp)
      if (R.L.side[i] == "L") {
        side_vec[idx %% 2 != 0] <- "L"
        side_vec[idx %% 2 == 0] <- "R"
      } else {
        side_vec[idx %% 2 != 0] <- "R"
        side_vec[idx %% 2 == 0] <- "L"
      }
      data_frame[[i]]$Side <- side_vec
    }

    # --- Trackway width estimation
    meanwidth <- numeric(length(data_frame))
    for (j in seq_along(data_frame)) {
      if (nrow(data_frame[[j]]) > 2) {
        vw <- rep(NA_real_, nrow(data_frame[[j]]) - 2)
        for (k in seq_len(length(vw))) {
          df <- data_frame[[j]][k:(k + 2), 1:2]
          x1 <- df[1, 1]; x2 <- df[2, 1]; x3 <- df[3, 1]
          y1 <- df[1, 2]; y2 <- df[2, 2]; y3 <- df[3, 2]
          Area <- 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
          Base <- dist(df[c(1, 3), c(1, 2)], method = "euclidean")
          vw[k] <- abs((Area * 2) / Base)
        }
        meanwidth[j] <- mean(vw, na.rm = TRUE)
      } else {
        meanwidth[j] <- NA_real_
      }
    }

    # --- Extrapolate missing footprints
    for (r in seq_len(nrow(NAs))) {
      trk <- NAs[r, 1]
      pos <- NAs[r, 2]
      x1 <- data_frame[[trk]][pos - 1, 1]; y1 <- data_frame[[trk]][pos - 1, 2]
      x2 <- data_frame[[trk]][pos + 1, 1]; y2 <- data_frame[[trk]][pos + 1, 2]
      dist_off <- if (data_frame[[trk]][pos, "Side"] == "R") -meanwidth[trk] else meanwidth[trk]

      x3 <- (x1 + x2) / 2; y3 <- (y1 + y2) / 2
      b <- x2 - x1; a <- y1 - y2
      norm <- sqrt(a * a + b * b); a <- a / norm; b <- b / norm
      x4 <- x3 + a * dist_off; y4 <- y3 + b * dist_off

      data_frame[[trk]][pos, 1] <- x4
      data_frame[[trk]][pos, 2] <- y4
    }
  }

  ## ---------- Medial trajectories (midpoints) ----------
  landmarks2 <- data_frame
  landmarks3 <- vector("list", length(landmarks2))

  for (i in seq_along(landmarks2)) {
    if (nrow(landmarks2[[i]]) < 2) {
      stop(sprintf("Trackway %d has fewer than 2 footprints after parsing; cannot compute trajectory.", i))
    }

    out <- as.data.frame(matrix(ncol = ncol(landmarks2[[i]]),
                                nrow = (nrow(landmarks2[[i]]) - 1)))
    colnames(out) <- colnames(landmarks2[[i]])

    out[, "IMAGE"] <- rep(landmarks2[[i]][1, "IMAGE"], nrow(out))
    out[, "ID"]    <- rep(landmarks2[[i]][1, "ID"],    nrow(out))
    out[, "Side"]  <- rep("Medial",                    nrow(out))

    out <- out[, c("X", "Y", "IMAGE", "ID", "Side")]

    for (j in seq_len(nrow(out))) {
      out[j, "X"] <- (landmarks2[[i]][j, "X"] + landmarks2[[i]][j + 1, "X"]) / 2
      out[j, "Y"] <- (landmarks2[[i]][j, "Y"] + landmarks2[[i]][j + 1, "Y"]) / 2
    }

    # scale midpoints
    out[, c("X", "Y")] <- out[, c("X", "Y")] * scale
    # convert to trajr object
    out <- trajr::TrajFromCoords(out)

    landmarks3[[i]] <- out
  }
  names(landmarks3) <- paste0("Trackway_", stringr::str_pad(seq_along(LM), nchar(length(LM)), pad = "0"))

  ## ---------- Mark missing & scale footprints ----------
  for (i in seq_along(data_frame)) {
    data_frame[[i]]$missing <- "Actual"
  }
  if (missing) {
    levelsnum <- as.numeric(levels(as.factor(NAs[, 1])))
    for (i in levelsnum) {
      idx_na <- c(NAs[which(NAs[, 1] == i), 2])
      data_frame[[i]][idx_na, "missing"] <- "Inferred"
    }
  }
  for (i in seq_along(data_frame)) {
    data_frame[[i]]$X <- data_frame[[i]]$X * scale
    data_frame[[i]]$Y <- data_frame[[i]]$Y * scale
  }

  ## ---------- Return ----------
  list(
    Trajectories = landmarks3,
    Footprints   = data_frame
  )
}
