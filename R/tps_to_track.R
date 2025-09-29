#' Transform a *.tps file into a \code{track} R object
#'
#' \code{tps_to_track()} reads a *.tps file containing footprint coordinates of one or several tracks and transforms it into a \code{track} R object.
#'
#' @param file A *.tps file containing (x,y) coordinates of footprints in tracks.
#' @param scale A numeric value specifying the scale in meters per pixel.
#' @param R.L.side A character vector specifying the side of the first footprint of each track. The length of the vector must be equal to the total number of tracks in the sample.
#' \itemize{
#'   \item \code{"L"}: first footprint corresponds to the left foot.
#'   \item \code{"R"}: first footprint corresponds to the right foot.
#' }
#' @param missing A logical value indicating whether there are missing footprints in any track to be interpolated: \code{TRUE}, or \code{FALSE} (the default).
#' @param NAs A matrix with two columns indicating which missing footprints will be interpolated.
#'           The first column gives the number of the track containing missing footprints, and the second column gives the number of the footprint that is missing within this track.
#'           The number of rows is equal to the total number of missing footprints in the sample.
#'
#' @details It is highly recommended that the *.tps file is built using the TPS software (Rohlf 2008, 2009).
#'          Tracks with a different number of footprints (i.e., landmarks) are allowed.
#'          This function transforms the coordinates of the footprints of each track into a set of trajectory coordinates.
#'          Each point of the trajectory is calculated as: \deqn{Point_i(x,y)= (Footprint_i(x,y) + Footprint_{i+1}(x,y)/2}
#'
#' @details The number of points of the resulting trajectory is \eqn{n_{footprints} - 1}.
#'
#' @details If \code{missing} is set to \code{TRUE}, missing footprints can be interpolated.
#'          This interpolation is based on adjacent footprints and the provided side information.
#'
#' @return A \code{track} R object, which is a list consisting of two elements:
#'    * \strong{\code{Trajectories}}: A list of trajectories (midpoints between consecutive footprints).
#'      Includes columns \code{X}, \code{Y}, \code{IMAGE}, \code{ID}, and \code{Side} (set to \code{"Medial"}).
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates and metadata,
#'      with a \code{Side} column (\code{"R"} or \code{"L"}) and a \code{missing} marker (\code{"Actual"} or \code{"Inferred"}).
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
#' # Example 1: Tracks without missing footprints.
#' # Based on the Paluxy River dinosaur chase sequence (Farlow et al., 2011).
#'
#' # Load the example TPS file provided in the QuAnTeTrack package.
#' tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")
#'
#' # Convert the TPS data into a track object.
#' # The 'scale' argument sets the scaling factor,
#' # 'R.L.side' specifies the starting side for each track,
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
#' # Example 2: Tracks with missing footprints.
#' # Based on dinosaur tracks from Mount Tom (Ostrom, 1972).
#'
#' # Load the example TPS file.
#' tpsMountTom <- system.file("extdata", "MountTom.tps", package = "QuAnTeTrack")
#'
#' # Define a matrix representing the missing footprints.
#' # Here, footprint 7 is missing in track 3.
#' NAs <- matrix(c(3, 7), nrow = 1, ncol = 2)
#'
#' # Convert the TPS data into a track object.
#' # The 'missing = TRUE' flag activates interpolation for missing footprints,
#' # 'NAs' specifies which footprints are missing,
#' # and 'R.L.side' indicates the starting side for each track.
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

  ## Errors and Warnings----

  # 'scale' required and must be positive numeric(1)
  if (is.null(scale)) {
    stop("The 'scale' argument is missing. Please provide the scale in pixels per meter.")
  }
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("The 'scale' argument must be a single positive numeric value.")
  }

  # 'missing' must be logical(1)
  if (!is.logical(missing) || length(missing) != 1) {
    stop("The 'missing' argument must be a single logical value: TRUE or FALSE.")
  }

  # Warn if 'NAs' is provided but 'missing' is FALSE
  if (!missing && !is.null(NAs)) {
    warning("The 'NAs' argument will be ignored because 'missing' is set to FALSE.")
  }

  # Validate 'NAs' if provided
  if (!is.null(NAs)) {
    if (!is.matrix(NAs) || ncol(NAs) != 2) {
      stop("The 'NAs' argument must be a matrix with two columns.")
    }
    if (any(NAs <= 0)) {
      stop("The 'NAs' matrix must contain positive integers.")
    }
  }

  # Validate 'R.L.side' basic presence and values
  if (missing(R.L.side) || is.null(R.L.side)) {
    stop("The 'R.L.side' argument is mandatory and must be provided.")
  }
  if (!all(R.L.side %in% c("R", "L"))) {
    stop("The 'R.L.side' vector must contain only 'R' or 'L' values.")
  }

  ## Code----

  # Read file lines
  a <- readLines(file)

  # Identify LM and ID lines
  LM <- grep("LM", a)
  ID.ind <- grep("ID", a)

  # Extract image names
  images <- basename(gsub("(IMAGE=)(.*)", "\\2", a[ID.ind - 1]))

  # Number of rows for each landmark set
  nrows <- as.numeric(gsub("(LM=)([0-9]+)", "\\2", grep("LM", a, value = TRUE)))
  l <- length(LM)

  # Now we can validate length of R.L.side against number of tracks
  if (length(R.L.side) != l) {
    stop(paste0("Length of 'R.L.side' (", length(R.L.side),
                ") must equal the number of tracks detected in the file (", l, ")."))
  }

  # Build landmarks list
  landmarks <- vector("list", l)
  for (i in 1:l) {
    landmarks[i] <- list(data.frame(
      read.table(
        file = file, header = FALSE, skip = LM[i],
        nrows = nrows[i], col.names = c("X", "Y")
      ),
      IMAGE = images[i],
      ID = read.table(
        file = file, header = FALSE, skip = ID.ind[i] - 1,
        nrows = 1, sep = "=", col.names = "ID"
      )[2, ]
    ))
  }

  # Start with landmarks as working data_frame
  data_frame <- landmarks

  ## Assign Side to ALL cases (mandatory R.L.side) ----
  for (i in 1:length(data_frame)) {
    nfp <- nrow(data_frame[[i]])
    idx <- seq_len(nfp)
    side_vec <- character(nfp)
    if (R.L.side[i] == "L") {
      side_vec[idx %% 2 != 0] <- "L"
      side_vec[idx %% 2 == 0] <- "R"
    } else { # R
      side_vec[idx %% 2 != 0] <- "R"
      side_vec[idx %% 2 == 0] <- "L"
    }
    data_frame[[i]]$Side <- side_vec
  }

  if (missing) {
    ## Inferring missing footprints ----

    # Levels from NAs
    levelsnum <- as.numeric(levels(as.factor(NAs[, 1])))

    # Include NAs rows
    for (i in levelsnum) {
      data_frame[[i]] <- berryFunctions::insertRows(
        data_frame[[i]],
        c(NAs[which(NAs[, 1] == i), 2]),
        new = NA
      )
    }

    # Refill IMAGE and ID in inserted rows
    for (i in levelsnum) {
      data_frame[[i]][c(NAs[which(NAs[, 1] == i), 2]), ]$IMAGE <- levels(as.factor(data_frame[[i]]$IMAGE))
      data_frame[[i]][c(NAs[which(NAs[, 1] == i), 2]), ]$ID    <- levels(as.factor(data_frame[[i]]$ID))
    }

    # Recompute Side in case of inserted rows (keep alternation)
    for (i in 1:length(data_frame)) {
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

    ### Calculating track width
    meanwidth <- numeric(length(data_frame))
    for (j in 1:length(data_frame)) {
      vectorwidth <- rep(NA_real_, max(0, nrow(data_frame[[j]]) - 2))
      if (length(vectorwidth) > 0) {
        for (i in 1:(nrow(data_frame[[j]]) - 2)) {
          df <- data_frame[[j]][i:(i + 2), 1:2]
          x1 <- df[1, 1]; x2 <- df[2, 1]; x3 <- df[3, 1]
          y1 <- df[1, 2]; y2 <- df[2, 2]; y3 <- df[3, 2]
          Area <- 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
          Base <- dist(df[c(1, 3), c(1, 2)], method = "euclidean")
          Height <- abs((Area * 2) / Base)
          vectorwidth[i] <- Height
        }
      }
      meanwidth[j] <- mean(vectorwidth[!is.na(vectorwidth)])
    }

    ### Extrapolating missing footprints
    for (i in 1:nrow(NAs)) {
      trk <- NAs[i, 1]
      pos <- NAs[i, 2]
      x1 <- data_frame[[trk]][pos - 1, 1]
      y1 <- data_frame[[trk]][pos - 1, 2]
      x2 <- data_frame[[trk]][pos + 1, 1]
      y2 <- data_frame[[trk]][pos + 1, 2]

      dist_off <- if (data_frame[[trk]][pos, "Side"] == "R") -meanwidth[trk] else meanwidth[trk]

      x3 <- (x1 + x2) / 2
      y3 <- (y1 + y2) / 2

      b <- x2 - x1
      a <- y1 - y2
      norm <- sqrt(a * a + b * b)
      a <- a / norm
      b <- b / norm

      x4 <- x3 + a * dist_off
      y4 <- y3 + b * dist_off

      data_frame[[trk]][pos, 1] <- x4
      data_frame[[trk]][pos, 2] <- y4
    }
  }

  ### Tracing medial tracks ----
  landmarks2 <- data_frame
  landmarks3 <- vector("list", length(landmarks2))

  for (i in 1:length(landmarks2)) {
    # create container with same columns (now includes Side)
    landmarks3[[i]] <- as.data.frame(matrix(ncol = ncol(landmarks2[[i]]),
                                            nrow = (nrow(landmarks2[[i]]) - 1)))
    colnames(landmarks3[[i]]) <- colnames(landmarks2[[i]])

    # carry IMAGE, ID; set Side = "Medial"
    landmarks3[[i]][, "IMAGE"] <- rep(landmarks2[[i]][1, "IMAGE"], nrow(landmarks3[[i]]))
    landmarks3[[i]][, "ID"]    <- rep(landmarks2[[i]][1, "ID"],    nrow(landmarks3[[i]]))
    landmarks3[[i]][, "Side"]  <- rep("Medial",                    nrow(landmarks3[[i]]))

    # Keep the first five columns: X, Y, IMAGE, ID, Side
    landmarks3[[i]] <- landmarks3[[i]][, c("X", "Y", "IMAGE", "ID", "Side")]

    # midpoint X, Y
    for (j in 1:nrow(landmarks3[[i]])) {
      landmarks3[[i]][j, "X"] <- (landmarks2[[i]][j, "X"] + landmarks2[[i]][j + 1, "X"]) / 2
      landmarks3[[i]][j, "Y"] <- (landmarks2[[i]][j, "Y"] + landmarks2[[i]][j + 1, "Y"]) / 2
    }
  }

  # Scale and convert to trajectories
  for (i in 1:length(landmarks3)) {
    landmarks3[[i]][, c("X", "Y")] <- landmarks3[[i]][, c("X", "Y")] * scale
    landmarks3[[i]] <- TrajFromCoords(landmarks3[[i]])
  }
  names(landmarks3) <- paste0("Track_", str_pad(1:length(LM), nchar(length(LM)), pad = "0"))

  # Mark missing footprints
  for (i in 1:length(data_frame)) {
    data_frame[[i]]$missing <- rep("Actual", nrow(data_frame[[i]]))
  }
  if (missing) {
    levelsnum <- as.numeric(levels(as.factor(NAs[, 1])))
    for (i in levelsnum) {
      data_frame[[i]][c(NAs[which(NAs[, 1] == i), 2]), ]$missing <- "Inferred"
    }
  }

  # Scale footprints X, Y
  for (i in 1:length(data_frame)) {
    data_frame[[i]]$X <- data_frame[[i]]$X * scale
    data_frame[[i]]$Y <- data_frame[[i]]$Y * scale
  }

  # Return
  list(
    Trajectories = landmarks3,
    Footprints = data_frame
  )
}
