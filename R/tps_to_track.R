#' Transform a *.tps file into a 'track' R object
#'
#' \code{tps_to_track()} reads a *.tps file containing footprint coordinates of one or several tracks and transforms it into a 'track' R object.
#'
#' @param file A *.tps file containing [x,y] coordinates of footprints in tracks.
#' @param scale A numeric value specifying the scale in meters per pixel.
#' @param missing A logical value indicating whether there are missing footprints in any track to be interpolated: \code{TRUE}, or \code{FALSE} (the default).
#' @param NAs A matrix with two columns indicating which missing footprints will be interpolated.
#'           The first column gives the number of the track containing missing footprints, and the second column gives the number of the footprint that is missing within this track.
#'           The number of rows is equal to the total number of missing footprints in the sample.
#' @param R.L.side A character vector specifying the side of the first footprint of each track.
#'           Only needed if \code{missing} is set to \code{TRUE}.
#'           The length of the vector must be equal to the total number of tracks in the sample.
#'    * \code{"L"}: first footprint corresponds to the left foot.
#'    * \code{"R"}: first footprint corresponds to the right foot.
#'
#' @details It is highly recommended that the *.tps file is built using the TPS software (Rohlf 2008, 2009).
#'          Tracks with a different number of footprints (i.e., landmarks) are allowed.
#'          This function transforms the coordinates of the footprints of each track into a set of trajectory coordinates.
#'          Each point of the trajectory is calculated as: \deqn{Point_i[x,y]= (Footprint_i[x,y] + Footprint_{i+1}[x,y])/2}
#'
#' @details The number of points of the resulting trajectory is \eqn{n_{footprints} - 1}.
#'
#' @details If \code{missing} is set to \code{TRUE}, missing footprints can be interpolated.
#'          This interpolation is based on adjacent footprints and the provided side information.
#'
#' @return A 'track' R object, which is a list consisting of two elements:
#'    * \strong{Trajectories}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
#'    * \strong{Footprints}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
#'
#' @section Logo:
#'\if{html}{\figure{Logo.png}{options: width=30\%}}
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
#' @examples
#' # Example 1: Tracks without missing footprints. Based on the Paluxy River
#' # dinosaur chase sequence (Farlow et al., 2011).
#'
#' # Create a temporary file to store the TPS data. This file will be
#' # automatically deleted when the R session ends.
#' tps.file <- tempfile(fileext = ".data")
#'
#' # Use the cat function to write TPS data directly into the temporary file.
#' # This TPS data includes footprint coordinates for two different tracks,
#' # along with associated metadata.
#' cat(
#'   "LM=30",
#'   "211.00000 146.00000", "137.00000 327.00000", "246.00000 447.00000",
#'   "201.00000 609.00000", "291.00000 704.00000", "246.00000 851.00000",
#'   "361.00000 959.00000", "327.00000 1129.00000", "438.00000 1226.00000",
#'   "402.00000 1383.00000", "518.00000 1482.00000", "440.00000 1584.00000",
#'   "528.00000 1696.00000", "469.00000 1852.00000", "563.00000 1967.00000",
#'   "497.00000 2112.00000", "597.00000 2211.00000", "550.00000 2365.00000",
#'   "647.00000 2492.00000", "576.00000 2643.00000", "656.00000 2804.00000",
#'   "560.00000 2922.00000", "621.00000 3046.00000", "530.00000 3135.00000",
#'   "588.00000 3291.00000", "496.00000 3430.00000", "566.00000 3570.00000",
#'   "470.00000 3704.00000", "537.00000 3845.00000", "436.00000 4014.00000",
#'   "IMAGE=Sauropod.png", "ID=0",
#'   "LM=25",
#'   "43.00000 317.00000", "125.00000 463.00000", "102.00000 587.00000",
#'   "162.00000 745.00000", "170.00000 897.00000", "249.00000 1044.00000",
#'   "253.00000 1206.00000", "344.00000 1349.00000", "331.00000 1494.00000",
#'   "386.00000 1612.00000", "358.00000 1753.00000", "413.00000 1886.00000",
#'   "391.00000 2027.00000", "470.00000 2139.00000", "486.00000 2341.00000",
#'   "491.00000 2538.00000", "550.00000 2657.00000", "485.00000 2794.00000",
#'   "516.00000 2921.00000", "460.00000 3065.00000", "506.00000 3231.00000",
#'   "507.00000 3388.00000", "566.00000 3574.00000", "539.00000 3752.00000",
#'   "495.00000 3945.00000", "IMAGE=Theropod.png", "ID=1",
#'   file = tps.file, sep = "\n"
#' )
#'
#' # Call the tps_to_track function to convert the TPS data in the file
#' # into a track object. The 'scale' argument sets the scaling factor
#' # for the coordinates, and 'missing=FALSE' indicates that no landmarks
#' # are missing in the dataset.
#' tps_to_track(tps.file, scale = 0.004341493, missing = FALSE, NAs = NULL)
#'
#'
#' # Example 2: Tracks with missing footprints. Based on dinosaur tracks from
#'  # the Mount Tom (Ostrom, 1972).
#'
#' # Create a temporary file to store the TPS data. This file will be
#' # automatically deleted when the R session ends.
#' tps.file <- tempfile(fileext = ".data")
#'
#' # Define a vector of lines representing the TPS data. Each line in this vector
#' # represents either the landmark coordinates for a specific image or metadata
#' # about the image (such as the image filename and ID).
#' data_lines <- c(
#'   "LM=10", "9328.00000 3441.00000", "9113.00000 3587.00000", "8982.00000 3782.00000",
#'   "8794.00000 3907.00000", "8637.00000 4080.00000", "8480.00000 4230.00000",
#'   "8361.00000 4427.00000", "8157.00000 4503.00000", "8048.00000 4659.00000",
#'   "7834.00000 4804.00000", "IMAGE=Track 01.png", "ID=0", "LM=10",
#'   "8957.00000 3299.00000", "8809.00000 3459.00000", "8674.00000 3603.00000",
#'   "8490.00000 3782.00000", "8339.00000 3945.00000", "8192.00000 4170.00000",
#'   "8051.00000 4415.00000", "7925.00000 4578.00000", "7753.00000 4754.00000",
#'   "7640.00000 4976.00000", "IMAGE=Track 02.png", "ID=1", "LM=6",
#'   "8731.00000 3412.00000", "8543.00000 3638.00000", "8355.00000 3804.00000",
#'   "8129.00000 4070.00000", "7944.00000 4290.00000", "7778.00000 4459.00000",
#'   "IMAGE=Track 03.png", "ID=2", "LM=6", "8148.00000 3518.00000",
#'   "7982.00000 3741.00000", "7822.00000 3967.00000", "7715.00000 4183.00000",
#'   "7549.00000 4355.00000", "7374.00000 4547.00000", "IMAGE=Track 04.png",
#'   "ID=3", "LM=3", "7383.00000 3161.00000", "7254.00000 3365.00000",
#'   "7073.00000 3578.00000", "IMAGE=Track 05.png", "ID=4",
#'   "LM=4", "7082.00000 3192.00000", "6932.00000 3374.00000", "6743.00000 3547.00000",
#'   "6577.00000 3741.00000", "IMAGE=Track 06.png", "ID=5",
#'   "LM=7", "7113.00000 3293.00000", "7032.00000 3496.00000", "6922.00000 3801.00000",
#'   "6790.00000 3970.00000", "6665.00000 4186.00000", "6508.00000 4343.00000",
#'   "6414.00000 4525.00000", "IMAGE=Track 07.png", "ID=6",
#'   "LM=7", "6653.00000 3324.00000", "6446.00000 3525.00000", "6320.00000 3779.00000",
#'   "6160.00000 4007.00000", "6082.00000 4243.00000", "5903.00000 4418.00000",
#'   "5762.00000 4719.00000", "IMAGE=Track 08.png", "ID=7",
#'   "LM=6", "6461.00000 3973.00000", "6775.00000 3775.00000", "7069.00000 3666.00000",
#'   "7383.00000 3543.00000", "7649.00000 3446.00000", "7878.00000 3374.00000",
#'   "IMAGE=Track 09.png", "ID=8", "LM=3", "5405.00000 3785.00000",
#'   "5543.00000 3534.00000", "5671.00000 3337.00000", "IMAGE=Track 10.png",
#'   "ID=9", "LM=3", "5182.00000 3860.00000", "5305.00000 3616.00000",
#'   "5411.00000 3362.00000", "IMAGE=Track 11.png", "ID=10",
#'   "LM=3", "5342.00000 3186.00000", "5276.00000 3427.00000", "5145.00000 3694.00000",
#'   "IMAGE=Track 12.png", "ID=11", "LM=6", "3339.00000 2221.00000",
#'   "3270.00000 2155.00000", "3183.00000 2062.00000", "3113.00000 2023.00000",
#'   "3066.00000 1913.00000", "3019.00000 1832.00000", "IMAGE=Track 13.png",
#'   "ID=12", "LM=4", "4535.00000 1652.00000", "4430.00000 1911.00000",
#'   "4415.00000 2113.00000", "4347.00000 2333.00000", "IMAGE=Track 14.png",
#'   "ID=13", "LM=6", "3962.00000 1950.00000", "3857.00000 2132.00000",
#'   "3765.00000 2387.00000", "3660.00000 2637.00000", "3512.00000 2845.00000",
#'   "3392.00000 3089.00000", "IMAGE=Track 15.png", "ID=14",
#'   "LM=9", "3711.00000 1603.00000", "3638.00000 1796.00000", "3527.00000 1971.00000",
#'   "3407.00000 2164.00000", "3289.00000 2382.00000", "3146.00000 2588.00000",
#'   "2992.00000 2742.00000", "2840.00000 2958.00000", "2664.00000 3061.00000",
#'   "IMAGE=Track 16.png", "ID=15", "LM=3", "2371.00000 2826.00000",
#'   "2332.00000 2770.00000", "2279.00000 2721.00000", "IMAGE=Track 17.png",
#'   "ID=16", "LM=6", "3069.00000 1687.00000", "2934.00000 1971.00000",
#'   "2895.00000 2230.00000", "2750.00000 2475.00000", "2664.00000 2755.00000",
#'   "2578.00000 2926.00000", "IMAGE=Track 18.png", "ID=17",
#'   "LM=3", "2214.00000 1098.00000", "2067.00000 1462.00000", "1962.00000 1725.00000",
#'   "IMAGE=Track 19.png", "ID=18", "LM=4", "1630.00000 1753.00000",
#'   "1570.00000 2027.00000", "1495.00000 2282.00000", "1416.00000 2500.00000",
#'   "IMAGE=Track 20.png", "ID=19", "LM=4", "1000.00000 1410.00000",
#'   "817.00000 1536.00000", "601.00000 1736.00000", "437.00000 1920.00000",
#'   "IMAGE=Track 21.png", "ID=20", "LM=4", "1265.00000 1335.00000",
#'   "1125.00000 1446.00000", "1089.00000 1661.00000", "975.00000 1866.00000",
#'   "IMAGE=Track 22.png", "ID=21", "LM=3", "1928.00000 2874.00000",
#'   "2392.00000 2837.00000", "2726.00000 2736.00000", "IMAGE=Track 23.png",
#'   "ID=22"
#' )
#'
#' # Write the data lines to the temporary TPS file. This TPS file now contains all
#' # the data necessary to define the landmark coordinates and associated metadata
#' # for each track.
#' writeLines(data_lines, con = tps.file)
#'
#' # Define a matrix representing the footprints that are missing from the dataset.
#' # In this example, the matrix 'NAs' specifies that footprint 7 is missing in track 3.
#' NAs <- matrix(c(7, 3), nrow = 1, ncol = 2)
#'
#' # Call the tps_to_track function, which will convert the TPS data in the file
#' # to a track object. The 'scale' argument sets the scaling factor for the coordinates,
#' # 'missing' specifies whether missing footprints should be handled, 'NAs' provides
#' # the missing footprints matrix, and 'R.L.side' specifies which side (Right or Left)
#' # is the first footprint of each track.
#' tps_to_track(tps.file, scale = 0.004411765, missing = TRUE, NAs = NAs,
#'              R.L.side = c("R", "L", "L", "L", "R", "L", "R", "R", "L", "L", "L",
#'                           "L", "L", "R", "R", "L", "R", "R", "L", "R", "R",
#'                           "R", "R"))
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


tps_to_track <- function(file, scale=NULL, missing=FALSE, NAs=NULL, R.L.side=NULL) {

  ## Errors and Warnings----

  # Check if the 'scale' argument is provided
  if(is.null(scale)) {
    stop("The 'scale' argument is missing. Please provide the scale in pixels per meter.")
  }

  # Validate 'scale': must be a positive numeric value
  if(!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("The 'scale' argument must be a single positive numeric value.")
  }

  # Warn if 'NAs' is provided but 'missing' is set to FALSE
  if(!missing && !is.null(NAs)) {
    warning("The 'NAs' argument will be ignored because 'missing' is set to FALSE.")
  }

  # Validate 'NAs': must be provided if 'missing' is TRUE
  if(missing && is.null(NAs)) {
    stop("The 'NAs' argument must be provided if 'missing' is set to TRUE.")
  }

  # Validate 'R.L.side': must be provided if 'missing' is TRUE
  if(missing && is.null(R.L.side)) {
    stop("The 'R.L.side' argument must be provided if 'missing' is set to TRUE.")
  }

  # Validate 'missing': must be a single logical value
  if(!is.logical(missing) || length(missing) != 1) {
    stop("The 'missing' argument must be a single logical value: TRUE or FALSE.")
  }

  # Validate 'NAs': must be a matrix with two columns containing positive integers
  if(!is.null(NAs)) {
    if(!is.matrix(NAs) || ncol(NAs) != 2) {
      stop("The 'NAs' argument must be a matrix with two columns.")
    }
    if(any(NAs <= 0)) {
      stop("The 'NAs' matrix must contain positive integers.")
    }
  }

  # Validate 'R.L.side': must contain only 'R' or 'L' values
  if(!is.null(R.L.side)) {
    if(!all(R.L.side %in% c("R", "L"))) {
      stop("The 'R.L.side' vector must contain only 'R' or 'L' values.")
    }
  }

  ##Code----

  # Read the lines from the file
  a = readLines(file)

  # Identify lines containing landmarks (LM) and IDs
  LM = grep("LM", a)
  ID.ind = grep("ID", a)

  # Extract image names from the lines preceding the ID lines
  images = basename(gsub("(IMAGE=)(.*)", "\\2", a[ID.ind - 1]))

  # Extract the number of rows for each landmark set
  skip = LM
  nrows = as.numeric(gsub("(LM=)([0-9])", "\\2", grep("LM", a, value=T)))
  l = length(LM)

  # Initialize a list to store landmark data frames
  landmarks = vector("list", l)

  for (i in 1:l) {
    # Read the landmark data into a data frame
    landmarks[i] = list(data.frame(
      read.table(file=file, header=F, skip=LM[i],
                 nrows=nrows[i], col.names=c("X", "Y")),
      IMAGE = images[i],
      ID = read.table(file=file, header=F, skip=ID.ind[i]-1,
                      nrows=1, sep="=", col.names="ID")[2,]
    ))
  }

  # Create a data frame from the landmarks list
  data_frame <- landmarks

  if(missing == TRUE){
    ##Inferring missing footprints

    ###Including NAs
    # Levels and corresponding numbers from the NAs matrix
    levels <- levels(as.factor(NAs[,1]))
    levelsnum <- as.numeric(levels)

    # Include missing footprints in the data frame
    data_frame <- landmarks

    for (i in levelsnum) {
      data_frame[[i]] <- berryFunctions::insertRows(landmarks[[i]], c(NAs[which(NAs[,1]==i),2]), new = NA)
    }

    # Assign the correct IMAGE and ID to missing footprints
    for (i in levelsnum) {
      data_frame[[i]][c(NAs[which(NAs[,1]==i),2]),]$IMAGE <- levels(as.factor(data_frame[[i]]$IMAGE))
    }

    for (i in levelsnum) {
      data_frame[[i]][c(NAs[which(NAs[,1]==i),2]),]$ID <- levels(as.factor(data_frame[[i]]$ID))
    }

    ###Calculating track width
    vectorwidth <- c()
    meanwidth <- c()

    for (j in 1:length(data_frame)) {
      vectorwidth <- c()

      # Calculate the width of the track using the area of triangles
      for (i in 1:(length(data_frame[[j]][,1])-2)) {
        df <- data_frame[[j]][i:(i+2),1:2]

        x1 <- df[1,1]
        x2 <- df[2,1]
        x3 <- df[3,1]
        y1 <- df[1,2]
        y2 <- df[2,2]
        y3 <- df[3,2]

        Area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        Base = dist(df[c(1,3), c(1,2)], method = "euclidean")
        Height = abs((Area * 2) / Base)

        vectorwidth[i] <- Height
      }

      # Compute the mean width of the track
      meanwidth[j] <- mean(vectorwidth[which(!is.na(vectorwidth))])
    }

    ###Extrapolating missing footprints
    for (i in 1:length(data_frame)) {
      vec <- c(1:length(data_frame[[i]][,1]))
      vec2 <- c(1:length(data_frame[[i]][,1]))

      # Assign "R" or "L" based on the R.L.side vector
      if (R.L.side[i] == "L") {
        vec2[which((vec %% 2) != 0)] <- rep("L", length(which((vec %% 2) != 0)))
        vec2[which((vec %% 2) == 0)] <- rep("R", length(which((vec %% 2) == 0)))
      }

      if (R.L.side[i] == "R") {
        vec2[which((vec %% 2) != 0)] <- rep("R", length(which((vec %% 2) != 0)))
        vec2[which((vec %% 2) == 0)] <- rep("L", length(which((vec %% 2) == 0)))
      }

      data_frame[[i]]$Side <- vec2
    }

    # Extrapolate missing footprints based on calculated track width
    for (j in 1:length(data_frame)) {
      for (i in 1:length(NAs[,1])) {
        x1 <- data_frame[[NAs[i,1]]][NAs[i,2]-1,1]
        y1 <- data_frame[[NAs[i,1]]][NAs[i,2]-1,2]
        x2 <- data_frame[[NAs[i,1]]][NAs[i,2]+1,1]
        y2 <- data_frame[[NAs[i,1]]][NAs[i,2]+1,2]

        # Calculate distance based on side (R or L)
        if(data_frame[[NAs[i,1]]][NAs[i,2],5] == "R") {
          dist <- -1 * meanwidth[j]
        }
        if(data_frame[[NAs[i,1]]][NAs[i,2],5] == "L") {
          dist <- meanwidth[j]
        }

        x3 = (x1 + x2) / 2
        y3 = (y1 + y2) / 2

        b = x2 - x1
        a = y1 - y2

        norm = sqrt(a * a + b * b)
        a = a / norm
        b = b / norm

        x4 = x3 + a * dist
        y4 = y3 + b * dist

        data_frame[[NAs[i,1]]][NAs[i,2],1] <- x4
        data_frame[[NAs[i,1]]][NAs[i,2],2] <- y4
      }
    }
  }

  ###Tracing medial tracks

  # Create medial tracks by averaging consecutive landmarks
  landmarks2 <- data_frame

  landmarks3 <- list()
  for (i in 1:length(landmarks2)) {
    landmarks3[[i]] <- as.data.frame(matrix(ncol=ncol(landmarks2[[i]]), nrow=(nrow(landmarks2[[i]])-1)))
    colnames(landmarks3[[i]]) <- colnames(landmarks2[[i]])

    landmarks3[[i]][,3] <- rep(landmarks2[[i]][1,3], length(landmarks3[[i]][,3]))
    landmarks3[[i]][,4] <- rep(landmarks2[[i]][1,4], length(landmarks3[[i]][,4]))
    landmarks3[[i]] <- landmarks3[[i]][,1:4]

    # Compute the average position for each consecutive pair of landmarks
    for (j in 1:length(landmarks3[[i]][,1])) {
      landmarks3[[i]][j,1] <- (landmarks2[[i]][j,1] + landmarks2[[i]][j+1,1]) / 2
    }

    for (j in 1:length(landmarks3[[i]][,2])) {
      landmarks3[[i]][j,2] <- (landmarks2[[i]][j,2] + landmarks2[[i]][j+1,2]) / 2
    }
  }

  # Scale the coordinates and create trajectories from the coordinates
  for (i in 1:length(landmarks3)) {
    landmarks3[[i]][,1:2] <- landmarks3[[i]][,1:2] * scale
    landmarks3[[i]] <- TrajFromCoords(landmarks3[[i]])
  }
  names(landmarks3) <- paste0("Track_", str_pad(1:length(LM), nchar(length(LM)), pad = "0"), sep = "")

  # Mark missing footprints
  for (i in 1:length(data_frame)) {
    data_frame[[i]]$missing <- rep("Actual", length(data_frame[[i]][,1]))
  }

  if(missing == TRUE) {
    levels <- levels(as.factor(NAs[,1]))
    levelsnum <- as.numeric(levels)
    for (i in levelsnum) {
      data_frame[[i]][c(NAs[which(NAs[,1]==i),2]),]$missing <- rep("Inferred", length(which(NAs[,1]==i)))
    }
  }

  # Scale the X and Y coordinates
  for (i in 1:length(data_frame)) {
    data_frame[[i]]$X <- data_frame[[i]]$X * scale
    data_frame[[i]]$Y <- data_frame[[i]]$Y * scale
  }

  # Create a list to return the results
  listdata <- list()
  listdata[[1]] <- landmarks3
  listdata[[2]] <- data_frame

  print(listdata)
}
