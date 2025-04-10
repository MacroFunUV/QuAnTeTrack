% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tps_to_track.R
\name{tps_to_track}
\alias{tps_to_track}
\title{Transform a *.tps file into a \code{track} R object}
\usage{
tps_to_track(file, scale = NULL, missing = FALSE, NAs = NULL, R.L.side = NULL)
}
\arguments{
\item{file}{A *.tps file containing (x,y) coordinates of footprints in tracks.}

\item{scale}{A numeric value specifying the scale in meters per pixel.}

\item{missing}{A logical value indicating whether there are missing footprints in any track to be interpolated: \code{TRUE}, or \code{FALSE} (the default).}

\item{NAs}{A matrix with two columns indicating which missing footprints will be interpolated.
The first column gives the number of the track containing missing footprints, and the second column gives the number of the footprint that is missing within this track.
The number of rows is equal to the total number of missing footprints in the sample.}

\item{R.L.side}{A character vector specifying the side of the first footprint of each track.
Only needed if \code{missing} is set to \code{TRUE}.
The length of the vector must be equal to the total number of tracks in the sample.
\itemize{
\item \code{"L"}: first footprint corresponds to the left foot.
\item \code{"R"}: first footprint corresponds to the right foot.
}}
}
\value{
A \code{track} R object, which is a list consisting of two elements:
\itemize{
\item \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
\item \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
}
}
\description{
\code{tps_to_track()} reads a *.tps file containing footprint coordinates of one or several tracks and transforms it into a \code{track} R object.
}
\details{
It is highly recommended that the *.tps file is built using the TPS software (Rohlf 2008, 2009).
Tracks with a different number of footprints (i.e., landmarks) are allowed.
This function transforms the coordinates of the footprints of each track into a set of trajectory coordinates.
Each point of the trajectory is calculated as: \deqn{Point_i(x,y)= (Footprint_i(x,y) + Footprint_{i+1}(x,y)/2}

The number of points of the resulting trajectory is \eqn{n_{footprints} - 1}.

If \code{missing} is set to \code{TRUE}, missing footprints can be interpolated.
This interpolation is based on adjacent footprints and the provided side information.
}
\section{Logo}{

\if{html}{\figure{Logo.png}{options: width=30\%}}
}

\examples{
# Example 1: Tracks without missing footprints. Based on the Paluxy River
# dinosaur chase sequence (Farlow et al., 2011).

# Load the example TPS file provided in the QuAnTeTrack package.
# This TPS data includes footprint coordinates for different tracks,
# along with associated metadata.
tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")

# Call the tps_to_track function to convert the TPS data in the file
# into a track object. The 'scale' argument sets the scaling factor
# for the coordinates, and 'missing=FALSE' indicates that no landmarks
# are missing in the dataset.
tps_to_track(tpsPaluxyRiver, scale = 0.004341493, missing = FALSE, NAs = NULL)


# Example 2: Tracks with missing footprints. Based on dinosaur tracks from
# the Mount Tom (Ostrom, 1972).

# Load the example TPS file provided in the QuAnTeTrack package.
# This TPS data includes footprint coordinates for different tracks,
# along with associated metadata.
tpsMountTom <- system.file("extdata", "MountTom.tps", package = "QuAnTeTrack")

# Define a matrix representing the footprints that are missing from the dataset.
# In this example, the matrix 'NAs' specifies that footprint 7 is missing in track 3.
NAs <- matrix(c(7, 3), nrow = 1, ncol = 2)

# Call the tps_to_track function, which will convert the TPS data in the file
# to a track object. The 'scale' argument sets the scaling factor for the coordinates,
# 'missing' specifies whether missing footprints should be handled, 'NAs' provides
# the missing footprints matrix, and 'R.L.side' specifies which side (Right or Left)
# is the first footprint of each track.
tps_to_track(tpsMountTom,
  scale = 0.004411765, missing = TRUE, NAs = NAs,
  R.L.side = c(
    "R", "L", "L", "L", "R", "L", "R", "R", "L", "L", "L",
    "L", "L", "R", "R", "L", "R", "R", "L", "R", "R",
    "R", "R"
  )
)

}
\references{
Farlow, J. O., O’Brien, M., Kuban, G. J., Dattilo, B. F., Bates, K. T., Falkingham, P. L., & Piñuela, L. (2012). Dinosaur Tracksites of the Paluxy River Valley (Glen Rose Formation, Lower Cretaceous), Dinosaur Valley State Park, Somervell County, Texas. In Proceedings of the V International Symposium about Dinosaur Palaeontology and their Environment (pp. 41-69). Burgos: Salas de los Infantes.

Ostrom, J. H. (1972). Were some dinosaurs gregarious?. Palaeogeography, Palaeoclimatology, Palaeoecology, 11(4), 287-301.

Rohlf, F. J. 2008. TPSUTIL. Version 1.40. Department of Ecology and Evolution, State University of New York.
Available at \url{https://sbmorphometrics.org/}.

Rohlf, F. J. 2009. tpsDig. Version 2.14. Department of Ecology and Evolution, State University of New York.
Available at \url{https://sbmorphometrics.org/}.
}
\author{
Humberto G. Ferrón

humberto.ferron@uv.es

Macroevolution and Functional Morphology Research Group (www.macrofun.es)

Cavanilles Institute of Biodiversity and Evolutionary Biology

Calle Catedrático José Beltrán Martínez, nº 2

46980 Paterna - Valencia - Spain

Phone: +34 (9635) 44477
}
