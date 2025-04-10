% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_direction.R
\name{plot_direction}
\alias{plot_direction}
\title{Plot direction data in tracks.}
\usage{
plot_direction(
  data,
  plot_type = "boxplot",
  angle_range = 30,
  y_labels_position = 90,
  y_breaks_manual = NULL
)
}
\arguments{
\item{data}{A \code{track} R object, which is a list consisting of two elements:
\itemize{
\item \strong{\code{Trajectories}}: A list of interpolated trajectories, where each trajectory is a series of midpoints between consecutive footprints.
\item \strong{\code{Footprints}}: A list of data frames containing footprint coordinates, metadata (e.g., image reference, ID), and a marker indicating whether the footprint is actual or inferred.
}}

\item{plot_type}{A character string indicating the type of plot to generate.
The options are \code{"boxplot"}, \code{"polar_steps"}, \code{"polar_average"}, and \code{"faceted"}.
Default is \code{"boxplot"}.}

\item{angle_range}{A numeric value specifying the width of the bins (in degrees)
used for polar plots. Default is \code{30} degrees.}

\item{y_labels_position}{A numeric value specifying the position (in degrees)
of the y-axis labels in the polar plots. Default is \code{90} degrees.}

\item{y_breaks_manual}{A numeric vector specifying manual breaks for the y-axis
in polar plots. If \code{NULL}, the breaks are calculated automatically. Default is \code{NULL}.}
}
\value{
A \code{ggplot} object that displays the specified plot type. The \pkg{ggplot2} package is used for plotting.
}
\description{
\code{plot_direction()}  generates different types of plots to visualize the direction data from \code{track} R objects. The function allows for the creation of boxplots, polar histograms of step directions, polar histograms of average directions per track, and faceted polar histograms.
}
\details{
The \code{plot_direction()} function provides four types of plots:

\itemize{
\item \code{"boxplot"}: A boxplot showing the distribution of step direction values
for each track.
\item \code{"polar_steps"}: A polar plot showing the frequency of step
in different direction bins.
\item \code{"polar_average"}: A polar plot showing the frequency of average
directions per track in different direction bins.
\item \code{"faceted"}: A polar plot similar to \code{"polar_steps"} but
faceted by track.
}

The \code{angle_range} parameter defines the bin width for the polar plots,
and \code{y_labels_position} allows for adjusting the position of y-axis labels.
The \code{y_breaks_manual} parameter lets users manually specify the breaks
on the y-axis for finer control over the appearance of the polar plots.
}
\section{Logo}{

\if{html}{\figure{Logo.png}{options: width=30\%}}
}

\examples{
# Example 1: Boxplot of Direction Data in MountTom Dataset
plot_direction(MountTom, plot_type = "boxplot")

# Example 2: Polar Plot of Step Directions in MountTom Dataset
plot_direction(MountTom, plot_type = "polar_steps")

# Example 3: Polar Plot of Average Directions Per Track in MountTom Dataset
plot_direction(MountTom, plot_type = "polar_average")

# Example 4: Faceted Polar Plot of Step Directions in MountTom Dataset
plot_direction(MountTom, plot_type = "faceted")

# Example 5: Polar Plot with Custom Angle Range in MountTom Dataset
plot_direction(MountTom, plot_type = "polar_steps", angle_range = 90)

# Example 6: Polar Plot with Custom Y-Axis Labels and Breaks in MountTom Dataset
plot_direction(MountTom,
  plot_type = "polar_steps", y_labels_position = 0,
  y_breaks_manual = c(0, 15, 30, 45, 60)
)

# Example 7: Boxplot of Direction Data in PaluxyRiver Dataset
plot_direction(PaluxyRiver, plot_type = "boxplot")

# Example 8: Polar Plot of Step Directions in PaluxyRiver Dataset
plot_direction(PaluxyRiver, plot_type = "polar_steps")

# Example 9: Polar Plot of Average Directions Per Track with Custom Breaks in PaluxyRiver Dataset
plot_direction(PaluxyRiver,
  plot_type = "polar_average",
  y_breaks_manual = c(1, 2)
)

# Example 10: Faceted Polar Plot of Step Directions in PaluxyRiver Dataset
plot_direction(PaluxyRiver, plot_type = "faceted")

# Example 11: Polar Plot of Average Directions Per Track with Custom Breaks in PaluxyRiver Dataset
plot_direction(PaluxyRiver,
  plot_type = "polar_average",
  y_breaks_manual = c(1, 2)
)

# Example 12: Polar Plot with Custom Y-Axis Labels in PaluxyRiver Dataset
plot_direction(PaluxyRiver,
  plot_type = "polar_steps",
  y_labels_position = -90
)

}
\seealso{
\code{\link{tps_to_track}}, \code{\link{test_direction}}
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
