% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mode_velocity.R
\name{mode_velocity}
\alias{mode_velocity}
\title{Test for steady, acceleration, or deceleration along trajectories}
\usage{
mode_velocity(trackvel)
}
\arguments{
\item{trackvel}{A \code{track velocity} object  where each element corresponds to a track and contains a vector of velocity or relative stride length data.}
}
\value{
A list where each element corresponds to a trajectory from the input \code{trackvel} and contains:
\itemize{
\item \strong{correlation:} The result of the Spearman correlation test, including the correlation coefficient and \emph{p}-value.
\item \strong{trend:} A classification of the trend as "Acceleration", "Deceleration", or "Steady" based on the \emph{p}-value and the correlation coefficient.
\item If a trajectory has fewer than 3 steps, the entry contains the message "Less than three steps."
}
}
\description{
\code{mode_velocity()} evaluates the trend in velocity along each trajectory by applying Spearman's rank correlation test. The function classifies the trend into "acceleration", "deceleration", or "steady" based on the correlation and the \emph{p}-value.
}
\details{
The \code{mode_velocity()} function performs the following operations:
\itemize{
\item \strong{Spearman's Rank Correlation Test:}
\itemize{
\item This non-parametric test assesses the strength and direction of a monotonic relationship between two variables. It does not require assumptions about the normality of data or a linear relationship between velocity and step number.
\item It uses ranks rather than raw values, making it robust to outliers and suitable for detecting general trends (acceleration or deceleration) in velocity data.
}
\item \strong{Function Operation:}
\itemize{
\item For each trajectory in the \code{trackvel} list, the function calculates the Spearman correlation coefficient and the associated \emph{p}-value between velocity and step number.
\item Based on the \emph{p}-value and correlation coefficient, it classifies the trend as "acceleration", "deceleration", or "steady".
\item If a trajectory contains fewer than 3 steps, the function returns a message indicating insufficient data for correlation analysis.
}
\item \strong{Advantages:}
\itemize{
\item The non-parametric nature allows flexibility with data distributions and reduced sensitivity to outliers compared to parametric tests.
\item Effective for detecting monotonic trends (either increasing or decreasing) when the correlation is statistically significant.
}
\item \strong{Limitations:}
\itemize{
\item May be unreliable with very small sample sizes (e.g., fewer than 3 steps), providing potentially non-informative results.
\item Does not capture the magnitude of change or provide detailed insights into the rate of acceleration or deceleration.
\item Identifies monotonic trends based on statistical significance but does not distinguish between different types of monotonic relationships (e.g., steady acceleration vs. abrupt changes).
}
}

\strong{Interpretation of Results:}
\itemize{
\item \strong{Acceleration:} If the \emph{p}-value is less than 0.05 and the Spearman correlation coefficient is positive.
\item \strong{Deceleration:} If the \emph{p}-value is less than 0.05 and the Spearman correlation coefficient is negative.
\item \strong{Steady:} If the \emph{p}-value is greater than or equal to 0.05, indicating no significant monotonic relationship.
}

\strong{Usage Considerations:}
\itemize{
\item Ensure that each trajectory in \code{trackvel} has a sufficient number of steps for meaningful analysis.
\item For more detailed analysis of velocity trends, consider complementary methods such as linear or non-linear regression, or specialized change point detection techniques.
}
}
\section{Logo}{

\if{html}{\figure{Logo.png}{options: width=30\%}}
}

\examples{
# Example 1: Test for Steady, Acceleration, or Deceleration in MountTom dataset.

# Hip heights for each track in the MountTom dataset
H_mounttom <- c(
  1.380, 1.404, 1.320, 1.736, 1.364, 1.432, 1.508, 1.768, 1.600,
  1.848, 1.532, 1.532, 0.760, 1.532, 1.688, 1.620, 0.636, 1.784,
  1.676, 1.872, 1.648, 1.760, 1.612
)

# Calculate velocities using the default Method "A"
V_mounttom <- velocity_track(MountTom, H = H_mounttom)

# Test for Steady, Acceleration, or Deceleration
mode_velocity(V_mounttom)

# Example 2: Test for Steady, Acceleration, or Deceleration in PaluxyRiver dataset.

# Hip heights for each track in the PaluxyRiver dataset
H_paluxyriver <- c(3.472, 2.200)

# Specify different methods for different tracks
Method_paluxyriver <- c("A", "B")

# Calculate velocities using specified methods
V_paluxyriver <- velocity_track(PaluxyRiver,
  H = H_paluxyriver,
  method = Method_paluxyriver
)

# Test for Steady, Acceleration, or Deceleration
mode_velocity(V_paluxyriver)

}
\seealso{
\code{\link{tps_to_track}}, \code{\link{velocity_track}}, \code{\link{plot_velocity}}
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
