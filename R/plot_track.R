#' Plot trajectories and footprints
#'
#' \code{plot_track()} visualizes trajectory and footprint data from code{track} R objects in various ways, allowing for the plotting of trajectories, footprints, or both combined, with customizable aesthetics.
#'
#' @param data A \code{trackway} R object, which is a list consisting of two elements:
#'
#'    * \strong{\code{Trajectories}}: A list of trajectories representing trackway midlines,
#'      interpolated by connecting the midpoints of successive left–right footprint pairs
#'      (i.e., footprints linked by pace lines). Includes columns \code{X}, \code{Y},
#'      \code{IMAGE}, \code{ID}, and \code{Side} (set to \code{"Medial"}).
#'
#'    * \strong{\code{Footprints}}: A list of data frames containing footprint coordinates
#'      and associated metadata, with a \code{Side} column (\code{"R"} or \code{"L"})
#'      and a \code{missing} marker (\code{"Actual"} or \code{"Inferred"}).
#' @param plot Type of plot to generate. Options are \code{"FootprintsTrajectories"} (default), \code{"Trajectories"}, or \code{"Footprints"}. Determines what elements are included in the plot.
#' @param colours A vector of colors to be used for different trajectories If \code{NULL}, defaults to black. The length of this vector should match the number of trackways in the data.
#' @param cex.f The size of the points representing footprints. Default is \code{2.5}. If \code{seq.foot = TRUE}, this controls the text size of footprint sequence numbers.
#' @param shape.f A vector of shapes to be used for representing footprints in different trackways. If \code{NULL}, defaults to \code{19} (solid circle). The length of this vector should match the number of trackways in the data.
#' @param alpha.f The transparency of the points representing footprints. Default is \code{0.5}.
#' @param cex.t The thickness of the trajectory lines. Default is \code{0.5}.
#' @param alpha.t The transparency of the trajectory lines. Default is \code{1}.
#' @param plot.labels Logical indicating whether to add labels to each trackway. Default is \code{FALSE}.
#' @param labels A vector of labels for each trackway. If \code{NULL}, labels are automatically generated from trackway names.
#' @param box.p Padding around label boxes, used only if \code{plot.labels} is \code{TRUE}. Adjusts the spacing around the label text.
#' @param cex.l The size of the labels. Default is \code{3.88}.
#' @param alpha.l The transparency of the labels. Default is \code{0.5}.
#' @param arrow.t Logical indicating whether to add an arrowhead at the end of each trajectory to verify direction. Default is \code{FALSE}.
#' @param arrow.size Numeric controlling the arrowhead size (in inches, passed to \code{grid::unit()}). Default is \code{0.15}.
#' @param seq.foot Logical indicating whether to display the footprint sequence numbers instead of footprint points. Default is \code{FALSE}. This is useful to verify that footprints are in the correct order along each trackway. When \code{TRUE}, the size of the numbers is controlled by \code{cex.f}.
#'
#' @details
#' The \code{plot_track()} function is designed as a diagnostic and exploratory tool.
#' Its primary purpose is to display the raw spatial data (footprint coordinates and
#' interpolated trajectories) that have been digitized, so that users can
#' visually confirm data integrity before conducting quantitative analyses. This
#' includes checking whether footprints are in the correct order, whether trackways are
#' oriented consistently, and whether interpolated trajectories align with the raw
#' footprint data.
#'
#' Importantly, these plots are not intended to replace traditional ichnological
#' illustrations. Hand-drawn maps and outlines often convey information that is not
#' captured in raw coordinate plots, such as tridactyl morphology, manus/pes
#' distinction, taxonomic attribution, or trackway orientation, and they frequently
#' provide clearer and more communicative visual summaries of ichnological material.
#'
#' By contrast, \code{plot_track()} focuses on plotting digitized data as they are,
#' without additional interpretation, stylization, or symbolic annotation. The goal is
#' to offer a reproducible, data-driven representation that complements, rather than
#' supplants, traditional methods. Users are encouraged to treat these plots as
#' quality-control visualizations that help detect potential errors or inconsistencies
#' in the raw data prior to downstream analyses, while continuing to rely on classical
#' ichnological illustrations for detailed morphological and taxonomic interpretation.
#'
#' @return A \code{ggplot} object that displays the specified plot type, including trajectories, footprints, or both, from \code{trackway} R objects. The \pkg{ggplot2} package is used for plotting.
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
#' # Example 1: Basic Plot with Default Settings - MountTom Dataset
#' plot_track(MountTom)
#'
#' # Example 2: Basic Plot with Default Settings - PaluxyRiver Dataset
#' plot_track(PaluxyRiver)
#'
#' # Example 3: Plot Trajectories Only - MountTom Dataset
#' plot_track(MountTom, plot = "Trajectories")
#'
#' # Example 4: Plot Footprints Only - PaluxyRiver Dataset
#' plot_track(PaluxyRiver, plot = "Footprints")
#'
#' # Example 5: Custom Colors for Trajectories - MountTom Dataset
#' custom_colors <- c(
#'   "#008000", "#0000FF", "#FF0000", "#800080", "#FFA500", "#FFC0CB", "#FFFF00",
#'   "#00FFFF", "#A52A2A", "#FF00FF", "#808080", "#000000", "#006400", "#00008B",
#'   "#8B0000", "#FF8C00", "#008B8B", "#A9A9A9", "#000080", "#808000", "#800000",
#'   "#008080", "#FFD700"
#' )
#' plot_track(MountTom, colours = custom_colors)
#'
#' # Example 6: Larger Footprints and wider Trajectories Lines - PaluxyRiver Dataset
#' plot_track(PaluxyRiver, cex.f = 5, cex.t = 2)
#'
#' # Example 7: Semi-Transparent Footprints and Trajectories - MountTom Dataset
#' plot_track(MountTom, alpha.f = 0.5, alpha.t = 0.5)
#'
#' # Example 8: Different Shapes for Footprints - PaluxyRiver Dataset
#' plot_track(PaluxyRiver, shape.f = c(16, 17))
#'
#' # Example 9: Plot with Labels for Trajectories - MountTom Dataset
#' labels <- paste("Trackway", seq_along(MountTom[[1]]))
#' plot_track(MountTom, plot.labels = TRUE, labels = labels, cex.l = 4, box.p = 0.3, alpha.l = 0.7)
#'
#' # Example 10: Custom Colors and Shapes for Footprints Only - PaluxyRiver Dataset
#' plot_track(PaluxyRiver, plot = "Footprints", colours = c("purple", "orange"), shape.f = c(15, 18))
#'
#' # Example 11: Wider Line Size & Custom Colors for Trajectories Only - MountTom Dataset
#' plot_track(MountTom, plot = "Trajectories", cex.t = 1.5, colours = custom_colors)
#'
#' # Example 12: Black Footprints and Trajectories with Labels - PaluxyRiver Dataset
#' plot_track(PaluxyRiver,
#'   colours = NULL, shape.f = c(16, 16), plot.labels = TRUE,
#'   labels = c("Saurpod", "Theropod"), cex.l = 2, alpha.l = 0.5
#' )
#'
#' # Example 13: Add arrowheads to show trajectory direction
#' plot_track(MountTom, plot = "Trajectories", arrow.t = TRUE, arrow.size = 0.2)
#'
#' # Example 14: Show footprint sequence numbers (quality control)
#' plot_track(MountTom, plot = "Footprints", seq.foot = TRUE)
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 alpha
#' @importFrom ggplot2 coord_fixed
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_light
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr bind_rows
#' @importFrom grid arrow unit
#'
#' @seealso \code{\link{tps_to_track}}
#'
#' @export
plot_track <- function(
    data,
    plot = "FootprintsTrajectories",
    colours = NULL,
    cex.f = NULL,
    shape.f = NULL,
    alpha.f = NULL,
    cex.t = NULL,
    alpha.t = NULL,
    plot.labels = NULL,
    labels = NULL,
    box.p = NULL,
    cex.l = NULL,
    alpha.l = NULL,
    arrow.t = FALSE,
    arrow.size = 0.15,
    seq.foot = FALSE
) {

  ## Set default values if arguments are NULL----
  if (is.null(cex.f)) cex.f <- 2.5
  if (is.null(cex.t)) cex.t <- 0.5
  if (is.null(cex.l)) cex.l <- 3.88

  if (is.null(alpha.f)) alpha.f <- 0.5
  if (is.null(alpha.t)) alpha.t <- 1
  if (is.null(alpha.l)) alpha.l <- 0.5

  if (is.null(colours)) colours <- c(rep("#000000", length(data[[1]])))

  if (is.null(shape.f)) shape.f <- c(rep(19, length(data[[1]])))

  if (is.null(plot.labels)) plot.labels <- FALSE

  if (is.null(box.p)) box.p <- 0.25

  ## Errors and Warnings----
  if (!is.list(data) || length(data) < 2) {
    stop("The 'data' argument must be a 'trackway' R object, which is a list consisting of two elements.")
  }
  if (!is.list(data[[1]]) || !is.list(data[[2]])) {
    stop("The two elements of 'data' must be lists.")
  }
  if (!is.null(colours) && length(colours) != length(data[[1]])) {
    stop("Error: The length of 'colours' must match the number of trackways in the data.")
  }
  if (!plot %in% c("FootprintsTrajectories", "Trajectories", "Footprints")) {
    stop("Error: The 'plot' argument must be one of 'FootprintsTrajectories', 'Trajectories', or 'Footprints'.")
  }
  if (!is.null(labels) && !plot.labels) {
    warning("Warning: 'labels' are provided but 'plot.labels' is set to FALSE. Labels will not be displayed.")
  }
  if (!is.null(shape.f) && length(shape.f) != length(data[[1]])) {
    stop("Error: The length of 'shape.f' must match the number of trackways in the data.")
  }
  if (plot.labels && is.null(box.p)) {
    warning("Warning: 'box.p' padding around label boxes is NULL. Default padding will be used.")
  }
  if (!is.null(labels) && length(labels) != length(data[[1]])) {
    stop("Error: The length of 'labels' must match the number of trackways in the data.")
  }
  if (!is.null(cex.f) && (!is.numeric(cex.f) || cex.f <= 0)) {
    stop("Error: 'cex.f' size of footprint points must be a positive numeric value.")
  }
  if (!is.null(cex.t) && (!is.numeric(cex.t) || cex.t <= 0)) {
    stop("Error: 'cex.t' width of trajectory lines must be a positive numeric value.")
  }
  if (!is.null(cex.l) && (!is.numeric(cex.l) || cex.l <= 0)) {
    stop("Error: 'cex.l' size of labels must be a positive numeric value.")
  }
  if (!is.null(alpha.f) && (!is.numeric(alpha.f) || alpha.f < 0 || alpha.f > 1)) {
    stop("Error: 'alpha.f' transparency of footprint points must be a numeric value between 0 and 1.")
  }
  if (!is.null(alpha.t) && (!is.numeric(alpha.t) || alpha.t < 0 || alpha.t > 1)) {
    stop("Error: 'alpha.t' transparency of trajectory lines must be a numeric value between 0 and 1.")
  }
  if (!is.null(alpha.l) && (!is.numeric(alpha.l) || alpha.l < 0 || alpha.l > 1)) {
    stop("Error: 'alpha.l' transparency of labels must be a numeric value between 0 and 1.")
  }
  if (plot.labels && is.null(labels)) {
    warning("Warning: 'plot.labels' is TRUE but 'labels' is NULL. Automatically generated labels will be used.")
  }
  if (!is.null(box.p) && (!is.numeric(box.p) || box.p <= 0)) {
    stop("Error: 'box.p' padding around label boxes must be a positive numeric value.")
  }
  if (!is.logical(arrow.t) || length(arrow.t) != 1) {
    stop("Error: 'arrow.t' must be a single logical value: TRUE or FALSE.")
  }
  if (!is.numeric(arrow.size) || length(arrow.size) != 1 || is.na(arrow.size) || arrow.size <= 0) {
    stop("Error: 'arrow.size' must be a single positive numeric value.")
  }
  if (!is.logical(seq.foot) || length(seq.foot) != 1 || is.na(seq.foot)) {
    stop("Error: 'seq.foot' must be a single logical value: TRUE or FALSE.")
  }

  ## Code----
  footprints <- bind_rows(data[[2]])
  tracks <- bind_rows(data[[1]])

  # Sequence number within each trackway (IMAGE)
  footprints$SeqFoot <- ave(footprints$X, footprints$IMAGE, FUN = seq_along)

  # Labels data frame
  mat <- data.frame(matrix(nrow = length(data[[1]]), ncol = 3))
  colnames(mat) <- c("X", "Y", "Name")
  for (i in 1:length(data[[1]])) {
    mat[i, 1] <- data[[1]][[i]]$x[[1]]
    mat[i, 2] <- data[[1]][[i]]$y[[1]]
  }
  if (is.null(labels) == FALSE) {
    mat$Name <- labels
  }
  if (is.null(labels) == TRUE) {
    for (i in 1:length(data[[1]])) {
      mat$Name[i] <- gsub("_", " ", names(data[[1]][i]))
    }
  }

  # Arrow specification
  arrow_spec <- if (arrow.t) grid::arrow(type = "closed", length = grid::unit(arrow.size, "inches")) else NULL

  # Footprint layer (points OR sequence numbers) + optional shape scale
  fp_layer <- if (!seq.foot) {
    geom_point(cex = cex.f, data = footprints,
               aes(x = X, y = Y, colour = IMAGE, shape = IMAGE),
               alpha = alpha.f)
  } else {
    geom_text(data = footprints,
              aes(x = X, y = Y, label = SeqFoot, colour = IMAGE),
              alpha = alpha.f, size = cex.f)
  }
  shape_scale <- if (!seq.foot) scale_shape_manual(values = shape.f) else NULL

  # Plot Footprints and Trajectories together----
  if (plot == "FootprintsTrajectories") {
    plotfig <- ggplot() +
      geom_path(data = tracks, linewidth = cex.t,
                aes(x = x, y = y, group = IMAGE, colour = IMAGE),
                alpha = alpha.t, arrow = arrow_spec) +
      fp_layer +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
      theme_light() +
      labs(y = "m", x = "m") +
      theme(legend.position = "none") +
      scale_colour_manual(values = colours) +
      shape_scale
  }

  # Plot only Trajectories----
  if (plot == "Trajectories") {
    plotfig <- ggplot() +
      geom_path(data = tracks, linewidth = cex.t,
                aes(x = x, y = y, group = IMAGE, colour = IMAGE),
                alpha = alpha.t, arrow = arrow_spec) +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
      theme_light() +
      labs(y = "m", x = "m") +
      theme(legend.position = "none") +
      scale_colour_manual(values = colours)
  }

  # Plot only Footprints----
  if (plot == "Footprints") {
    plotfig <- ggplot() +
      fp_layer +
      coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
      theme_light() +
      labs(y = "m", x = "m") +
      theme(legend.position = "none") +
      scale_colour_manual(values = colours) +
      shape_scale
  }

  # Add labels if requested----
  if (plot.labels == TRUE) {
    options(ggrepel.max.overlaps = Inf)
    plotfig <- plotfig + geom_label_repel(
      data = mat, aes(x = X, y = Y, label = Name),
      box.padding = box.p, point.padding = 0,
      segment.color = "grey50",
      fill = alpha(c("white"), alpha.l),
      cex = cex.l
    )
  }

  return(plotfig)
}
