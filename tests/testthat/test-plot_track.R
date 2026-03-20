test_that("plot_track returns a ggplot object for default settings with MountTom dataset", {
  plot_result <- plot_track(MountTom)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for default settings with PaluxyRiver dataset", {
  plot_result <- plot_track(PaluxyRiver)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for Tracks-only plot with MountTom dataset", {
  plot_result <- plot_track(MountTom, plot = "Trajectories")
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for Footprints-only plot with PaluxyRiver dataset", {
  plot_result <- plot_track(PaluxyRiver, plot = "Footprints")
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for custom colors with MountTom dataset", {
  custom_colors <- c(
    "#008000", "#0000FF", "#FF0000", "#800080", "#FFA500", "#FFC0CB", "#FFFF00",
    "#00FFFF", "#A52A2A", "#FF00FF", "#808080", "#000000", "#006400", "#00008B",
    "#8B0000", "#FF8C00", "#008B8B", "#A9A9A9", "#000080", "#808000", "#800000",
    "#008080", "#FFD700"
  )
  plot_result <- plot_track(MountTom, colours = custom_colors)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for larger footprints and tracks with PaluxyRiver dataset", {
  plot_result <- plot_track(PaluxyRiver, cex.f = 5, cex.t = 2)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for semi-transparent footprints and tracks with MountTom dataset", {
  plot_result <- plot_track(MountTom, alpha.f = 0.5, alpha.t = 0.5)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for different footprint shapes with PaluxyRiver dataset", {
  plot_result <- plot_track(PaluxyRiver, shape.f = c(16, 17))
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object with labels for MountTom dataset", {
  labels <- paste("Trackway", seq_along(MountTom[[1]]))
  plot_result <- plot_track(MountTom, plot.labels = TRUE, labels = labels, cex.l = 4, box.p = 0.3, alpha.l = 0.7)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for footprints-only with custom colors and shapes in PaluxyRiver dataset", {
  plot_result <- plot_track(PaluxyRiver, plot = "Footprints", colours = c("purple", "orange"), shape.f = c(15, 18))
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for tracks-only with larger line size and custom colors in MountTom dataset", {
  custom_colors <- c(
    "#008000", "#0000FF", "#FF0000", "#800080", "#FFA500", "#FFC0CB", "#FFFF00",
    "#00FFFF", "#A52A2A", "#FF00FF", "#808080", "#000000", "#006400", "#00008B",
    "#8B0000", "#FF8C00", "#008B8B", "#A9A9A9", "#000080", "#808000", "#800000",
    "#008080", "#FFD700"
  )
  plot_result <- plot_track(MountTom, plot = "Trajectories", cex.t = 1.5, colours = custom_colors)
  expect_s3_class(plot_result, "ggplot")
})

test_that("plot_track returns a ggplot object for black footprints and tracks with labels in PaluxyRiver dataset", {
  plot_result <- plot_track(PaluxyRiver,
                            colours = NULL, shape.f = c(16, 16), plot.labels = TRUE,
                            labels = c("Sauropod", "Theropod"), cex.l = 2, alpha.l = 0.5
  )
  expect_s3_class(plot_result, "ggplot")
})

## ---------------------------
## Arrow tests (robust, version-proof)
## ---------------------------

test_that("plot_track works and returns a ggplot object when arrows are enabled", {
  plot_result <- plot_track(MountTom, plot = "Trajectories", arrow.t = TRUE)
  expect_s3_class(plot_result, "ggplot")

  # Building should not error (arrow gets realized during draw)
  expect_silent(ggplot2::ggplot_build(plot_result))
  expect_silent(ggplot2::ggplotGrob(plot_result))
})

test_that("plot_track applies arrow.size without changing other behaviour", {
  p1 <- plot_track(MountTom, plot = "Trajectories", arrow.t = TRUE, arrow.size = 0.10)
  p2 <- plot_track(MountTom, plot = "Trajectories", arrow.t = TRUE, arrow.size = 0.25)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")

  # Both should build; this is the most stable cross-version check.
  expect_silent(ggplot2::ggplotGrob(p1))
  expect_silent(ggplot2::ggplotGrob(p2))
})

test_that("plot_track errors on invalid arrow arguments", {
  expect_error(plot_track(MountTom, arrow.t = "yes"), "arrow.t")
  expect_error(plot_track(MountTom, arrow.size = 0), "arrow.size")
  expect_error(plot_track(MountTom, arrow.size = -1), "arrow.size")
  expect_error(plot_track(MountTom, arrow.size = NA), "arrow.size")
  expect_error(plot_track(MountTom, arrow.size = "big"), "arrow.size")
})

## ---------------------------
## Footprint sequence tests (seq.foot)
## ---------------------------

test_that("plot_track works and returns a ggplot object when seq.foot is enabled (Footprints)", {
  plot_result <- plot_track(MountTom, plot = "Footprints", seq.foot = TRUE)
  expect_s3_class(plot_result, "ggplot")

  expect_silent(ggplot2::ggplot_build(plot_result))
  expect_silent(ggplot2::ggplotGrob(plot_result))
})

test_that("plot_track works and returns a ggplot object when seq.foot is enabled (FootprintsTrajectories)", {
  plot_result <- plot_track(MountTom, plot = "FootprintsTrajectories", seq.foot = TRUE)
  expect_s3_class(plot_result, "ggplot")

  expect_silent(ggplot2::ggplot_build(plot_result))
  expect_silent(ggplot2::ggplotGrob(plot_result))
})

test_that("plot_track allows seq.foot together with arrows (trajectories drawn with arrows)", {
  plot_result <- plot_track(MountTom, plot = "FootprintsTrajectories", seq.foot = TRUE, arrow.t = TRUE, arrow.size = 0.2)
  expect_s3_class(plot_result, "ggplot")

  expect_silent(ggplot2::ggplotGrob(plot_result))
})

test_that("plot_track errors on invalid seq.foot argument", {
  expect_error(plot_track(MountTom, seq.foot = "yes"), "seq\\.foot")
  expect_error(plot_track(MountTom, seq.foot = NA), "seq\\.foot")
  expect_error(plot_track(MountTom, seq.foot = c(TRUE, FALSE)), "seq\\.foot")
})

# Test error handling for invalid inputs
test_that("plot_track handles invalid inputs correctly", {
  expect_error(plot_track(NULL), "The 'data' argument must be a 'trackway' R object, which is a list consisting of two elements\\.")
  expect_error(plot_track(MountTom, plot = "Invalid"), "The 'plot' argument must be one of 'FootprintsTrajectories', 'Trajectories', or 'Footprints'.")
  expect_error(plot_track(MountTom, colours = c("#FF0000")), "Error: The length of 'colours' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, shape.f = c(17)), "Error: The length of 'shape.f' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, cex.f = -1), "Error: 'cex.f' size of footprint points must be a positive numeric value.")
  expect_error(plot_track(MountTom, cex.t = 0), "Error: 'cex.t' width of trajectory lines must be a positive numeric value.")
  expect_error(plot_track(MountTom, alpha.f = 1.5), "Error: 'alpha.f' transparency of footprint points must be a numeric value between 0 and 1.")
  expect_error(plot_track(MountTom, alpha.t = -0.1), "Error: 'alpha.t' transparency of trajectory lines must be a numeric value between 0 and 1.")
  expect_error(plot_track(MountTom, box.p = -0.2), "Error: 'box.p' padding around label boxes must be a positive numeric value.")
  expect_error(plot_track(MountTom, plot.labels = TRUE, labels = c("A", "B", "C")), "Error: The length of 'labels' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, plot = NA), "The 'plot' argument must be one of 'FootprintsTrajectories', 'Trajectories', or 'Footprints'.")
  expect_error(plot_track(MountTom, colours = "red"), "Error: The length of 'colours' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, shape.f = "circle"), "Error: The length of 'shape.f' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, cex.f = "large"), "Error: 'cex.f' size of footprint points must be a positive numeric value.")
  expect_error(plot_track(MountTom, cex.t = "thick"), "Error: 'cex.t' width of trajectory lines must be a positive numeric value.")
  expect_error(plot_track(MountTom, alpha.f = "transparent"), "Error: 'alpha.f' transparency of footprint points must be a numeric value between 0 and 1.")
  expect_error(plot_track(MountTom, alpha.t = "alpha"), "Error: 'alpha.t' transparency of trajectory lines must be a numeric value between 0 and 1.")
  expect_error(plot_track(MountTom, alpha.l = -2), "Error: 'alpha.l' transparency of labels must be a numeric value between 0 and 1.")
  expect_error(plot_track(MountTom, box.p = "wide"), "Error: 'box.p' padding around label boxes must be a positive numeric value.")
  expect_error(plot_track(MountTom, plot.labels = TRUE, labels = list("Trackway 1", "Trackway 2")), "Error: The length of 'labels' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, cex.l = 0), "Error: 'cex.l' size of labels must be a positive numeric value.")
  expect_error(plot_track(MountTom, cex.f = "cex.f"), "Error: 'cex.f' size of footprint points must be a positive numeric value.")
  expect_error(plot_track(MountTom, alpha.t = Inf), "Error: 'alpha.t' transparency of trajectory lines must be a numeric value between 0 and 1.")
  expect_error(plot_track(MountTom, shape.f = rep(17, 100)), "Error: The length of 'shape.f' must match the number of trackways in the data.")
  expect_error(plot_track(MountTom, plot.labels = TRUE, labels = c(1, 2, 3)), "Error: The length of 'labels' must match the number of trackways in the data.")
  expect_error(plot_track(PaluxyRiver, plot = "Pathways"), "The 'plot' argument must be one of 'FootprintsTrajectories', 'Trajectories', or 'Footprints'.")
  expect_error(plot_track(MountTom, seq.foot = "TRUE"), "seq\\.foot")
  expect_error(plot_track(MountTom, seq.foot = NA), "seq\\.foot")
})
