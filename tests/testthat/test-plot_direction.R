test_that("plot_direction runs without errors for all plot types with MountTom dataset", {
  expect_no_error(plot_direction(MountTom, plot_type = "boxplot"))
  expect_no_error(plot_direction(MountTom, plot_type = "polar_steps"))
  expect_no_error(plot_direction(MountTom, plot_type = "polar_average"))
  expect_no_error(plot_direction(MountTom, plot_type = "faceted"))
})

test_that("plot_direction handles custom parameters correctly with MountTom dataset", {
  expect_no_error(plot_direction(MountTom, plot_type = "polar_steps", angle_range = 90))
  expect_no_error(plot_direction(MountTom,
    plot_type = "polar_steps", y_labels_position = 0,
    y_breaks_manual = c(0, 15, 30, 45, 60)
  ))
})

test_that("plot_direction runs without errors for all plot types with PaluxyRiver dataset", {
  expect_no_error(plot_direction(PaluxyRiver, plot_type = "boxplot"))
  expect_no_error(plot_direction(PaluxyRiver, plot_type = "polar_steps"))
  expect_no_error(plot_direction(PaluxyRiver, plot_type = "polar_average", y_breaks_manual = c(1, 2)))
  expect_no_error(plot_direction(PaluxyRiver, plot_type = "faceted"))
})

test_that("plot_direction handles custom parameters correctly with PaluxyRiver dataset", {
  expect_no_error(plot_direction(PaluxyRiver, plot_type = "polar_average", y_breaks_manual = c(1, 2)))
  expect_no_error(plot_direction(PaluxyRiver, plot_type = "polar_steps", y_labels_position = -90))
})

test_that("plot_direction handles invalid inputs correctly", {
  expect_error(plot_direction(NULL), "The 'data' argument must be a 'track' R object, which is a list consisting of two elements: 'Trajectories' and 'Footprints'.")
  expect_error(plot_direction(MountTom, plot_type = "invalid"), "The 'plot_type' must be one of 'boxplot', 'polar_steps', 'polar_average', or 'faceted'.")
  expect_error(plot_direction(MountTom, plot_type = "polar_steps", angle_range = -10), "The 'angle_range' must be a numeric value between 1 and 180 degrees.")
  expect_error(plot_direction(MountTom, plot_type = "polar_steps", y_labels_position = 200), "The 'y_labels_position' must be a numeric value between -180 and 180 degrees.")
})
