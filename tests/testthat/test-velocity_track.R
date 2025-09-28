test_that("velocity_track gives a deprecation warning", {
  # Define hip heights for MountTom dataset
  H_mounttom <- c(
    1.380, 1.404, 1.320, 1.736, 1.364, 1.432, 1.508, 1.768, 1.600,
    1.848, 1.532, 1.532, 0.760, 1.532, 1.688, 1.620, 0.636, 1.784,
    1.676, 1.872, 1.648, 1.760, 1.612
  )

  expect_warning(
    result <- velocity_track(MountTom, H = H_mounttom),
    "deprecated"
  )

  # Check output structure
  expect_true(is.list(result))
  expect_equal(length(result), length(MountTom[[1]]))

  # Check that each track has expected parameter names
  expected_names <- c(
    "Step_velocities", "Mean_velocity", "Standard_deviation_velocity",
    "Maximum_velocity", "Minimum_velocity", "Step_relative_stride",
    "Mean_relative_stride", "Standard_deviation_relative_stride",
    "Maximum_relative_stride", "Minimum_relative_stride"
  )

  for (i in seq_along(result)) {
    expect_true(is.list(result[[i]]))
    expect_named(result[[i]], expected_names)
  }
})

test_that("velocity_track gives a deprecation warning for PaluxyRiver dataset", {
  H_paluxyriver <- c(3.472, 2.200)

  expect_warning(
    result <- velocity_track(PaluxyRiver, H = H_paluxyriver),
    "deprecated"
  )

  # Check output structure
  expect_true(is.list(result))
  expect_equal(length(result), length(PaluxyRiver[[1]]))

  expected_names <- c(
    "Step_velocities", "Mean_velocity", "Standard_deviation_velocity",
    "Maximum_velocity", "Minimum_velocity", "Step_relative_stride",
    "Mean_relative_stride", "Standard_deviation_relative_stride",
    "Maximum_relative_stride", "Minimum_relative_stride"
  )

  for (i in seq_along(result)) {
    expect_true(is.list(result[[i]]))
    expect_named(result[[i]], expected_names)
  }
})

test_that("velocity_track handles different velocity calculation methods with deprecation warning", {
  H_paluxyriver <- c(3.472, 2.200)
  method_paluxyriver <- c("A", "B")

  expect_warning(
    result <- velocity_track(PaluxyRiver, H = H_paluxyriver, method = method_paluxyriver),
    "deprecated"
  )

  expect_true(is.list(result))
  expect_equal(length(result), length(PaluxyRiver[[1]]))

  expected_names <- c(
    "Step_velocities", "Mean_velocity", "Standard_deviation_velocity",
    "Maximum_velocity", "Minimum_velocity", "Step_relative_stride",
    "Mean_relative_stride", "Standard_deviation_relative_stride",
    "Maximum_relative_stride", "Minimum_relative_stride"
  )

  for (i in seq_along(result)) {
    expect_named(result[[i]], expected_names)
  }
})
