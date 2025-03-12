test_that("tps_to_track processes valid files correctly", {
  # Load an example TPS file
  tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")

  # Run function with no missing footprints
  result <- tps_to_track(tpsPaluxyRiver, scale = 0.004341493, missing = FALSE, NAs = NULL)

  # Check if output is a list
  expect_type(result, "list")
  expect_named(result, c("Trajectories", "Footprints"))

  # Check structure of Trajectories and Footprints
  expect_true(all(sapply(result$Trajectories, is.matrix)))
  expect_true(all(sapply(result$Footprints, is.data.frame)))
})


# Test for missing footprints handling

test_that("tps_to_track interpolates missing footprints correctly", {
  tpsMountTom <- system.file("extdata", "MountTom.tps", package = "QuAnTeTrack")
  NAs <- matrix(c(7, 3), nrow = 1, ncol = 2)
  R.L.side <- rep("L", 10) # Assume all tracks start with Left

  result <- tps_to_track(tpsMountTom, scale = 0.004411765, missing = TRUE, NAs = NAs, R.L.side = R.L.side)

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("Trajectories", "Footprints"))

  # Ensure missing footprints have been interpolated
  expect_true(any(result$Footprints[[3]]$inferred))
})


# Error handling tests

test_that("tps_to_track throws errors for incorrect inputs", {
  # Missing scale
  expect_error(tps_to_track("invalid.tps", scale = NULL, missing = FALSE),
               "The 'scale' argument is missing")

  # Invalid file
  expect_error(tps_to_track("nonexistent.tps", scale = 0.01, missing = FALSE),
               "File does not exist")

  # Missing footprints but no NAs provided
  expect_error(tps_to_track("valid.tps", scale = 0.01, missing = TRUE, NAs = NULL),
               "NAs must be provided if missing = TRUE")
})
