test_that("tps_to_track processes valid files correctly", {
  # Load an example TPS file
  tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")

  # Run function with no missing footprints (R.L.side now mandatory)
  result <- tps_to_track(
    tpsPaluxyRiver,
    scale = 0.004341493,
    R.L.side = c("R", "L"),
    missing = FALSE,
    NAs = NULL
  )

  # Check if output is a list
  expect_type(result, "list")

  # Check structure of Trajectories and Footprints
  expect_true(all(sapply(result$Trajectories, is.data.frame)))
  expect_true(all(sapply(result$Footprints,   is.data.frame)))

  # Footprints must include a Side column (R/L) now in all cases
  expect_true(all(vapply(result$Footprints, function(df) "Side" %in% names(df), logical(1))))
  expect_true(all(unlist(lapply(result$Footprints, function(df) all(df$Side %in% c("R","L"))))))
})


# Test for missing footprints handling
test_that("tps_to_track interpolates missing footprints correctly", {
  tpsMountTom <- system.file("extdata", "MountTom.tps", package = "QuAnTeTrack")

  # Here: first column = track, second = footprint (track 3, footprint 7)
  NAs <- matrix(c(3, 7), nrow = 1, ncol = 2)

  R.L.side <- c(
    "R", "L", "L", "L", "R", "L", "R", "R", "L", "L", "L",
    "L", "L", "R", "R", "L", "R", "R", "L", "R", "R",
    "R", "R"
  )

  result <- tps_to_track(
    tpsMountTom,
    scale = 0.004411765,
    R.L.side = R.L.side,
    missing = TRUE,
    NAs = NAs
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("Trajectories", "Footprints"))

  # Ensure missing footprints have been interpolated and marked in the correct track
  expect_true(any(result$Footprints[[3]]$missing == "Inferred"))

  # Footprints include Side column
  expect_true(all(vapply(result$Footprints, function(df) "Side" %in% names(df), logical(1))))
})


# Error handling tests
test_that("tps_to_track throws informative errors for incorrect inputs", {
  # Missing scale (still caught before file read)
  expect_error(
    tps_to_track("invalid.tps", scale = NULL, R.L.side = c("R"), missing = FALSE),
    "The 'scale' argument is missing"
  )

  # Missing R.L.side (now mandatory)
  expect_error(
    tps_to_track("invalid.tps", scale = 0.01, missing = FALSE),
    "The 'R.L.side' argument is mandatory"
  )
})

test_that("tps_to_track validates R.L.side length vs. number of tracks", {
  tpsPaluxyRiver <- system.file("extdata", "PaluxyRiver.tps", package = "QuAnTeTrack")

  # Provide wrong-length R.L.side vector to trigger length check
  expect_error(
    tps_to_track(
      tpsPaluxyRiver,
      scale = 0.004341493,
      R.L.side = c("R"),  # too short; Paluxy example has 2 tracks in docs
      missing = FALSE
    ),
    "Length of 'R.L.side'"
  )
})
