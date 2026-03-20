# tests/testthat/test-test_direction_circular.R

test_that("test_direction correctly performs Watson–Williams test on MountTom dataset", {
  suppressWarnings(result <- test_direction(MountTom, analysis = "Watson-Williams"))

  # Check for expected output structure
  expect_true("assumption_results" %in% names(result))
  expect_true("global_test" %in% names(result))
  expect_true("pairwise" %in% names(result))

  # Check that assumption results contain required elements
  expect_true(all(c("rayleigh", "kappa", "kappa_range", "kappa_ratio") %in% names(result$assumption_results)))

  # Check that the global test is an htest object
  expect_s3_class(result$global_test, "htest")

  # Check pairwise structure if present
  if (!is.null(result$pairwise) && nrow(result$pairwise) > 0) {
    expect_true(all(c("track1","track2","statistic","p_value","method","p_adj") %in% names(result$pairwise)))
  }
})

test_that("test_direction correctly performs Watson–Wheeler test with permutation (B=5) on MountTom dataset", {
  suppressWarnings(result <- test_direction(MountTom, analysis = "Watson-Wheeler",
                                            permutation = TRUE, B = 5, seed = 123))

  # Check for expected output structure
  expect_true("assumption_results" %in% names(result))
  expect_true("global_test" %in% names(result))
  expect_true("pairwise" %in% names(result))

  # Check that the global test is permutation-based
  expect_s3_class(result$global_test, "htest")
  expect_true(grepl("Watson-Wheeler \\(permutation, B=5\\)", result$global_test$method))

  # Check pairwise structure if present
  if (!is.null(result$pairwise) && nrow(result$pairwise) > 0) {
    expect_true(all(c("track1","track2","statistic","p_value","method","p_adj") %in% names(result$pairwise)))
  }
})

test_that("test_direction correctly performs Watson–Williams test on PaluxyRiver dataset", {
  suppressWarnings(result <- test_direction(PaluxyRiver, analysis = "Watson-Williams"))

  # Check for expected output structure
  expect_true("assumption_results" %in% names(result))
  expect_true("global_test" %in% names(result))
  expect_true("pairwise" %in% names(result))

  # Check that the global test is an htest object
  expect_s3_class(result$global_test, "htest")

  # Check pairwise structure if present
  if (!is.null(result$pairwise) && nrow(result$pairwise) > 0) {
    expect_true(all(c("track1","track2","statistic","p_value","method","p_adj") %in% names(result$pairwise)))
  }
})

test_that("test_direction correctly performs Watson–Wheeler test with permutation (B=5) on PaluxyRiver dataset", {
  suppressWarnings(result <- test_direction(PaluxyRiver, analysis = "Watson-Wheeler",
                                            permutation = TRUE, B = 5, seed = 123))

  # Check for expected output structure
  expect_true("assumption_results" %in% names(result))
  expect_true("global_test" %in% names(result))
  expect_true("pairwise" %in% names(result))

  # Check that the global test is permutation-based
  expect_s3_class(result$global_test, "htest")
  expect_true(grepl("Watson-Wheeler \\(permutation, B=5\\)", result$global_test$method))

  # Check pairwise structure if present
  if (!is.null(result$pairwise) && nrow(result$pairwise) > 0) {
    expect_true(all(c("track1","track2","statistic","p_value","method","p_adj") %in% names(result$pairwise)))
  }
})

test_that("test_direction gives an error for invalid analysis type", {
  expect_error(
    test_direction(MountTom, analysis = "InvalidMethod"),
    "Invalid 'analysis' argument. Choose 'Watson-Williams' \\(default\\) or 'Watson-Wheeler'."
  )
})

test_that("test_direction gives an error when data is not a valid track object", {
  expect_error(
    test_direction(NULL, analysis = "Watson-Williams"),
    "The 'data' argument must be a 'track' R object, which is a list consisting of two elements: 'Trajectories' and 'Footprints'."
  )
  expect_error(
    test_direction(list(1, 2, 3), analysis = "Watson-Williams"),
    "Both elements of 'data' must be lists. Ensure that 'Trajectories' and 'Footprints' are provided."
  )
})

test_that("test_direction emits expected warnings (removal + assumptions) when applicable", {
  warnings <- capture_warnings(test_direction(MountTom, analysis = "Watson-Williams"))

  # Check removal warnings
  expect_true(any(grepl("removed from the analysis due to having 3 or fewer footprints", warnings)))

  # Check at least one assumption warning is triggered
  expect_true(any(grepl("near-uniform directions \\(Rayleigh p > 0\\.05\\)", warnings)) ||
                any(grepl("unstable concentration \\(kappa\\)", warnings)) ||
                any(grepl("kappa.*heterogeneous.*ratio > 2", warnings)))
})

test_that("test_direction gives an error when there are not enough valid tracks", {
  small_track_data <- subset_track(MountTom, tracks = c(1)) # only one track
  expect_error(
    test_direction(small_track_data, analysis = "Watson-Williams"),
    "Not enough tracks with more than 3 footprints for meaningful analysis"
  )
})
