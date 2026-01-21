test_that("subset_track correctly subsets the first three trackways of MountTom dataset", {
  subset_data <- subset_track(MountTom, tracks = c(1:3))

  # Check that the number of trackways is correct
  expect_equal(length(subset_data[[1]]), 3)
  expect_equal(length(subset_data[[2]]), 3)
})

test_that("subset_track correctly subsets specific trackways by indices", {
  subset_data <- subset_track(MountTom, tracks = c(5, 7, 10))

  # Check that the correct number of trackways is returned
  expect_equal(length(subset_data[[1]]), 3)
  expect_equal(length(subset_data[[2]]), 3)
})

test_that("subset_track defaults to all trackways when 'tracks' is NULL", {
  subset_data <- subset_track(MountTom, tracks = NULL)

  # Check that all tracks are included
  expect_equal(length(subset_data[[1]]), length(MountTom[[1]]))
  expect_equal(length(subset_data[[2]]), length(MountTom[[2]]))
})

test_that("subset_track handles out-of-bounds trackway indices gracefully", {
  expect_warning(
    subset_data <- subset_track(MountTom, tracks = c(1, 2, 9999)),
    "The 'tracks' argument contains indices that exceed the length of the data. These indices will be ignored."
  )

  # Check that only the valid trackways are returned
  expect_equal(length(subset_data[[1]]), 2)
  expect_equal(length(subset_data[[2]]), 2)
})

test_that("subset_track returns an error for invalid trackway indices", {
  expect_error(
    subset_track(MountTom, tracks = c(-1, 0, "a")),
    "The 'tracks' argument must be a numeric vector of positive indices."
  )
})

test_that("subset_track correctly subsets trackways for PaluxyRiver dataset", {
  subset_data <- subset_track(PaluxyRiver, tracks = c(1, 2))

  # Check that the correct number of trackways is returned
  expect_equal(length(subset_data[[1]]), 2)
  expect_equal(length(subset_data[[2]]), 2)
})

test_that("subset_track returns an error for invalid data input", {
  expect_error(
    subset_track(NULL, tracks = c(1, 2)),
    "The 'data' argument must be a 'trackway' R object, which is a list consisting of two elements."
  )
  expect_error(
    subset_track(list(1, 2, 3), tracks = c(1, 2)),
    "The two elements of 'data' must be lists."
  )
})
