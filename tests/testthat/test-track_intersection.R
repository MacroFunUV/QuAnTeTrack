test_that("track_intersection correctly calculates intersection metrics without testing", {
  M <- track_intersection(PaluxyRiver, test = FALSE)

  # Basic structure
  expect_true(is.data.frame(M))
  expect_equal(nrow(M), length(PaluxyRiver[[1]]))
  expect_equal(ncol(M), length(PaluxyRiver[[1]]))

  # Symmetry + NA diagonal
  MM <- as.matrix(M)
  expect_true(all(is.na(diag(MM))))
  expect_true(isTRUE(all.equal(MM, t(MM), check.attributes = FALSE, tolerance = 0)))
})

test_that("track_intersection correctly calculates intersection metrics with testing (None)", {
  set.seed(1)
  suppressWarnings(s1 <- simulate_track(PaluxyRiver, nsim = 10))
  res <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1, origin.permutation = "None")

  # Output pieces present
  expect_true(is.list(res))
  expect_true(all(c("Intersection_metric",
                    "Intersection_metric_p_values",
                    "Intersection_metric_p_values_combined",
                    "Intersection_metric_simulations") %in% names(res)))

  # Intersection metric symmetric with NA diag
  MM <- as.matrix(res$Intersection_metric)
  expect_true(all(is.na(diag(MM))))
  expect_true(isTRUE(all.equal(MM, t(MM), check.attributes = FALSE, tolerance = 0)))
})

test_that("track_intersection correctly calculates intersection metrics with testing (Min.Box)", {
  set.seed(2)
  suppressWarnings(s2 <- simulate_track(PaluxyRiver, nsim = 10))
  res <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s2, origin.permutation = "Min.Box")

  expect_true(is.list(res))
  expect_true("Intersection_metric" %in% names(res))
  expect_true("Intersection_metric_p_values" %in% names(res))
})

test_that("track_intersection correctly calculates intersection metrics with testing (Conv.Hull)", {
  set.seed(3)
  suppressWarnings(s3 <- simulate_track(PaluxyRiver, nsim = 10))
  res <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s3, origin.permutation = "Conv.Hull")

  expect_true(is.list(res))
  expect_true("Intersection_metric" %in% names(res))
})

test_that("track_intersection correctly calculates intersection metrics with custom coordinate permutation", {
  sbMountTom <- subset_track(MountTom, tracks = c(1, 2, 3, 4, 7, 8, 9, 13, 15, 16, 18))
  set.seed(4)
  suppressWarnings(s5 <- simulate_track(sbMountTom, nsim = 10))
  area_origin <- matrix(c(50, 5, 10, 5, 10, 20, 50, 20), ncol = 2, byrow = TRUE)

  res <- track_intersection(sbMountTom,
                            test = TRUE, H1 = "Higher", sim = s5,
                            origin.permutation = "Custom", custom.coord = area_origin)

  expect_true(is.list(res))
  expect_true("Intersection_metric" %in% names(res))
})

test_that("track_intersection gives an error for invalid data input", {
  expect_error(track_intersection(NULL, test = FALSE),
               "The 'data' argument must be a 'track' R object, which is a list consisting of two elements.")
})

test_that("track_intersection gives an error for invalid origin permutation method", {
  set.seed(5)
  suppressWarnings(simdat <- simulate_track(PaluxyRiver, nsim = 10))
  expect_error(
    track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = simdat, origin.permutation = "InvalidMethod"),
    "Invalid 'origin.permutation'. Valid options are: None, Min.Box, Conv.Hull, Custom"
  )
})

test_that("track_intersection gives an error when test = TRUE but no sim data is provided", {
  expect_error(track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower"),
               "A 'sim' argument must be provided when 'test' is TRUE.")
})

test_that("track_intersection errors when custom.coord is missing for 'Custom' permutation", {
  set.seed(6)
  suppressWarnings(simdat <- simulate_track(PaluxyRiver, nsim = 10))
  expect_error(
    track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = simdat, origin.permutation = "Custom"),
    "If 'origin.permutation' is set to 'Custom', the 'custom.coord' must be provided."
  )
})

test_that("track_intersection errors when custom.coord is not a matrix or dataframe", {
  set.seed(7)
  suppressWarnings(simdat <- simulate_track(PaluxyRiver, nsim = 10))
  expect_error(
    track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = simdat, origin.permutation = "Custom", custom.coord = list(1, 2, 3)),
    "The 'custom.coord' must be a matrix or a data frame."
  )
})

test_that("track_intersection errors when custom.coord does not have two columns", {
  set.seed(8)
  suppressWarnings(simdat <- simulate_track(PaluxyRiver, nsim = 10))
  expect_error(
    track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = simdat, origin.permutation = "Custom", custom.coord = matrix(1:9, ncol = 3)),
    "The 'custom.coord' must have exactly two columns."
  )
})

test_that("track_intersection returns BH-adjusted p-values (None)", {
  set.seed(9)
  suppressWarnings(s1 <- simulate_track(PaluxyRiver, nsim = 10))
  res <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1, origin.permutation = "None")

  expect_true("Intersection_metric_p_values_BH" %in% names(res))

  # Coerce to matrices for numeric ops
  PrawM <- as.matrix(res$Intersection_metric_p_values)
  PBHM  <- as.matrix(res$Intersection_metric_p_values_BH)
  MM    <- as.matrix(res$Intersection_metric)

  # Shape + names
  expect_equal(dim(PrawM), dim(PBHM))
  expect_equal(rownames(PrawM), rownames(PBHM))
  expect_equal(colnames(PrawM), colnames(PBHM))

  # Diagonals must be NA
  expect_true(all(is.na(diag(PrawM))))
  expect_true(all(is.na(diag(PBHM))))
  expect_true(all(is.na(diag(MM))))

  # Intersection metric symmetric
  expect_true(isTRUE(all.equal(MM, t(MM), check.attributes = FALSE, tolerance = 0)))

  # P-value matrices are lower-triangular; check only that part
  lt <- lower.tri(PrawM)
  expect_true(all(PrawM[lt] >= 0 & PrawM[lt] <= 1, na.rm = TRUE))
  expect_true(all(PBHM[lt]  >= 0 & PBHM[lt]  <= 1, na.rm = TRUE))

  # BH equals independent recomputation on the same set (lower triangle)
  padj_check <- p.adjust(PrawM[lt], method = "BH")
  expect_equal(PBHM[lt], padj_check, tolerance = 1e-12)

  # Combined p-value in [0,1]
  expect_true(is.numeric(res$Intersection_metric_p_values_combined))
  expect_true(res$Intersection_metric_p_values_combined >= 0 &&
                res$Intersection_metric_p_values_combined <= 1)
})

test_that("track_intersection BH-adjusted p-values (Higher, Min.Box)", {
  set.seed(10)
  suppressWarnings(s2 <- simulate_track(PaluxyRiver, nsim = 10))
  res <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Higher", sim = s2, origin.permutation = "Min.Box")

  expect_true("Intersection_metric_p_values_BH" %in% names(res))

  PrawM <- as.matrix(res$Intersection_metric_p_values)
  PBHM  <- as.matrix(res$Intersection_metric_p_values_BH)

  lt <- lower.tri(PrawM)
  expect_true(all(PrawM[lt] >= 0 & PrawM[lt] <= 1, na.rm = TRUE))
  expect_true(all(PBHM[lt]  >= 0 & PBHM[lt]  <= 1, na.rm = TRUE))

  # Exact BH recomputation check (lower triangle only)
  padj_check <- p.adjust(PrawM[lt], method = "BH")
  expect_equal(PBHM[lt], padj_check, tolerance = 1e-12)
})

test_that("track_intersection (test=FALSE) returns symmetric data.frame with NA diagonal", {
  M <- track_intersection(PaluxyRiver, test = FALSE)
  expect_true(is.data.frame(M))
  expect_equal(nrow(M), length(PaluxyRiver[[1]]))
  expect_equal(ncol(M), length(PaluxyRiver[[1]]))

  MM <- as.matrix(M)
  expect_true(all(is.na(diag(MM))))
  expect_true(isTRUE(all.equal(MM, t(MM), check.attributes = FALSE, tolerance = 0)))
})

test_that("p-value outputs are data.frames with NA diagonal", {
  set.seed(11)
  suppressWarnings(s1 <- simulate_track(PaluxyRiver, nsim = 5))
  res <- track_intersection(PaluxyRiver, test = TRUE, H1 = "Lower", sim = s1, origin.permutation = "None")
  expect_true(is.data.frame(res$Intersection_metric_p_values))
  expect_true(is.data.frame(res$Intersection_metric_p_values_BH))
  expect_true(all(is.na(diag(as.matrix(res$Intersection_metric_p_values)))))
  expect_true(all(is.na(diag(as.matrix(res$Intersection_metric_p_values_BH)))))
})
