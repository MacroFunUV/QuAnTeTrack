test_that("anatomical_error_partitioning basic structure and contents", {
  skip_if_not_installed("trajr")

  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }

  vars <- c("Distance","Length","Straightness","Sinuosity","StLength")

  trk <- make_data()
  res <- anatomical_error_partitioning(trk,
                                       error_radius = 0.05,
                                       variables = vars,
                                       n_sim = 20,
                                       distribution = "uniform",
                                       seed = 999)

  expect_s3_class(res, "error_partitioning")
  expect_true(is.data.frame(res$summary))
  expect_true(is.data.frame(res$snr))
  expect_true(is.data.frame(res$qc))
  expect_true(is.data.frame(res$analysis_table))
  expect_true(all(c("variable","component","variance","percent","R2_c") %in% names(res$summary)))
  expect_true(all(res$summary$component %in% c("track","anatomical","Residual")))
  expect_equal(sort(unique(res$summary$variable)), sort(vars))
  expect_true(all(c("variable","bio_component","bio_var","error_var","SNR") %in% names(res$snr)))
  expect_true(all(res$snr$bio_component == "track"))
  expect_true(all(c("variable","SNR","SNR_rating","top_component","top_percent","R2_c","observer_estimable") %in% names(res$qc)))
  expect_true(all(res$qc$observer_estimable == "no (anatomical-only, sim-based)"))
  expect_equal(nrow(res$analysis_table), 20 * 2)
  expect_true(is.factor(res$analysis_table$track))
  expect_true(is.factor(res$analysis_table$anatomical))
  expect_equal(length(levels(res$analysis_table$anatomical)), 20)
})

test_that("anatomical_error_partitioning respects distribution and seed", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }
  vars <- c("Distance","Length","Straightness")

  trk <- make_data()
  a <- anatomical_error_partitioning(trk, 0.05, variables = vars, n_sim = 25, distribution = "gaussian", seed = 101)
  b <- anatomical_error_partitioning(trk, 0.05, variables = vars, n_sim = 25, distribution = "gaussian", seed = 101)
  c <- anatomical_error_partitioning(trk, 0.05, variables = vars, n_sim = 25, distribution = "uniform",  seed = 101)

  expect_equal(a$summary$variance, b$summary$variance, tolerance = 1e-12)
  av <- a$summary$variance[a$summary$component == "anatomical"]
  cv <- c$summary$variance[c$summary$component == "anatomical"]
  expect_false(isTRUE(all.equal(av, cv)))
})

test_that("error_radius handling: scalar recycling and per-track vector", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }
  vars <- c("Distance","Length")

  trk <- make_data()
  r1 <- anatomical_error_partitioning(trk, 0.05, variables = vars, n_sim = 15, distribution = "uniform", seed = 7)
  r2 <- anatomical_error_partitioning(trk, c(0.05, 0.05), variables = vars, n_sim = 15, distribution = "uniform", seed = 7)
  expect_equal(r1$summary$variance, r2$summary$variance, tolerance = 1e-12)

  r3 <- anatomical_error_partitioning(trk, c(0.02, 0.08), variables = vars, n_sim = 15, distribution = "uniform", seed = 7)
  a1 <- r1$summary$variance[r1$summary$component == "anatomical"]
  a3 <- r3$summary$variance[r3$summary$component == "anatomical"]
  expect_false(isTRUE(all.equal(a1, a3)))
})

test_that("angular variables are included and finite", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }

  trk <- make_data()
  res <- anatomical_error_partitioning(trk, 0.05,
                                       variables = c("TurnAng","PaceAng"),
                                       n_sim = 20, distribution = "gaussian", seed = 2025)
  expect_true(all(c("TurnAng","PaceAng") %in% unique(res$summary$variable)))
  sub <- subset(res$summary, component != "Residual")
  expect_true(all(is.finite(sub$variance)))
})

test_that("edge regimes: radius=0 gives zero anatomical; identical tracks give near-zero biological", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }

  trk <- make_data()
  res0 <- anatomical_error_partitioning(trk, 0, variables = c("Distance","Length"),
                                        n_sim = 10, distribution = "uniform", seed = 1)
  anat0 <- subset(res0$summary, component == "anatomical")
  expect_true(all(abs(anat0$variance) < 1e-12))

  trk2 <- make_data()
  trk2$Footprints$B  <- trk2$Footprints$A
  trk2$Trajectories$B <- trk2$Trajectories$A
  resI <- anatomical_error_partitioning(trk2, 0.05, variables = c("Distance","Length"),
                                        n_sim = 30, distribution = "gaussian", seed = 2)
  bioI  <- subset(resI$summary, component == "track")
  anatI <- subset(resI$summary, component == "anatomical")
  expect_true(all(anatI$variance > 0))
  expect_true(mean(bioI$variance, na.rm = TRUE) <= 0.1 * mean(anatI$variance, na.rm = TRUE))
})

test_that("percent sums ~100 and R2_c in [0,1]", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }

  vars <- c("Distance","Length","Straightness")
  trk <- make_data()
  res <- anatomical_error_partitioning(trk, 0.05, variables = vars,
                                       n_sim = 25, distribution = "gaussian", seed = 314)

  for (v in vars) {
    rows <- subset(res$summary, variable == v)
    expect_equal(sum(rows$percent, na.rm = TRUE), 100, tolerance = 1e-7)
  }
  expect_true(all(res$summary$R2_c >= 0 | is.na(res$summary$R2_c)))
  expect_true(all(res$summary$R2_c <= 1 | is.na(res$summary$R2_c)))
})

test_that("track_names preserved or sensibly defaulted", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }

  trk <- make_data()
  res <- anatomical_error_partitioning(trk, 0.05, variables = "Distance", n_sim = 10, distribution = "uniform", seed = 10)
  expect_equal(res$track_names, c("A","B"))

  trk2 <- trk
  names(trk2$Trajectories) <- c("", "")
  res2 <- anatomical_error_partitioning(trk2, 0.05, variables = "Distance", n_sim = 10, distribution = "uniform", seed = 10)
  expect_equal(res2$track_names, c("Track_01","Track_02"))
})

test_that("input validation errors and warnings follow current behavior", {
  skip_if_not_installed("trajr")
  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2, Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2, Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B), Trajectories = list(A = trajA, B = trajB))
  }
  trk <- make_data()

  expect_error(anatomical_error_partitioning("notalist", 0.05))
  expect_error(anatomical_error_partitioning(list(Footprints = trk$Footprints), 0.05))
  expect_error(anatomical_error_partitioning(list(Trajectories = trk$Trajectories), 0.05))
  expect_error(anatomical_error_partitioning(trk, NULL))
  expect_error(anatomical_error_partitioning(trk, c(0.05,0.05,0.05)))
  expect_error(anatomical_error_partitioning(trk, "abc"))

  # variables: mixed unknown warns but runs; all unknown -> error (no warning)
  expect_warning(
    anatomical_error_partitioning(trk, 0.05, variables = c("Distance","NopeVar"),
                                  n_sim = 5, distribution = "uniform", seed = 1)
  )
  expect_error(
    anatomical_error_partitioning(trk, 0.05, variables = c("Foo","Bar"),
                                  n_sim = 5, distribution = "uniform", seed = 1)
  )

  # n_sim invalid -> errors (character first warns, then errors)
  expect_error(anatomical_error_partitioning(trk, 0.05, n_sim = 0))
  expect_error(anatomical_error_partitioning(trk, 0.05, n_sim = -2))
  expect_warning(expect_error(anatomical_error_partitioning(trk, 0.05, n_sim = "ten")))

  # distribution invalid -> error (match.arg)
  expect_error(anatomical_error_partitioning(trk, 0.05, distribution = "weird"))

  # seed non-numeric -> error from set.seed(); suppress base coercion warning
  expect_error(suppressWarnings(anatomical_error_partitioning(trk, 0.05, n_sim = 5, seed = "seed")))
})

test_that("metadata: models is NULL and formula note present", {
  skip_if_not_installed("trajr")

  make_data <- function() {
    A <- data.frame(X = c(0.00, 0.50, 1.10, 1.60),
                    Y = c(0.00, 0.10, 0.35, 0.85),
                    IMAGE = "imgA", ID = "A", Side = c("L","R","L","R"))
    B <- data.frame(X = c(0.00, 0.80, 1.70, 2.70),
                    Y = c(0.00, 0.05, 0.10, 0.15),
                    IMAGE = "imgB", ID = "B", Side = c("L","R","L","R"))
    trajA <- trajr::TrajFromCoords(data.frame(X = (A$X[-1]+A$X[-4])/2,
                                              Y = (A$Y[-1]+A$Y[-4])/2))
    trajB <- trajr::TrajFromCoords(data.frame(X = (B$X[-1]+B$X[-4])/2,
                                              Y = (B$Y[-1]+B$Y[-4])/2))
    list(Footprints = list(A = A, B = B),
         Trajectories = list(A = trajA, B = trajB))
  }

  trk <- make_data()
  res <- anatomical_error_partitioning(trk, 0.05, variables = "Distance",
                                       n_sim = 6, distribution = "uniform", seed = 1)
  expect_null(res$models)
  expect_true(any(grepl("simulation-only decomposition via law of total variance",
                        unlist(res$formulae), fixed = TRUE)))
})
