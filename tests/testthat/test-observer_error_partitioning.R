test_that("observer_error_partitioning: setup from extdata works", {
  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  expect_true(file.exists(tps))

  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)

  trackway_list <- tps_to_track(
    file = tps, scale = 0.004341493, R.L.side = RL, missing = FALSE, NAs = NULL
  )

  expect_true(is.list(trackway_list))
  expect_true(length(trackway_list$Trajectories) == length(trackway_list$Footprints))
  expect_gt(length(trackway_list$Trajectories), 4)

  n <- length(trackway_list$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  vars <- c("PathLen","Straightness","TurnAng","PaceAng")

  res <- suppressWarnings(
    observer_error_partitioning(data = trackway_list, metadata = md, variables = vars)
  )

  expect_s3_class(res, "error_partitioning")
  expect_true(is.data.frame(res$summary))
  expect_true(is.data.frame(res$snr))
  expect_true(is.data.frame(res$qc))
  expect_true(is.data.frame(res$analysis_table))
  expect_true(is.list(res$formulae))
  expect_true(is.character(res$trackway_names) || is.null(res$trackway_names))

  expect_true(all(c("variable","component","variance","percent","R2_c") %in% names(res$summary)))
  expect_true(all(res$summary$component %in% c("trackway","observer","observer:trackway","Residual")))
  expect_equal(sort(unique(res$summary$variable)), sort(vars))

  expect_true(all(c("variable","bio_component","bio_var","error_var","SNR") %in% names(res$snr)))
  expect_true(all(res$snr$bio_component == "trackway"))

  expect_true(all(c("variable","SNR","SNR_rating","top_component",
                    "top_percent","R2_c","observer_estimable") %in% names(res$qc)))
  expect_true(any(res$qc$observer_estimable %in% c("yes","no (1 level)")))
})

test_that("observer_error_partitioning: includes interaction when estimable", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  vars <- c("PathLen","BeelineLen")
  res <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md, variables = vars)
  )

  comp_by_var <- split(res$summary$component, res$summary$variable)
  expect_true("observer:trackway" %in% comp_by_var[["PathLen"]])
})

test_that("observer_error_partitioning: intra-only and inter-only regimes behave", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  vars <- c("PathLen","Straightness")

  # Intra-only: force single observer level -> 'observer' & 'observer:trackway' drop
  md_intra <- md
  md_intra$observer <- "obs1"
  res_intra <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md_intra, variables = vars)
  )
  comps_intra <- unique(res_intra$summary$component)
  expect_false("observer" %in% comps_intra)
  expect_false("observer:trackway" %in% comps_intra)
  expect_true("trackway" %in% comps_intra)

  # Inter-only: keep multiple observers but set one replicate per {observer,trackway}
  grp <- interaction(md$observer, md$trackway, drop = TRUE)
  keep_idx <- ave(seq_len(n), grp, FUN = function(i) i[1])  # first row per cell
  keep <- sort(unique(keep_idx))

  trks_inter <- list(
    Trajectories = trks$Trajectories[keep],
    Footprints   = trks$Footprints[keep]
  )
  md_inter <- md[keep, , drop = FALSE]
  md_inter$replica <- 1L

  res_inter <- suppressWarnings(
    observer_error_partitioning(data = trks_inter, metadata = md_inter, variables = vars)
  )
  comps_inter <- unique(res_inter$summary$component)
  expect_true("observer" %in% comps_inter)
  expect_false("observer:trackway" %in% comps_inter)  # not estimable when nlevels(obs:trackway) == nrow
})

test_that("observer_error_partitioning: angular handling present and finite", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  res <- suppressWarnings(
    observer_error_partitioning(
      data = trks, metadata = md, variables = c("TurnAng","PaceAng")
    )
  )
  expect_true(all(c("TurnAng","PaceAng") %in% unique(res$summary$variable)))

  sub <- subset(res$summary, component != "Residual")
  expect_true(all(is.finite(sub$variance)))

  expect_true(is.list(res$models$TurnAng$fits))
  expect_true(all(c("sin","cos") %in% names(res$models$TurnAng$fits)))
})

test_that("observer_error_partitioning: percent sums ~100 and R2_c in [0,1]", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  vars <- c("PathLen","BeelineLen","Straightness")
  res <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md, variables = vars)
  )

  for (v in vars) {
    rows <- subset(res$summary, variable == v)
    expect_equal(sum(rows$percent, na.rm = TRUE), 100, tolerance = 1e-7)
  }
  expect_true(all(res$summary$R2_c >= 0 | is.na(res$summary$R2_c)))
  expect_true(all(res$summary$R2_c <= 1 | is.na(res$summary$R2_c)))
})

test_that("observer_error_partitioning: analysis table and formula capture", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  vars <- c("PathLen","BeelineLen")
  res <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md, variables = vars)
  )

  expect_true(all(c("replica","trackway","observer") %in% names(res$analysis_table)))
  fm <- res$formulae[["PathLen"]]
  expect_true(!is.null(fm))
})

test_that("observer_error_partitioning: input validation and warnings", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_trackways_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_trackways_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    trackway    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(
    seq_len(n),
    interaction(md$observer, md$trackway, drop = TRUE),
    FUN = seq_along
  )
  md <- md[, c("replica","trackway","observer")]

  # Missing required columns (avoid default vars triggering Gauge requirement)
  expect_error(
    observer_error_partitioning(trks, md[, c("trackway","observer")], variables = "PathLen"),
    "must include columns"
  )

  # Row mismatch (avoid default vars triggering Gauge requirement)
  md_bad <- md[-1, ]
  expect_error(
    observer_error_partitioning(trks, md_bad, variables = "PathLen"),
    "must have one row per trackway entry"
  )

  # Mixed known/unknown -> warn, but run (also avoid Gauge requirement)
  expect_warning(
    observer_error_partitioning(trks, md, variables = c("PathLen","NopeVar"))
  )

  # All unknown -> error
  expect_error(
    observer_error_partitioning(trks, md, variables = c("Foo","Bar")),
    "No valid variables"
  )
})
