test_that("observer_error_partitioning: setup from extdata works", {
  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  expect_true(file.exists(tps))

  lines_tps <- readLines(tps, warn = FALSE)
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)

  track_list <- tps_to_track(
    file = tps, scale = 0.004341493, R.L.side = RL, missing = FALSE, NAs = NULL
  )

  expect_true(is.list(track_list))
  expect_true(length(track_list$Trajectories) == length(track_list$Footprints))
  expect_gt(length(track_list$Trajectories), 4)

  n <- length(track_list$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

  vars <- c("Distance","Straightness","TurnAng","PaceAng")

  res <- suppressWarnings(
    observer_error_partitioning(data = track_list, metadata = md, variables = vars)
  )

  expect_s3_class(res, "error_partitioning")
  expect_true(is.data.frame(res$summary))
  expect_true(is.data.frame(res$snr))
  expect_true(is.data.frame(res$qc))
  expect_true(is.data.frame(res$analysis_table))
  expect_true(is.list(res$formulae))
  expect_true(is.character(res$track_names) || is.null(res$track_names))

  expect_true(all(c("variable","component","variance","percent","R2_c") %in% names(res$summary)))
  expect_true(all(res$summary$component %in% c("track","observer","observer:track","Residual")))
  expect_equal(sort(unique(res$summary$variable)), sort(vars))

  expect_true(all(c("variable","bio_component","bio_var","error_var","SNR") %in% names(res$snr)))
  expect_true(all(res$snr$bio_component == "track"))

  expect_true(all(c("variable","SNR","SNR_rating","top_component",
                    "top_percent","R2_c","observer_estimable") %in% names(res$qc)))
  expect_true(any(res$qc$observer_estimable %in% c("yes","no (1 level)")))
})

test_that("observer_error_partitioning: includes interaction when estimable", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

  vars <- c("Distance","Length")
  res <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md, variables = vars)
  )

  comp_by_var <- split(res$summary$component, res$summary$variable)
  expect_true(any("observer:track" %in% comp_by_var[["Distance"]]))
})

test_that("observer_error_partitioning: intra-only and inter-only regimes behave", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

  vars <- c("Distance","Straightness")

  # Intra-only: force single observer level -> 'observer' & 'observer:track' drop
  md_intra <- md; md_intra$observer <- "obs1"
  res_intra <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md_intra, variables = vars)
  )
  comps_intra <- unique(res_intra$summary$component)
  expect_false("observer" %in% comps_intra)
  expect_false("observer:track" %in% comps_intra)
  expect_true("track" %in% comps_intra)

  # Inter-only: keep multiple observers but **collapse to one replicate per {observer,track}**
  grp <- interaction(md$observer, md$track, drop = TRUE)
  keep_idx <- ave(seq_len(n), grp, FUN = function(i) i[1])  # first row per cell
  keep <- sort(unique(keep_idx))

  trks_inter <- list(
    Trajectories = trks$Trajectories[keep],
    Footprints  = trks$Footprints[keep]
  )
  md_inter <- md[keep, , drop = FALSE]
  md_inter$replica <- 1L

  res_inter <- suppressWarnings(
    observer_error_partitioning(data = trks_inter, metadata = md_inter, variables = vars)
  )
  comps_inter <- unique(res_inter$summary$component)
  expect_true("observer" %in% comps_inter)
  expect_false("observer:track" %in% comps_inter)  # not estimable with 1 rep/cell
})

test_that("observer_error_partitioning: angular handling present and finite", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

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
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

  vars <- c("Distance","Length","Straightness")
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
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

  vars <- c("Distance","Length")
  res <- suppressWarnings(
    observer_error_partitioning(data = trks, metadata = md, variables = vars)
  )

  expect_true(all(c("replica","track","observer") %in% names(res$analysis_table)))
  fm <- res$formulae[["Distance"]]
  expect_true(!is.null(fm))
})

test_that("observer_error_partitioning: input validation and warnings", {
  skip_if_not_installed("trajr"); skip_if_not_installed("lme4"); skip_if_not_installed("MuMIn")

  tps <- system.file("extdata", "PaluxyRiverObsRep.tps", package = "QuAnTeTrack")
  lines_tps <- readLines(tps, warn = FALSE)
  n_tracks_file <- length(grep("^\\s*LM\\s*=\\s*", lines_tps))
  RL <- rep(c("R","L"), length.out = n_tracks_file)
  trks <- tps_to_track(tps, 0.004341493, RL, FALSE, NULL)

  n  <- length(trks$Trajectories)
  md <- data.frame(
    track    = rep(c("T01","T02"), length.out = n),
    observer = ifelse(seq_len(n) <= ceiling(n/2), "obs1", "obs2"),
    stringsAsFactors = FALSE
  )
  md$replica <- ave(seq_len(n),
                    interaction(md$observer, md$track, drop = TRUE),
                    FUN = seq_along)
  md <- md[, c("replica","track","observer")]

  # Missing required columns
  expect_error(observer_error_partitioning(trks, md[, c("track","observer")]),
               "must include columns")

  # Row mismatch
  md_bad <- md[-1, ]
  expect_error(observer_error_partitioning(trks, md_bad),
               "must have one row per track entry")

  # Mixed known/unknown -> warn, but run
  expect_warning(
    observer_error_partitioning(trks, md, variables = c("Distance","NopeVar"))
  )

  # All unknown -> error (and a warning about unknowns is acceptable)
  expect_error(
    observer_error_partitioning(trks, md, variables = c("Foo","Bar")),
    "No valid variables"
  )
})
