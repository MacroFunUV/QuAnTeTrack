test_that("cluster_track works correctly with MountTom and PaluxyRiver", {

  ## ---------------------------
  ## Helper: validate cluster labels (hclust or mclust)
  ## ---------------------------
  expect_class_labels <- function(x, n) {
    expect_true(is.numeric(x))
    expect_equal(length(x), n)
    expect_true(all(is.finite(x)))
    expect_true(all(x >= 1))
    expect_true(all(x == floor(x)))
  }

  ## ---------------------------
  ## MountTom dataset
  ## ---------------------------
  H_mounttom <- c(
    1.380, 1.404, 1.320, 1.736, 1.364, 1.432, 1.508, 1.768, 1.600,
    1.848, 1.532, 1.532, 0.760, 1.532, 1.688, 1.620, 0.636, 1.784,
    1.676, 1.872, 1.648, 1.760, 1.612
  )
  veltrack_MountTom <- velocity_track(MountTom, H = H_mounttom)

  expect_warning(
    cluster_track(
      MountTom,
      veltrack_MountTom,
      variables = c("TurnAng", "Velocity")
    ),
    "12 tracks were discarded for having fewer than 4 footprints. Discarded track indices: 5, 6, 10, 11, 12, 14, 17, 19, 20, 21, 22, 23"
  )

  ## hclust (default)
  result1 <- suppressWarnings(cluster_track(
    MountTom,
    veltrack_MountTom,
    variables = c("TurnAng", "Velocity")
  ))
  expect_class_labels(result1$clust$classification, nrow(result1$matrix))

  ## mclust (explicit)
  result2 <- suppressWarnings(cluster_track(
    MountTom,
    veltrack_MountTom,
    variables = c("Sinuosity", "StLength"),
    analysis = "mclust"
  ))
  expect_s3_class(result2$clust, "Mclust")
  expect_class_labels(result2$clust$classification, nrow(result2$clust$data))

  ## hclust
  result3 <- suppressWarnings(cluster_track(
    MountTom,
    veltrack_MountTom,
    variables = c("MaxVelocity", "MinVelocity")
  ))
  expect_class_labels(result3$clust$classification, nrow(result3$matrix))

  result4 <- suppressWarnings(cluster_track(
    MountTom,
    veltrack_MountTom,
    variables = "Straightness"
  ))
  expect_class_labels(result4$clust$classification, nrow(result4$matrix))

  ## ---------------------------
  ## PaluxyRiver dataset
  ## ---------------------------
  H_paluxyriver <- c(3.472, 2.200)
  Method_paluxyriver <- c("A", "B")
  veltrack_PaluxyRiver <- velocity_track(
    PaluxyRiver,
    H = H_paluxyriver,
    method = Method_paluxyriver
  )

  result5 <- cluster_track(
    PaluxyRiver,
    veltrack_PaluxyRiver,
    variables = c("Distance", "Straightness")
  )
  expect_true(is.data.frame(result5$matrix))
  expect_class_labels(result5$clust$classification, nrow(result5$matrix))

  result6 <- cluster_track(
    PaluxyRiver,
    veltrack_PaluxyRiver,
    variables = c("Length", "sdVelocity")
  )
  expect_class_labels(result6$clust$classification, nrow(result6$matrix))

  result7 <- cluster_track(
    PaluxyRiver,
    veltrack_PaluxyRiver,
    variables = c("TurnAng", "sdTurnAng")
  )
  expect_class_labels(result7$clust$classification, nrow(result7$matrix))

  result8 <- cluster_track(
    PaluxyRiver,
    veltrack_PaluxyRiver,
    variables = c("Sinuosity")
  )
  expect_class_labels(result8$clust$classification, nrow(result8$matrix))
})
