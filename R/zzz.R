.onLoad <- function(libname, pkgname) {
  # Specify the paths to the 'MountTom.rds' and 'PaluxyRiver.rds' files relative to the package's root directory
  mount_tom_file <- system.file("data", "MountTom.rds", package = pkgname)
  paluxy_river_file <- system.file("data", "PaluxyRiver.rds", package = pkgname)
  
  # Check if the 'MountTom.rds' file exists and load it into the environment
  if (file.exists(mount_tom_file)) {
    mount_tom <- readRDS(mount_tom_file)
    assign("MountTom", mount_tom, envir = parent.env(environment()))
  } else {
    warning("MountTom.rds file not found in package data directory.")
  }
  
  # Check if the 'PaluxyRiver.rds' file exists and load it into the environment
  if (file.exists(paluxy_river_file)) {
    paluxy_river <- readRDS(paluxy_river_file)
    assign("PaluxyRiver", paluxy_river, envir = parent.env(environment()))
  } else {
    warning("PaluxyRiver.rds file not found in package data directory.")
  }
}

