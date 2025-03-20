if (Sys.getenv("CI") == "true" || !interactive()) {
  rgl::rgl.useNULL(TRUE)
}
