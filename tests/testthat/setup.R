<<<<<<< HEAD
if (requireNamespace("rgl", quietly = TRUE)) {
  rgl::rgl.useNULL()  # Ensure rgl works on headless systems (e.g., Linux CI)
=======
if (Sys.getenv("CI") == "true" || !interactive()) {
  rgl::rgl.useNULL()
>>>>>>> 2a23fe4079698aaf6c39532ab5fc23dca3674457
}
