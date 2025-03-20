if (requireNamespace("rgl", quietly = TRUE)) {
  rgl::rgl.useNULL()  # Ensure rgl works on headless systems (e.g., Linux CI)
}
