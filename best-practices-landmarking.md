# Introduction

This vignette complements the main QuAnTeTrack vignette by providing
practical guidance for quantifying, controlling, and reporting
uncertainty throughout the landmarking-to-analysis pipeline. The focus
is on decisions that directly affect computed trajectories and movement
parameters (e.g., turning angles, step lengths, sinuosity, velocity) and
on transparent practices that make results reproducible and comparable.

While QuAnTeTrack implements analysis, simulation, and visualization
tools, the automation of uncertainty quantification and propagation
(e.g., automated repeat-digitization audits or error-aware resampling)
is outside the scope of the current package. Nonetheless, this document
outlines concrete procedures that users can follow before running
downstream analyses. If community feedback indicates this would help new
users, an automated uncertainty-quantification workflow can be
considered for future versions.

# Why Uncertainty Matters and Sources of Uncertainty

Objective, quantitative methods for trackway analysis do not imply false
precision if the inherent uncertainty is explicitly measured, reported,
and managed. Landmarking and geometric morphometrics have proven highly
useful across biology and paleontology, but their reliability depends on
transparent handling of error. Without careful consideration,
uncertainty can accumulate from multiple sources and bias biological
interpretations. Below, the main sources of uncertainty are outlined,
together with examples of how they may affect results.

- **Anatomical definition of landmarks**  
  If the anatomical basis for placing a landmark is vague or
  inconsistently applied, different samplers may digitize slightly
  different points. For example, defining a landmark at the “tip of the
  digit” may be ambiguous in a poorly preserved track where the digit
  outline is rounded. Such ambiguities inflate variance and reduce
  reproducibility.

- **Sampler/observer effects**  
  Experience, training, and even fatigue of the sampler influence
  results. Intra-observer error (the same person digitizing multiple
  times) can lead to small but systematic differences, while
  inter-observer error (different people digitizing the same footprint)
  may reveal much larger discrepancies. Both affect the precision of
  derived metrics such as step length or pace angulation.

- **Preservation and sediment effects**  
  Substrate conditions strongly influence footprint shape. Erosion can
  erase edges, infill can distort depth, and plastic deformation can
  shift toe positions. Landmarks digitized in these areas may not
  reflect true foot anatomy, and comparisons across trackways with
  different preservation states can be misleading if not carefully
  controlled.

- **Image acquisition and scaling**  
  Digital recording methods—photographs, photogrammetry, laser
  scanning—introduce their own biases. Lens distortion, uneven
  resolution, or shadows can obscure features. If scale bars are tilted
  or not aligned with the footprint plane, all subsequent measurements
  (distances, velocities) are biased. Even small scaling errors
  propagate through the entire pipeline.

- **Georeferencing and orientation**  
  When trackways are digitized from field images or slabs, orientation
  relative to north, slope, or bedding planes must be set. Misalignment
  alters calculated movement directions and turning angles. Inconsistent
  or incorrect reference systems can make comparisons across trackways
  invalid.

- **Interpolation of missing footprints**  
  When parts of a trackway are reconstructed (e.g., a missing footprint
  inferred by stride length), uncertainty arises because the position is
  model-based. Overconfidence in interpolated points can bias trajectory
  smoothness or sinuosity measures.

- **Coordinate handling and unit conversion**  
  Errors in coordinate systems (e.g., mixing cm and mm, or incorrect CRS
  transformations) propagate into all derived parameters. Even rounding
  during resampling may bias small-scale metrics such as step
  variability.

These uncertainties highlight the importance of integrating error
evaluation into the workflow rather than treating digitized landmarks as
exact data.

# Best Practices for Dealing with Uncertainty in Landmarking

## Anatomical definition of landmarks

- Provide a precise, objective, replicable anatomical definition for
  each landmark.  
- Maintain written protocols and illustrated guides.  
- Run multiple rounds of digitization by the same sampler to evaluate
  internal consistency.

## Sampler and experience

- Involve multiple independent samplers with varied expertise.  
- Each sampler performs repeat digitization rounds to quantify
  inter-observer variability.  
- Use training/calibration sessions to harmonize criteria and reduce
  bias.

## Preservation and sediment effects

- Identify zones most affected by erosion, deformation, or infill.  
- Prioritize best-preserved anatomical regions when possible.  
- Where ambiguity exists, repeat digitization on alternative anatomical
  locations (e.g., heel vs. toe) to assess robustness.  
- Simulate uncertainty by randomly jittering landmarks within a
  realistic tolerance envelope (e.g., ±2–5 mm) and re-running analyses.

## Image acquisition and scaling

- Prefer high-resolution, calibrated methods (controlled photography,
  photogrammetry, laser scanning).  
- Ensure scale bars are in-plane with the surface and clearly visible.  
- Correct optical distortion; avoid oblique angles.  
- Regularly verify scale against known distances on the slab.

## Georeferencing and orientation of trackways

- Establish orientation relative to external references (north, bedding
  planes, slope).  
- Apply consistent coordinate systems across trackways to make turning
  angles comparable.  
- Where orientation is uncertain, digitize alternative alignments and
  test sensitivity of derived metrics.

## Data management and documentation

- Keep a digitization log (dates, observers, devices, versions,
  settings).  
- Save raw images/models, `.TPS` files, and scaling notes; record any
  interpolations done.

# Recommended Workflow with QuAnTeTrack (for quantifying and controlling uncertainty)

1.  **Protocol definition**
    - Finalize anatomical criteria and imaging protocols (camera/scanner
      setup, scale placement, lighting).
2.  **Digitization rounds**
    - Perform repeated digitization of the same footprints by the same
      sampler (intra-observer error).  
    - Repeat with multiple independent samplers (inter-observer error).
3.  **Preservation assessment**
    - Identify best-preserved footprints for primary analyses.  
    - Where preservation is poor, create multiple plausible landmark
      sets (e.g., toe-focused vs. heel-focused) and compare outcomes.
4.  **Image and scaling verification**
    - Validate scaling against known in-field/slab distances.  
    - From photographs, test distortion by re-digitizing scale-bar
      points and comparing results.
5.  **Orientation and georeferencing**
    - Standardize all trackways to a common orientation framework.  
    - When orientation is uncertain, run parallel analyses using
      alternative plausible orientations and compare impacts on turning
      angles and trajectories.
6.  **Uncertainty simulation**
    - Apply random jitter to landmark coordinates within a plausible
      error range (e.g., ±2–5 mm); re-run parameter extraction and
      analyses.  
    - Compare standard deviations of outputs (sinuosity, turning angles,
      velocities) across replicates to quantify robustness.
7.  **Integration into analyses**
    - Use averaged parameters (means across repeated rounds) in the main
      analysis.  
    - Report variability (SD, range, or confidence bounds) alongside
      point estimates.
8.  **Transparent reporting**
    - Document who digitized what, when, and how, and provide
      uncertainty summaries in Methods/Appendices.  
    - When possible, share raw data (images or 3D meshes), `.TPS`, and
      protocol files.

# How This Relates to the Package Documentation

The help page for
[`tps_to_track()`](https://macrofunuv.github.io/QuAnTeTrack/reference/tps_to_track.md)
includes a detailed checklist of potential uncertainty sources and
corresponding mitigation/quantification steps, aligned with the best
practices above.

Throughout the vignette examples (simulation, similarity metrics,
intersections, clustering), users are encouraged to repeat analyses over
replicated or jittered datasets and to report dispersion of results.

# Notes on Downstream Analyses

- **Direction and velocity tests**: run on repeated/jittered
  digitizations; compare the spread of p-values/effect sizes.  
- **Similarity metrics (DTW/Fréchet)**: test sensitivity to alignment
  choices (origin/centroid) and replicate datasets.  
- **Intersection counts**: evaluate whether conclusions (greater/fewer
  intersections than expected) are stable under scaling/orientation
  variants and jitter.  
- **Clustering**: when using
  [`cluster_track()`](https://macrofunuv.github.io/QuAnTeTrack/reference/cluster_track.md),
  consider running the pipeline on replicate matrices (from
  repeated/jittered digitizations) and summarizing cluster stability
  (e.g., co-assignment frequencies).

# Limitations and Future Directions

Automated, end-to-end uncertainty audits (e.g., built-in
repeated-digitization managers, automated jitter-resampling with summary
reports) are not part of the current release. If users and reviewers
consider such tooling beneficial for non-specialists, a future version
could include utilities to automate repeated-digitization aggregation,
jitter-based sensitivity analyses, and standardized uncertainty reports.

# References and Selected Reading

- Zelditch, M. L., Swiderski, D. L., & Sheets, H. D. (2021). *Geometric
  Morphometrics for Biologists: A Primer* (3rd ed.). Academic Press.  
- Rohlf, F. J., & Slice, D. (1990). Extensions of the Procrustes method
  for the optimal superimposition of landmarks. *Systematic Zoology,
  39*(1), 40–59.  
- Klingenberg, C. P., & McIntyre, G. S. (1998). Geometric morphometrics
  of developmental instability: analyzing patterns of fluctuating
  asymmetry with Procrustes methods. *Evolution, 52*(5), 1363–1375.  
- Bookstein, F. L. (1991). *Morphometric Tools for Landmark Data*.
  Cambridge University Press.  
- Dryden, I. L., & Mardia, K. V. (2016). *Statistical Shape Analysis*
  (2nd ed.). Wiley.  
- Slice, D. E. (2007). Geometric morphometrics. *Annual Review of
  Anthropology, 36*, 261–281.  
- Batschelet, E. (1981). *Circular Statistics in Biology*. Academic
  Press.  
- Benhamou, S. (2004). How to reliably estimate the tortuosity of an
  animal’s path: straightness, sinuosity, or fractal dimension? *Journal
  of Theoretical Biology, 229*(2), 209–220.  
- Alexander, R. M. (1976). Estimates of speeds of dinosaurs. *Nature,
  261*(5556), 129–130.  
- Ruiz, J., & Torices, A. (2013). Humans running at stadiums and beaches
  and the accuracy of speed estimations from fossil trackways. *Ichnos,
  20*(1), 31–35.

**Additional references**  
- Fruciano, C. (2016). Measurement error in geometric morphometrics.
*Development Genes and Evolution, 226*(3), 139–158.
<https://doi.org/10.1007/s00427-016-0537-4>  
- Robinson, C., & Terhune, C. E. (2017). Error in landmark-based
geometric morphometrics: Testing the assumption of homology in fossils.
*American Journal of Physical Anthropology, 132*(2), 301–310.
<https://doi.org/10.1002/ajpa.20616>

# Acknowledgments

These guidelines were shaped by feedback from reviewers and users
seeking clear, field-ready practices for managing uncertainty from
digitization to inference. They are intended to help users implement
robust pipelines with transparent error reporting in QuAnTeTrack.
