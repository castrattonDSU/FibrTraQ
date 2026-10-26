# FibrTraQ Publication Statistical Testing Plan

This plan separates two different claims that should not be mixed in a
manuscript:

1. Software validation: FibrTraQ gives the correct answer when the expected
   tract/ROI relationship is known.
2. Biological inference: a tractography-derived connectivity pattern differs
   across regions, subjects, groups, datasets, or hypotheses.

The current validation evidence can support the software-validation claim. A
strong biological claim needs subject-level replication or a clearly defined
spatial-null/permutation design.

## Publication Claims To Support

Recommended primary claims:

- FibrTraQ performs exact voxel-resolved tract-to-target counting from raw
  tractograms and ROI volumes.
- The workflow preserves important non-target information through
  `leave_no_target`, so atlas coverage failures and unresolved exits are visible
  instead of silently discarded.
- The software generalizes from simple one-bundle known-answer tests to a
  many-target exploratory analysis in a large MNI-space normative connectome.
- FibrTraQ distinguishes full-path target intersection from anatomically
  stronger paired focal-target endpoint connectivity and quantifies the number
  of exact target voxels visited by each target-intersection-positive
  streamline.
- Exploratory figures can localize focal-voxel connectivity patterns and show
  sampled real streamlines in 2D and 3D anatomical context.

Do not claim disease, developmental, treatment, or population biology from a
single tractogram. Use the single-tractogram dTOR analysis as a utility and
face-validity demonstration.

## Software Validation Tests

### 1. Known-Answer Intersection And Endpoint Recovery

Use ISMRM 2015 bundle/sub-bundle ground-truth tractograms and ROI masks.

Unit of analysis:

- One bundle/sub-bundle run per tractogram submission or ground-truth bundle.

Primary metrics:

- Full-path target-intersection count.
- Total focal-hit tracks.
- Target-intersection rate: `target_intersection_tracks / focal_hit_tracks`.
- Exact paired focal-target endpoint count and rate among full-path target
  intersections.
- Exact paired focal-target endpoint rate among focal-endpoint streamlines.
- One-voxel paired endpoint sensitivity, explicitly labeled as 26-neighbor
  voxel-grid sensitivity rather than a physical-distance tolerance.
- Target penetration, measured as the exact number of unique target voxels
  visited by each full-path target intersection.
- `leave_no_target` count and fraction.
- False-target rate when non-matching targets are included.
- Reproducibility: identical output tables across repeated runs.

Recommended tests:

- Exact known-answer comparison for the ground-truth bundle manifest.
- Bootstrap confidence intervals across bundle/sub-bundle units.
- Wrong-target negative control: paired comparison of correct-target hit rate
  versus anatomically wrong target hit rate.
- Shifted-mask negative control: paired comparison of correct target/focal ROIs
  versus translated ROIs.
- Exact reconciliation of the five mutually exclusive endpoint classes to each
  target's full-path intersection count.
- Exact reconciliation of each sparse penetration histogram to its target's
  full-path intersection count.

Recommended reporting:

- Median and range of correct target-hit rate across bundles.
- Median and range of wrong-target false-positive rate.
- Number of manifests passing preflight, run, and validation.
- Exact reproducibility status from rerun output equality.
- Paired endpoint percentages with the denominator named explicitly; do not
  combine percentages of target intersections with percentages of
  focal-endpoint streamlines.
- One-, two-, three-, five-, and ten-voxel target-penetration thresholds plus
  the 25th, 50th, 75th, and 95th percentiles.

### 2. Negative Controls

Controls should be built from the same input files so only the biological
relationship changes.

Controls:

- Wrong endpoint: keep the focal ROI, replace the target with a different
  anatomical endpoint.
- Shifted endpoint: translate the endpoint mask by a fixed voxel offset that
  preserves mask size but breaks anatomy.
- Shifted focal ROI: translate the focal mask away from the true bundle entry.
- Left-right swapped label where anatomy makes this inappropriate.

Recommended tests:

- Paired Wilcoxon signed-rank test or permutation sign-flip test across bundle
  units: correct hit rate versus wrong/shifted hit rate.
- Report effect sizes as paired median difference and bootstrap confidence
  interval.
- Report false positive rate under each negative control.

Interpretation:

- The correct condition should have materially higher target-hit rate than the
  negative controls.
- `leave_no_target` should increase when target masks are wrong or incomplete.
- These tests validate counting and ROI handling; they do not validate that any
  one tractography method reconstructs biological axons.

### 3. Coordinate and Discretization Robustness

Run the same manifest under accepted coordinate settings.

Checks:

- `track_space = "trackvis"` for TrackVis `.trk` files with valid header affine.
- `track_space = "world"` only when streamlines are already in image-world
  coordinates.
- Compare `round`, `floor`, `ceil`, and `any` discretization modes when the
  scientific question requires sensitivity analysis.

Recommended tests:

- Equivalence/tolerance check on total focal-hit and target-hit counts.
- Confirm that validated coordinate mode is stable across reruns.
- Confirm that intentionally wrong coordinate settings fail preflight or produce
  clear negative-control behavior.

### 4. Deterministic Reproducibility

Run the same manifest twice from raw `.trk` and ROI files.

Required result:

- Per-voxel output tables are identical.
- RSC summaries are identical.
- Endpoint summary, penetration, and context tables are identical when endpoint
  analysis is enabled.
- Spatial maps are regenerated without missing panels or color-scale
  inconsistencies.

Recommended reporting:

- Hashes or exact equality status for output CSV files.
- Runtime and peak-memory notes for streamed versus non-streamed runs.

### 5. Biological Context And Provenance

Endpoint pairing strengthens interpretation but does not by itself establish
that tractogram coordinates and anatomy share a valid biological reference.
Report the separate `biological_qc_status` as `pass`, `review`, or
`unavailable` and preserve the reason for any non-pass status.

Required reporting:

- Reference-grid endpoint counts, including per-axis low/high boundary exits
  and separation of one-voxel from greater-than-one-voxel outliers.
- Endpoint containment for every externally validated anatomical context mask,
  including the mask path, checksum, dimensions, affine, and voxel count.
- Focal-target mask overlap and shared-mask endpoint ambiguity.
- The exact tractogram reference geometry and provenance when available.

Do not infer an anatomical mask by thresholding an intensity image, and do not
adopt a formal context-containment pass threshold without an external
scientific justification.

## Biological Hypothesis Tests

These are appropriate when the dataset contains multiple independent subjects,
independent tractograms, or a prespecified spatial-null framework.

### 1. Ipsilateral Enrichment

Hypothesis:

- A left focal ROI has more left-hemisphere than right-hemisphere target associations
  after excluding unlateralized targets and `leave_no_target`.

Subject-level metric:

- `logit((left_target_hits + 0.5) / (left_target_hits + right_target_hits + 1))`
  or a left/right normalized proportion.

Recommended tests:

- One-sample test against 0.5 using subject-level values.
- Mixed-effects model with subject as the independent unit and covariates such
  as age, sex, site, tractogram size, and scanner/protocol when available.
- Non-parametric subject-level bootstrap confidence interval.

### 2. Focal-Voxel Topography

Hypothesis:

- Target specificity varies systematically across the focal ROI, for example
  anterior/posterior or medial/lateral gradients.

Voxel-level metrics:

- Target fraction per focal voxel.
- Most-supported target class per focal voxel.
- Shannon entropy of target distribution per focal voxel.
- `leave_no_target` fraction per focal voxel.

Recommended tests:

- Subject-level regression slopes of target fraction on MNI x/y/z coordinates.
- Mixed-effects model: `target_fraction ~ coordinate + (1 | subject)`.
- Cluster or permutation correction over focal voxels when producing maps.
- False discovery rate correction across target classes.

Interpretation:

- This makes FibrTraQ useful for testing spatial organization inside an ROI,
  not just whether a region connects somewhere.

### 3. Target Specificity and Diversity

Hypothesis:

- Some focal voxels have concentrated target support while others have broader
  target diversity.

Metrics:

- Entropy per focal voxel.
- Number of nonzero target classes per focal voxel.
- Top-target fraction per focal voxel.
- ROI-level median entropy per subject.

Recommended tests:

- Paired ROI comparison when two focal ROIs are available in each subject.
- Group comparison using subject-level summaries.
- Mixed-effects model if multiple ROIs, sites, or acquisition protocols are
  included.

### 4. Atlas Coverage and `leave_no_target`

Hypothesis:

- The fraction of focal-exiting streamlines without a target reflects atlas
  coverage, target definition, or tractogram method differences.

Metrics:

- ROI-level `leave_no_target` fraction.
- Per-voxel `leave_no_target` fraction map.
- Ratio of `leave_no_target` to total focal-hit tracks.

Recommended tests:

- Compare atlas choices or target sets in a paired design.
- Compare tractography pipelines in the same subjects.
- Use as a quality-control endpoint unless biological interpretation is
  independently justified.

### 5. Group or Condition Effects

Hypothesis:

- Disease, treatment, genotype, age, or acquisition group changes normalized
  target connectivity from the focal ROI.

Subject-level outcomes:

- Target fractions normalized by total focal-hit tracks.
- Compositional target vectors.
- Entropy, specificity, and `leave_no_target` fractions.

Recommended tests:

- Beta-binomial or logistic mixed models for target-specific fractions.
- Multinomial or Dirichlet-multinomial models for target-composition vectors.
- Permutation tests at the subject level when sample size is modest.
- FDR correction across targets and focal ROIs.

## Current dTOR-985 Demonstration Statistics

The current left-thalamus Harvard-Oxford full run should be reported as a
single-tractogram exploratory utility demonstration, not as a population
inference study.

Recommended descriptive statistics to report:

- Total assigned focal-voxel tracks.
- Total target-hit events.
- Number of nonzero atlas targets.
- `leave_no_target` count and fraction.
- Top target classes by total hits.
- Focal-voxel entropy distribution.
- Per-voxel target class map.
- Sampled real streamlines shown in 2D and 3D anatomical context.

Current run values from the generated full report:

- Assigned focal-voxel tracks: 347,052.
- Target-hit events: 826,890.
- `leave_no_target`: 360.
- Nonzero targets: 116 of 117 atlas target entries.

These values show that FibrTraQ is not just matching one focal ROI to one
endpoint. It is summarizing a many-target connectivity distribution across
individual focal voxels.

## Recommended Implementation Additions

To make publication statistics routine rather than ad hoc, add these helpers in
a future pass:

- `ftq_stat_subject_summary()`: collapse per-voxel FibrTraQ output into
  subject-level target fractions, entropy, specificity, and `leave_no_target`
  metrics.
- `ftq_stat_spatial_gradients()`: fit focal-voxel coordinate-gradient summaries
  per subject.
- `ftq_stat_negative_controls()`: compare correct manifests with wrong-target
  and shifted-ROI manifests.
- `ftq_stat_compare_groups()`: run target-level group comparisons on
  subject-level summaries with FDR correction.
- `ftq_null_shifted_rois()`: generate spatially shifted focal/target masks for
  negative-control manifests.

## Manuscript-Ready Statistical Language

Suggested methods language:

> Validation was performed with known-answer tractography benchmarks and
> matched negative controls. For each manifest, FibrTraQ was run from raw
> tractogram and ROI files. We quantified target-hit rate, false-target rate,
> and the fraction of focal-exiting streamlines that did not intersect any
> target ROI (`leave_no_target`). Correct-target performance was compared with
> wrong-target and shifted-mask controls using paired non-parametric tests across
> bundle-level validation units, with bootstrap confidence intervals for effect
> sizes.

Suggested exploratory language:

> The dTOR-985 left-thalamus analysis was treated as an exploratory utility
> demonstration. We summarized target-hit distributions, focal-voxel target
> specificity, entropy, and `leave_no_target` patterns, and visualized sampled
> real streamlines in anatomical context. Because this demonstration uses a
> single normative tractogram, it supports workflow utility and biological face
> validity but not population-level inference.
