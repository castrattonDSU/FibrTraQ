# Native TCK validation

This evidence concerns the 0.3.2.9000 development implementation. Native MRtrix
input uses the existing analysis and report workflow. Software checks do not
establish anatomical accuracy for arbitrary input datasets.

## Synthetic correctness and regression coverage

Reader tests cover Float32/Float64 in both byte orders, padded MRtrix headers,
chunk boundaries, optional counts, malformed records, truncation and EOF.
Orientation tests use independently constructed asymmetric fixtures and known
counts: translated, rotated, reflected, oblique and sheared grids; millimetres,
metres and micrometres; voxel boundaries; all three target-assignment modes;
both voxelization modes; and point-order reversal. Equivalent synthetic TCK
and TRK inputs reproduce the same scientific results in R and C++. Complete
report tests compare tables, maps, density context and 3-D scene data.

TCK points retain their defined RAS+ world-millimetre convention. Declared
NIfTI units and valid, consistent spatial transforms are required; orientation
is not inferred by maximizing mask overlap. Existing TRK mapping remains
available. The Windows browser workflow has separate TRK and TCK fixtures;
check its result for the exact revision being released.

## Public human demonstration

The optional validation uses the native 10,000-streamline MRtrix tractogram
published with the [archived Neurodesk tutorial](https://doi.org/10.5281/zenodo.23201854).
Its public source is BTC_preop control participant `sub-CON02`. Exact processed
files, revision, checksums, licenses and the source-version caveat are in
[tck-public-sources.json](tck-public-sources.json).

The TCK contains 272,137 points. NiBabel supplies an independent coordinate and
point-count oracle. Both FibrTraQ readers reproduce every decoded coordinate
and point count byte for byte under the declared canonical encodings. Native
analysis tables and maps also match exactly between streaming/C++ and memory/R,
and after reversing stored point order. The final native execution/invariance
summary passes 58 of 58 checks, including 23 exact backend/memory comparisons
and 23 exact reversal comparisons.

MRtrix converts the coregistered tissue images to NIfTI. Preparation removes
only a trailing singleton dimension from cortical GM, preserving voxel values
and spatial transforms. The focal mask is GM/WM interface above 0.5 at world
x < 0. Targets are cortical GM above 0.5 at x >= 0 and x < 0, with focal voxels
excluded from the latter. These are geometric tissue partitions, not atlas
parcels or validated hemisphere labels. Report backgrounds show cortical-GM
probability, not T1 intensity. No registration is inferred.

Full-streamline assignment, all focal voxels, segment interpolation and round
discretization identify 4,958 focal-intersecting streamlines in 76,215 focal
voxels, producing 21,345 voxel–streamline associations. Positive-x and
negative-x targets cover 1,225 and 6,913 focal voxels: 1.6073% and 9.0704% RSC.
These values describe the downsampled demonstration, not a complete human
connectome. Computational repeats remain one biological subject.

## Conversion precision is a separate result

Strict equality with a conventional NiBabel-written Float32 TRK conversion
does **not** hold. The NiBabel round trip changes coordinates by at most
0.000023 mm. Three streamline voxel paths change at discretization boundaries;
one affects focal results. Streamline 9982 gains voxel (49, 142, 93), changing
the focal count from 4,958 to 4,959 and negative-x coverage from 6,913 to 6,914
voxels. The RSC difference is 100 / 76,215 = 0.00131208 percentage points.
Voxel coordinates and streamline indices here are one-based.

The strict comparison retains failed checks and a nonzero exit status.
Separate summaries distinguish native execution/invariance from cross-format
equality. A precision audit reconstructs the changed paths and checks their
predicted effects against complete tables and density maps without relaxing
count tolerances. Successful causal accounting does not make the original
comparison an exact match. The final causal audit passes all 72 checks while
retaining `exact_format_parity = false`; the strict representation comparison
records 19 passed and 20 failed checks. The failed checks share the documented
voxel-boundary effects across their dependent maps, tables and report exports.

Prefer the original native file. Record conversion, registration, rounding
choices and ROI definitions when comparing results. Matching tables establish
computational consistency, not biological validity.

See the [validation instructions](README.md) for acquisition and analysis.
Archive provenance, strict comparisons and precision-audit records together
for manuscript evidence. Imaging files are not bundled with the package.
