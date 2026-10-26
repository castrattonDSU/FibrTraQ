# dTOR-985 Backend Benchmark

This record summarizes the local production-size backend benchmark run on
2026-06-05. It compares the reference R backend with the compiled C++ backend
on the same dTOR-985 Harvard-Oxford left-thalamus workflows.

## Configuration

- Dataset: dTOR-985 normative connectome validation data.
- Tractogram: `raw/dTOR_full_tractogram.trk`.
- Streamlines processed per run: 11,820,000.
- Focal ROI: Harvard-Oxford left thalamus.
- Focal voxels: 1,321.
- Workflows:
  - `pilot`: 12 Harvard-Oxford targets.
  - `full`: full Harvard-Oxford cortical/subcortical target set.
- Maps and HTML reports were disabled to isolate analysis/counting time.
- Core outputs compared:
  - per-focal-voxel hits CSV,
  - RSC summary CSV,
  - fate summary CSV.

## Results

| Workflow | R elapsed seconds | C++ elapsed seconds | Elapsed speedup | R count-stage seconds | C++ count-stage seconds | Count-stage speedup | Core outputs identical |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Pilot | 3305.973 | 469.772 | 7.04x | 3305.182 | 469.080 | 7.05x | yes |
| Full | 2820.164 | 495.199 | 5.70x | 2814.884 | 490.641 | 5.74x | yes |

Both workflows produced identical per-voxel, RSC, and fate outputs between the
R and C++ backends. The C++ backend is therefore used by `backend = "auto"` when
compiled bindings are available, while `backend = "r"` remains available as the
explicit reference implementation.

## Output Files

The local benchmark artifacts were written under the dTOR validation data root:

```text
public-validation/dtor985/outputs/harvard-oxford-left-thalamus/backend-benchmarks/20260605-153301
```

Key files:

```text
dtor985_backend_timing.csv
dtor985_backend_profiles.csv
dtor985_backend_equivalence.csv
```

The benchmark can be regenerated on the validation workstation with:

```r
source("inst/scripts/fibrtraq-dtor985-backend-benchmark.R")
```
