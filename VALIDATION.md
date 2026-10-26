# FibrTraQ Validation Evidence

This document records the evidence supporting FibrTraQ and the acceptance gates
for the **0.3.1 research preview candidate**. It is intentionally separate from
the [README](README.md): the README is the user guide, while this file describes
what has been verified, what each check establishes, and what remains outside
the scope of software testing.

## Evidence Model

FibrTraQ distinguishes three questions that require different evidence:

1. **Software correctness:** does the implementation reproduce known counts,
   obey its documented semantics, and produce internally consistent outputs?
2. **Computational reproducibility:** do the same versioned inputs and
   parameters produce the same scientific tables across reruns, readers,
   chunk sizes, backends, and supported operating systems?
3. **Biological validity:** are a particular tractogram, anatomical image,
   masks, and coordinate mapping valid representations of the specimen and
   scientific question?

The first two questions can be tested by the package. The third depends on
dataset provenance, image registration, tractography reconstruction, ROI
provenance, and expert anatomical review. A successful FibrTraQ run does not,
by itself, prove those external inputs are biologically valid.

Use the [spatial-verification guide](docs/spatial-verification.md) to check and
record dataset-specific alignment before interpreting a new analysis. An
automatic mapping decision or a computational validation pass cannot replace
independent anatomical evidence.

## Current Evidence

| Evidence | Status | Scope | Artifact |
| --- | --- | --- | --- |
| Bundled known-answer analysis | Complete and user-runnable | Exact end-to-end reproduction of seven scientific tables and 18 independently declared values, plus report and 3-D QC generation | [`inst/extdata/fibrtraq-known-answer`](inst/extdata/fibrtraq-known-answer) |
| Known-answer integrity inventory | Complete and checked before every example run | File size and MD5 verification for every bundled input and expected table | [`checksums.csv`](inst/extdata/fibrtraq-known-answer/checksums.csv) |
| Unit and integration tests | Continuously tested | Counting semantics, coordinate mappings, streamed input, endpoint classes, penetration, validation reconciliation, report generation, and compatibility behavior | [`tests/testthat`](tests/testthat) |
| R and C++ backend parity | Continuously tested on synthetic cases | Exact agreement between the reference and compiled counting implementations across supported analysis modes | [`tests/testthat`](tests/testthat) |
| Cross-platform package checks | Continuously tested | Package installation, tests, examples, and `R CMD check` on Linux, macOS, and Windows | [GitHub Actions](https://github.com/castrattonDSU/FibrTraQ/actions/workflows/R-CMD-check.yaml) |
| Windows report-browser acceptance | Workflow available; inspect the result and screenshots for each candidate | Native Windows Chromium opens relocated synthetic reports, checks local links and downloads, and exercises rendered 3-D controls; this is separate from real-brain anatomical review | [Browser acceptance procedure](.github/scripts/browser-acceptance/README.md) |
| Automatic-reflection regression | Implemented and passed in the preceding development revision | An asymmetric synthetic case makes an incorrect reflection score more focal hits; the engine and file workflow must reject that unconfirmed choice before counting, while explicit mappings preserve known counts | [`test_orientation_guard.R`](tests/testthat/test_orientation_guard.R) |
| dTOR-985 production-scale backend comparison | Reported historical local benchmark for pre-endpoint core outputs | Exact R/C++ agreement was recorded for per-voxel, RSC, and fate tables on 11,820,000 streamlines; performance characterization only, without a self-contained archived evidence bundle | [`backend-benchmark-dtor985.md`](inst/validation/backend-benchmark-dtor985.md) |
| Public benchmark setup assets | Available, not a completed release claim | Reusable manifests and scripts for independent external validation | [`inst/validation`](inst/validation) |

The dTOR-985 benchmark predates the additive endpoint tables. It supports
compiled-backend equivalence and scale for the established per-voxel, RSC, and
fate outputs; it is not evidence that every newer output has been validated on
that dataset.

The last completed pre-candidate check was for version **0.3.0.9002**, commit
`ef8fc0abb7fffee53af0e38b581a66d817d062b1`. Its recorded package checks passed on
Linux, macOS and Windows, with 2,246 passed test expectations on Linux/macOS and
2,240 on Windows, and no test failures, warnings or skips. These are historical
results for that exact revision. They do not establish that a later candidate,
tag, installation, or browser review has passed. The 0.3.1 candidate requires
fresh evidence under the [release checklist](RELEASE.md).

The reflection safeguard is deliberately limited: it blocks an automatically
selected axis reflection when `axis_flip` was omitted. It does not certify an
unreflected coordinate convention or an explicitly supplied transformation.
The regression establishes that the reproduced failure cannot silently proceed
through the tested entry points; it does not establish universal anatomical
accuracy.

## Reproduce The Bundled Analysis

The package-owned fixture is the primary user-facing reproducibility example.
It contains 14 synthetic TrackVis streamlines on a 9 x 9 x 9 grid, focal and
local masks, two targets, a reference image, and an endpoint context mask.

```r
library(FibrTraQ)

example <- ftq_run_example(
  output_dir = file.path(tempdir(), "fibrtraq-known-answer"),
  overwrite = TRUE
)

stopifnot(example$validation$passed)
example$validation$tables
example$validation$claims
browseURL(example$report_index)
```

The run must exactly reproduce:

- the per-focal-voxel result table;
- the RSC summary;
- the mutually exclusive fate summary;
- the ROI-level unique-streamline summary;
- the exact and one-voxel endpoint summary;
- the target-penetration histogram;
- the reference-grid and context-mask endpoint table.

It must also recover all 18 values declared in
[`known-answer-claims.csv`](inst/extdata/fibrtraq-known-answer/known-answer-claims.csv)
from both the versioned expected tables and the newly generated tables. These
claims cover processed and focal-positive streamlines, focal associations,
target intersections, RSC associations, fate categories, endpoint pairs, and
reference-grid endpoint containment.

The fixture is generated by
[`data-raw/generate-known-answer-example.R`](data-raw/generate-known-answer-example.R).
That script constructs the inputs, runs the reference R backend, checks
independently specified counts, and snapshots the seven expected tables.
Regeneration is a maintainer action, not part of normal package use.

## Acceptance Gates

A release candidate is acceptable only when all applicable gates pass:

- the checksum inventory for the bundled fixture is valid;
- `ftq_run_example()` reproduces all seven expected scientific tables exactly;
- all 18 independently declared known answers match both stored and generated
  tables;
- standard, exploratory, and interactive 3-D reports are generated and
  non-empty;
- the complete automated test suite passes;
- `R CMD build` and `R CMD check` pass from the built source archive;
- a fresh installed-package session can locate and run the bundled example;
- reports, local downloads, and 3-D interactions are reviewed in the intended
  browser, including a recorded native Windows review of an extracted bundle;
- the package archive contains the intended example assets and no private
  validation data, credentials, or machine-specific absolute paths;
- all pre-existing output schemas remain compatible unless an intentional
  versioned change is documented.

The GitHub Actions workflow applies the package checks on all three supported
desktop operating-system families. Local release review additionally inspects
the built archive and exercises the installed-package example. Browser
generation or nonempty-file checks are not a substitute for rendering and
interaction. The native Windows browser gate for the 0.3.1 candidate remains
unverified until its evidence is recorded; an unavailable check is not a pass.

## What The Evidence Supports

The current evidence supports these claims:

- FibrTraQ can execute a complete file-based analysis from a fresh R session
  using package-owned inputs.
- The bundled analysis reproduces its seven versioned scientific tables and 18
  independently declared values exactly.
- Existing full-streamline metrics and additive endpoint metrics are reconciled
  by explicit validation rules.
- The reference R and compiled C++ implementations agree on the tested
  analysis modes.
- The package can generate its principal tables, maps, reports, and interactive
  QC artifact without private data or hidden setup code.

## What The Evidence Does Not Support

The current software evidence does not establish that:

- every `.trk` file is in the same physical coordinate space as its MRI;
- a tractography reconstruction is anatomically accurate;
- a supplied focal, target, or anatomical context mask is biologically valid;
- streamline intersections or endpoints are equivalent to axonal connections;
- one synthetic example establishes performance on all species, scanners,
  acquisition protocols, or tractography pipelines;
- a candidate benchmark listed in the validation registry has been completed.

Those claims require versioned external data, acquisition and reconstruction
provenance, geometry audits, appropriate controls, biological review, and
study-specific statistical analysis.

## Extended Validation

The installed [`inst/validation/README.md`](inst/validation/README.md) documents
the extended validation machinery, evidence-bundle policy, acceptance criteria,
and public benchmark setup. The accompanying
[`public-datasets.yml`](inst/validation/public-datasets.yml) identifies whether
each entry is completed evidence, a runnable template, or a registry-only
candidate.

Large tractograms, downloaded atlases, private datasets, and generated reports
are deliberately excluded from the package source. A publication evidence
bundle should archive final inputs or stable source identifiers, checksums,
manifests, software versions, validation records, machine-readable tables,
figures, and reports under a persistent public identifier when licensing and
data-governance rules permit.

## Updating This Record

Add an analysis to the completed-evidence table only after its final artifacts
have been inspected and its reproducibility commands, parameters, provenance,
and limitations are recorded. Setup scripts, manifest templates, partial runs,
and planned analyses must remain labeled as such.
