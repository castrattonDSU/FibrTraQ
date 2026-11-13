# Preparing the CRAN resubmission

Version **0.3.5** is a resubmission candidate. CRAN's incoming functional tests
passed for version 0.3.4, but the Windows overall check time was 11 minutes,
exceeding the ten-minute guidance. Version 0.3.5 caches the fixed output-column
dictionary after its first use, removing repeated construction during table
and report generation. The existing scientific tests and report-rendering
settings are retained. The earlier 0.3.3 submission required Windows
report-path and README-link corrections, which remain in the package.

FibrTraQ has not been accepted by CRAN. The GitHub installation in the README
remains the supported route until CRAN publication is verified. The older
`RELEASE.md` records the 0.3.1 research-preview procedure; use this checklist
for the CRAN resubmission.

## 1. Review the candidate

- Reconcile `DESCRIPTION`, `CITATION.cff`, `inst/CITATION` and `NEWS.md`.
- Check package-name availability against current and archived CRAN packages
  and current Bioconductor packages, including differences in capitalization.
- Review the distributable archive, licenses, bundled synthetic data, examples,
  links and help pages. Keep private analyses and manuscript files excluded.
- Require the ordinary Linux, macOS and Windows package checks and both Windows
  report-browser checks to pass for the same candidate commit. Review the
  report screenshots as well as their automated interaction results.
- Require Windows checks with a deeper temporary root comparable to CRAN's.
  Set `TMPDIR`, `TMP` and `TEMP` before starting R and record the effective
  `tempdir()`; retain the ordinary Windows checks as well.
- Check README link destinations against the source archive. Links to files
  intentionally excluded from it must use complete public URLs.

## 2. Check the exact submission archive

Once the reviewed candidate and workflow have been authorized for GitHub,
open **Actions → CRAN readiness → Run workflow** and select the candidate
branch. This manual workflow does not submit anything to CRAN.

The [CRAN readiness workflow](.github/workflows/cran-readiness.yaml) builds one
source archive with current R-release, then checks that same archive with
`R CMD check --as-cran` under R-devel on Linux and Windows. The PDF reference
manual is included in those checks. The check jobs verify the archive's
SHA-256 checksum and source commit before running it. The check must exit
successfully, finish its log, and complete the tests and PDF/HTML manual checks.
An incomplete check fails even if its parsed error and warning counts are zero.

Every completed check must also meet the internal 360-second elapsed-time
limit: Linux, ordinary Windows and Windows with the deeper temporary root.
This six-minute limit provides a margin below CRAN's ten-minute overall-check
guidance; it is a project release gate, not a CRAN rule. Measure the complete
`R CMD check` process, including installation, examples, tests and manuals,
not only the test suite or selected stages. Keep detailed stage timings using
`_R_CHECK_TIMINGS_=0` and the test logs. A timeout or an incomplete run cannot
satisfy the runtime gate. Record the actual time for each environment in the
submission evidence.

The cache retains the original dictionary builder and field definitions.
Regression tests check complete equivalence, construction once per cache and
isolation from caller changes to returned values. Require the full existing
suite and report-rendering coverage, together with those new tests. Review
actual assertion, failure, warning and skip counts in each final check. See
[CRAN runtime guidance](https://contributor.r-project.org/cran-cookbook/docs_issues.html#overall-checktime)
and [repository policy](https://cran.r-project.org/web/packages/policies.html).

The check artifacts include `check-runtime.dcf` with both reported and measured
elapsed times, `check-stage-timings.csv` with the original stage timing strings,
and `test-timings.csv` with per-test timings and assertion counts. Both elapsed
time measurements must meet the 360-second limit. The longer process timeout
is retained so a slow check can finish and preserve its diagnostics; it does
not relax the release gate.

Errors and warnings also fail the job. Every NOTE requires review, and an
unexpected NOTE fails the release gate. Incoming-package, remote-resource,
file-URI and system-clock checks remain enabled, all suggested packages are
required, and tests run with `NOT_CRAN=false`. The additional incoming spelling
check runs on Linux and requires a working Aspell installation and English
dictionary; record those tools and review the flagged terms.

Download and preserve the artifacts before their 30-day retention expires:

- `cran-source-<commit>`: `FibrTraQ_0.3.5.tar.gz`, `source.dcf`, build log,
  source inventory and build-session details.
- `cran-check-<platform>-<temp-profile>-<commit>`: complete check output, test results,
  reference manual, check summary and check-session details.

Inspect the full logs, including timings and any skipped optional tests.
Resolve errors, warnings and significant notes. Record and explain unavoidable
notes, including a new-submission note, without suppressing them. The previous
submission flagged the legitimate terms `ROI` (region of interest),
`Tractograms`/`tractograms` (collections of reconstructed streamlines) and
`tractography` (the reconstruction method); explain any repeated flags in the
resubmission comments. Local checks on an older R version are useful additional
evidence, not a replacement for these current-release and R-devel checks. Any
source change requires a new archive and repeat checks.

## 3. Assemble and authorize submission

After the final checks, prepare `cran-comments.md` with the actual check
platforms, R versions, results and explanations for remaining notes. Explain
the runtime correction and retained coverage, including any intentional CRAN/CI
differences. Retain the history of the Windows path and README-link corrections.
Review the reference manual and archive inventory. Confirm that the source
archive matches `source.dcf`; submit that archive rather than rebuilding it
afterward.

The listed maintainer is **Chase Stratton <cstratton@desu.edu>**. CRAN sends a
confirmation email to that address. Submission and email confirmation require
the maintainer's authorization. Use the [CRAN submission form](https://cran.r-project.org/submit.html)
and avoid a second submission while the first is pending. A GitHub workflow
pass establishes check results, not CRAN acceptance.

## 4. Verify publication and simplify installation

After CRAN acceptance and publication, verify the package page and installation
from CRAN in a fresh R session. Check availability of the Windows and macOS
binaries, which can follow the source publication. Run the bundled example
with the installed CRAN package. Then make this the primary README command:

```r
install.packages("FibrTraQ")
```

Retain GitHub installation as an alternative for development versions, and
record the final published version and source checksum in release evidence.

## Official requirements

Recheck the [CRAN Repository Policy](https://cran.r-project.org/web/packages/policies.html),
[submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
and [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
before the actual submission. This checklist was updated on November 13, 2026;
CRAN's review and acceptance remain separate from local or GitHub checks.
