# Preparing the first CRAN submission

Version **0.3.3** is a submission candidate. It has not been submitted to or
accepted by CRAN. The GitHub installation in the README remains the supported
installation route until CRAN publication is verified. The older `RELEASE.md`
records the 0.3.1 research-preview procedure; use this checklist for the first
CRAN submission.

## 1. Review the candidate

- Reconcile `DESCRIPTION`, `CITATION.cff`, `inst/CITATION` and `NEWS.md`.
- Check package-name availability against current and archived CRAN packages
  and current Bioconductor packages, including differences in capitalization.
- Review the distributable archive, licenses, bundled synthetic data, examples,
  links and help pages. Keep private analyses and manuscript files excluded.
- Require the ordinary Linux, macOS and Windows package checks and both Windows
  report-browser checks to pass for the same candidate commit. Review the
  report screenshots as well as their automated interaction results.

## 2. Check the exact submission archive

Once the reviewed candidate and workflow have been authorized for GitHub,
open **Actions → CRAN readiness → Run workflow** and select the candidate
branch. This manual workflow does not submit anything to CRAN.

The [CRAN readiness workflow](.github/workflows/cran-readiness.yaml) builds one
source archive with current R-release, then checks that same archive with
`R CMD check --as-cran` under R-devel on Linux and Windows. The PDF reference
manual is included in those checks. Both check jobs verify the archive's
SHA-256 checksum and source commit before running it. Errors and warnings fail
the job; every NOTE remains available for individual review. Incoming-package,
remote-resource and system-clock checks remain enabled, all suggested packages
are required, and tests run with `NOT_CRAN=false`.

Download and preserve the artifacts before their 30-day retention expires:

- `cran-source-<commit>`: `FibrTraQ_0.3.3.tar.gz`, `source.dcf`, build log,
  source inventory and build-session details.
- `cran-check-<platform>-<commit>`: complete check output, test results,
  reference manual, check summary and check-session details.

Inspect the full logs, including timings and any skipped optional tests.
Resolve errors, warnings and significant notes. Record and explain unavoidable
notes, including a new-submission note, without suppressing them. Local checks
on an older R version are useful additional evidence, not a replacement for
these current-release and R-devel checks. Any source change requires a new
archive and repeat checks.

## 3. Assemble and authorize submission

After the final checks, prepare `cran-comments.md` with the actual check
platforms, R versions, results and explanations for remaining notes. Review
the reference manual and archive inventory. Confirm that the source archive
matches `source.dcf`; submit that archive rather than rebuilding it afterward.

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
before the actual submission. This checklist was prepared on November 12, 2026;
CRAN's review and acceptance remain separate from local or GitHub checks.
