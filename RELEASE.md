# FibrTraQ release review

This document preserves the **0.3.1 research-preview checklist**. For the
active **0.3.4 CRAN resubmission candidate**, follow [CRAN preparation](CRAN.md).

The checklist below was prepared for R package version **0.3.1**, licensed
under **GPL-2.0-or-later**, with the proposed tag **v0.3.1-rc1**. It is not
evidence that a tag or release exists or that a candidate passed every gate.

The preceding tagged candidate, **v0.3.0-rc2**, contains package version
**0.3.0** under its original MIT license. Preserve it unchanged. Do not
relabel old releases or imply that the current license applies retroactively.

This is a maintainer checklist. For installation and a first analysis, start
with the [README](README.md). Repository visibility and third-party archival
are separate decisions.

## 0. Confirm the public source and history

Review the exact material intended to become public, including tracked files,
commit history, branches, tags, release assets, discussions, and Actions logs.
An ignored or deleted file can remain retrievable from history. Package build
exclusions do not control what GitHub exposes.

If using a clean public source snapshot, preserve the development repository
and its research history privately. Inventory the snapshot's files and hashes,
and verify that its code, tests, fixtures, and documentation match the reviewed
candidate. Build and check the snapshot itself before publication. Do not copy
private validation artifacts or old Git metadata into the public repository.
Creating or changing a remote repository, rewriting history, publishing a tag,
and changing visibility require the maintainer's explicit authorization.

## 1. Check the exact source

Start with a clean checkout of the commit intended for release. Inspect both
tracked and untracked files: `R CMD build .` reads the working directory, not
only the Git index. Local study inputs, manuscript files and private evidence
must remain outside the distributable source archive.

```sh
git status --short
R CMD build .
R CMD check --no-manual --no-build-vignettes FibrTraQ_0.3.1.tar.gz
tar -tzf FibrTraQ_0.3.1.tar.gz
```

Require zero errors and zero warnings. Inspect notes individually. The package
check runs the full test suite; source documentation generation is a separate
development step and should not silently change release sources. Inspect the
archive for required synthetic fixtures, scripts and output guides, and for the
absence of credentials, machine-specific paths, local study audits, manuscripts,
private data and compiled objects.

## 2. Exercise an installed package

Use a temporary library for maintainer checks of the built source archive.
This archive check is separate from the user-facing GitHub installation.

```sh
release_lib=$(mktemp -d /tmp/fibrtraq-release-library.XXXXXX)
R CMD INSTALL -l "$release_lib" FibrTraQ_0.3.1.tar.gz
R_LIBS_USER="$release_lib" Rscript --vanilla -e 'library(FibrTraQ); for (backend in c("cpp", "r")) { x <- ftq_run_example(backend = backend, overwrite = TRUE, verbose = FALSE); stopifnot(x$validation$passed, nrow(x$validation$tables) == 7L, all(x$validation$tables$passed), nrow(x$validation$claims) == 18L, all(x$validation$claims$passed), all(x$validation$reports$passed)) }'
```

Require seven exact scientific-table matches and all 18 declared scientific
claims for both backends. The four report-file checks establish existence and
nonzero size; they do not establish browser interaction or biological validity.
Inspect the figures, connection table, denominators and complete-voxel exports.
Check every local HTML link in the original folder and a copied folder whose
name contains spaces, percent signs and other reserved URL characters.

## 3. Check all supported systems

After the candidate is authorized for a remote check, require the configured
Linux, macOS and Windows `R CMD check` jobs to pass for the exact release
commit. A successful earlier commit does not clear the candidate. If a pull
request is used, recheck the merged commit before attaching a release tag.

Separately, record a native Windows review of a complete extracted results
folder: open the index, Results, companion reports and the 3-D viewer, and
follow download links. Record the operating system, browser/version, source
revision, input fixture, extracted folder name, inspected views and any
failures. Exercise rotation, roll, zoom, reset, and the report's return links.
Repeat navigation from a copied folder whose name contains spaces and a
percent sign. Automated path tests and Windows package checks do not substitute
for this browser review. Automated browser evidence must state the actual
operating system and interactions tested; it must not be labeled a manual
review. An unavailable check remains unverified.

The [Windows report-browser workflow](.github/workflows/report-browser-check.yaml)
provides repeatable browser evidence using the package-owned synthetic example.
Inspect its `windows-report-browser-<commit>` artifact, including
`browser-evidence/browser-check.json` and every saved screenshot. The
[browser acceptance procedure](.github/scripts/browser-acceptance/README.md)
documents its scope and reproduction commands. A green workflow alone does
not establish that its screenshots were visually reviewed.

For a real-brain acceptance run, also record independent spatial evidence and
the explicit mapping using the [spatial-verification guide](docs/spatial-verification.md).
Report rendering and anatomical correctness are separate acceptance questions.

## 4. Align metadata and release evidence

Before tagging, reconcile `DESCRIPTION`, `CITATION.cff`, `inst/CITATION` and
`NEWS.md`: all identify package version **0.3.1**, and the new tag must point to
that checked source. Add the actual release date and release URL only when
they exist. Preserve the analysis commit in scientific outputs and manuscripts;
record later software/report validation separately rather than implying that
all biological analyses were rerun with the candidate.

Keep validation evidence outside the R package: source commit, platform-job
URLs, test and example outcomes, artifact inventory/checksums, report inspection
results and unresolved limitations. Do not include private inputs or credentials
in release notes or source archives. Add a software DOI only after an archival
record exists and its version and license have been verified.

## 5. Create and verify the authorized candidate

After all applicable checks pass and publication is authorized, tag the
approved release commit **v0.3.1-rc1**, provided that name is still unused, and
mark the GitHub release as a **prerelease**. Verify that the tag resolves to the
same source revision recorded in the release evidence. Keep existing access
settings until the public-repository decision is explicitly approved.

Test the actual GitHub installation after the new tag exists. Use the
[README installation](README.md#install-and-try-the-rudimentary-example), replacing
`ref = "main"` with `ref = "v0.3.1-rc1"`. Confirm the installed version is
`0.3.1` and its `License` field is `GPL (>= 2)`.

Restart R, run the README's bundled example and verify the recorded `RemoteSha`
against the tagged commit. Private repository installations require an authorized
GitHub account; reuse the workstation's credential store and never put access
tokens in scripts or reports.

Record this post-tag installation result before announcing the candidate to
users. A tag name written in documentation is not proof that installation from
that tag works.
