# Windows report browser acceptance

The `Windows report browser check` workflow builds and installs the checked-out
candidate on a native Windows runner. It uses only the package-owned synthetic
example in a four-case matrix (`trk`, `tck` × `normal`, `near-limit`), copies each complete result
bundle into a folder containing spaces, `#`, `%` and Unicode, then opens the
copied HTML directly with `file://` URLs in Chromium.
The original outputs remain intact. No web server or external data are needed.

The TCK case encodes the known synthetic points as native world-millimeter
coordinates and uses copies of the fixture images with declared NIfTI geometry.
Only format-specific provenance and copied-image checksums differ in its
expected tables; the bundled inputs and scientific values remain unchanged.

Each format verifies the same seven scientific tables and 18 declared values before
testing local links, fragment targets, loaded report images, navigation from
the index, CSV downloads, PNG opening, and actual canvas rendering during drag, edge
roll, zoom, reset, and a layer toggle. It checks advanced path/endpoint controls
and captures the main views and a narrow viewer. Pixel checks inspect the centre
of the scene so a moving compass alone cannot satisfy an interaction check.
Playwright is pinned here as a development-only dependency; it is not installed
with the R package.

The image check opens report sections containing figures and scrolls each image
into view before requiring successful decoding. This covers Chromium's deferred
loading of collapsed content while still rejecting missing or corrupt images.
The evidence records initial and decoded image states and an expanded Results
screenshot so the advanced figures can also be inspected.

The normal cases retain the existing browser assertions. Each near-limit case
first creates a short-root baseline, then runs the complete analysis again in a
nested root sized so that its longest generated file path is exactly 255 UTF-16
code units. It checks the same scientific known answers and file inventory, then
copies that deep-generated bundle to a reserved-character root whose longest
file path is exactly 259 units. Every browser assertion runs on the relocated
bundle, including image loading, downloads, report links, and viewer controls.

The near-limit cases also request an output root whose planned filenames exceed
the supported limit. The public workflow must return the specific path-length
error before creating its output directory, and existing baseline files and a
sentinel must remain intact. `over-limit-evidence.json` records the observed
outcome before the assertions, including any unexpected failure.

Use a short acceptance output root so the harness can add the measured padding.
The browser evidence records both absolute roots, nesting depth, longest
relative filenames, planned and actual maximum lengths, and the rejected
request. JavaScript string lengths and the R fixture use UTF-16 code units;
the accented character in these paths occupies one unit. Some Windows
applications still use the traditional
[260-character path limit](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation).
Creating a file successfully in R does not establish that a browser can open
the same long path. These cases test direct generation under roots at least as
long as the temporary prefix seen in the failed CRAN check and browser opening
near the supported boundary. They do not claim support for arbitrary CRAN
directory layouts or unlimited Windows path lengths.

CSV links must produce an actual download whose bytes match the linked file.
For a PNG, Chromium may either download the file or display it directly. The
check accepts a byte-matching download, or an observed navigation to the exact
linked URL with `image/png` content type and decoded dimensions matching the
PNG header. It records which behavior occurred and captures an opened image;
file existence alone cannot satisfy this check.

Each format and path mode's workflow artifact contains the synthetic report folders, input manifest,
screenshots, and `browser-evidence/browser-check.json`. Review the screenshots
before accepting visual layout. A successful automated check establishes the
listed browser behaviors on that recorded Windows/browser/package revision;
it is not manual user acceptance, a real-brain Windows analysis, or anatomical
validation. Keep those evidence items separate in the release record.

For an authorized native Windows reproduction, install the intended R package
in a fresh session, install Node.js 24, then from this directory run:

```powershell
npm ci --ignore-scripts
npx playwright install chromium
Rscript --vanilla generate-fixture.R "C:/fibrtraq-browser-review-trk" trk
node report-browser-check.cjs "C:/fibrtraq-browser-review-trk/browser-fixture.json"
Rscript --vanilla generate-fixture.R "C:/fibrtraq-browser-review-tck" tck
node report-browser-check.cjs "C:/fibrtraq-browser-review-tck/browser-fixture.json"
Rscript --vanilla generate-fixture.R "C:/ftq-trk-long" trk near-limit
node report-browser-check.cjs "C:/ftq-trk-long/browser-fixture.json"
Rscript --vanilla generate-fixture.R "C:/ftq-tck-long" tck near-limit
node report-browser-check.cjs "C:/ftq-tck-long/browser-fixture.json"
```

The optional second `generate-fixture.R` argument selects `trk` or `tck`;
omitting it retains the `trk` default. The optional third argument selects
`normal` (the default) or `near-limit`. Each output directory must be new.
The browser script refuses non-Windows hosts
so a different operating system cannot be mistaken for Windows evidence.
