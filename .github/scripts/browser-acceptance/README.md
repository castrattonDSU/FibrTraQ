# Windows report browser acceptance

The `Windows report browser check` workflow builds and installs the checked-out
candidate on a native Windows runner. It uses only the package-owned synthetic
example in a two-format matrix (`trk`, `tck`), copies each complete result
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

Use a short output root. The browser check records the longest file path and
requires fewer than 260 characters, retaining spaces, `#`, `%` and Unicode in
the copied folder. Some Windows applications still use the traditional
[260-character path limit](https://learn.microsoft.com/en-us/windows/win32/fileio/maximum-file-path-limitation).
Creating a file successfully in R does not establish that a browser can open
the same long path. This check does not claim support for arbitrary path lengths.

CSV links must produce an actual download whose bytes match the linked file.
For a PNG, Chromium may either download the file or display it directly. The
check accepts a byte-matching download, or an observed navigation to the exact
linked URL with `image/png` content type and decoded dimensions matching the
PNG header. It records which behavior occurred and captures an opened image;
file existence alone cannot satisfy this check.

Each format's workflow artifact contains the synthetic report folders, input manifest,
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
```

The optional second `generate-fixture.R` argument selects `trk` or `tck`;
omitting it retains the `trk` default. Each output directory must be new.
The browser script refuses non-Windows hosts
so a different operating system cannot be mistaken for Windows evidence.
