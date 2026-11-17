"use strict";
// Windows is required for acceptance. An explicit --local-smoke run exercises
// the same browser interactions for development on the local operating system.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const os = require("node:os");
const crypto = require("node:crypto");
const {pathToFileURL, fileURLToPath} = require("node:url");
const {chromium} = require("playwright");

assert(process.argv[2], "Supply browser-fixture.json from generate-fixture.R");
const localSmoke = process.argv[3] === "--local-smoke";
assert(process.argv.length === 3 || (process.argv.length === 4 && localSmoke));
const manifestPath = path.resolve(process.argv[2]);
const manifest = JSON.parse(fs.readFileSync(manifestPath, "utf8"));
const root = fs.realpathSync(manifest.relocated_root);
const evidenceDir = path.join(path.dirname(manifestPath), "browser-evidence");
fs.mkdirSync(evidenceDir, {recursive: true});
const result = {
  passed: false,
  platform: process.platform,
  os_release: os.release(),
  package_version: manifest.package_version,
  source_commit: manifest.source_commit,
  fixture: manifest.fixture,
  path_mode: manifest.path_mode,
  checks: [], pages: [], downloads: [], screenshots: [], page_errors: [], image_loading: [],
  scope: localSmoke ? "Local browser smoke test" : "Native Windows browser acceptance",
  native_windows_acceptance: !localSmoke
};
let browser, context, page, viewer;
const hash = bytes => crypto.createHash("sha256").update(bytes).digest("hex");
function check(value, description, details) {
  result.checks.push({description, passed: Boolean(value), ...(details === undefined ? {} : {details})});
  assert(value, description);
}
function localFile(url, bundleRoot) {
  const parsed = new URL(url);
  assert.equal(parsed.protocol, "file:", "Expected a local file URL");
  const file = fileURLToPath(parsed);
  const real = fs.realpathSync(file);
  const relative = path.relative(bundleRoot, real);
  assert(relative && !relative.startsWith(`..${path.sep}`) && relative !== ".." && !path.isAbsolute(relative),
    `Local link must remain inside the relocated bundle: ${url}`);
  assert(fs.statSync(real).isFile(), `Local link must identify a file: ${url}`);
  return real;
}
const singleURL = pathToFileURL(path.join(root, manifest.report)).href;
async function settled() {
  await page.evaluate(async () => {
    if (document.fonts) await document.fonts.ready;
    await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)));
  });
}
async function open(url) {
  await page.goto(url, {waitUntil: "load"});
  await settled();
}
async function screenshot(name, fullPage = false) {
  const filename = `${name}.png`;
  await page.screenshot({path: path.join(evidenceDir, filename), fullPage});
  result.screenshots.push(filename);
}
async function exposeReportImages(relative) {
  const initiallyUnloaded = await page.locator("img").evaluateAll(images => images
    .filter(image => !image.complete || image.naturalWidth === 0)
    .map(image => ({url: image.src, complete: image.complete,
      inside_closed_details: Boolean(image.closest("details:not([open])"))})));
  // Chromium can defer images inside collapsed details. Exercise the same
  // disclosure and scrolling actions a reader uses before requiring decoding.
  // The assertions below still reject missing, corrupt or incorrectly linked images.
  const disclosures = page.locator("details");
  let expanded = 0;
  for (let i = 0; i < await disclosures.count(); i++) {
    const section = disclosures.nth(i);
    if (await section.locator("img").count() > 0 &&
        !await section.evaluate(element => element.open)) {
      await section.locator(":scope > summary").click();
      expanded++;
    }
  }
  const images = page.locator("img");
  for (let i = 0; i < await images.count(); i++) {
    await images.nth(i).scrollIntoViewIfNeeded();
  }
  await page.waitForFunction(() => [...document.images].every(image => image.complete));
  const decoded = await images.evaluateAll(elements => elements.map(image => ({
    url: image.src, complete: image.complete,
    natural_width: image.naturalWidth, natural_height: image.naturalHeight
  })));
  result.image_loading.push({relative_path: relative, initially_unloaded: initiallyUnloaded,
    expanded_sections: expanded, decoded});
  await page.evaluate(() => window.scrollTo(0, 0));
  await settled();
}
async function observeFileClick(link, targetURL, allowImageNavigation) {
  const startingURL = page.url();
  let finish;
  // Resolve every observation, including timeouts, so a failed click cannot
  // leave a rejecting event-wait promise behind. Remove both listeners once
  // either supported browser behavior has been observed.
  const observation = new Promise(resolve => {
    let completed = false, timer;
    const onDownload = download => finish({behavior: "download", download});
    const onNavigation = frame => {
      if (frame !== page.mainFrame() || frame.url() === startingURL) return;
      if (allowImageNavigation && frame.url() === targetURL) {
        finish({behavior: "image_navigation"});
      } else {
        finish({error: new Error(`Unexpected navigation after file click: ${frame.url()}`)});
      }
    };
    finish = outcome => {
      if (completed) return;
      completed = true;
      clearTimeout(timer);
      page.off("download", onDownload);
      page.off("framenavigated", onNavigation);
      resolve(outcome);
    };
    page.on("download", onDownload);
    page.on("framenavigated", onNavigation);
    timer = setTimeout(() => finish({error: new Error(
      `No supported file-open behavior observed within 15000ms: ${targetURL}`
    )}), 15000);
  });
  try {
    await link.click();
  } catch (error) {
    finish({error});
  }
  const outcome = await observation;
  if (outcome.error) throw outcome.error;
  return outcome;
}

async function inspectResultsFirst(bundleRoot, run, isKnownExample) {
  const viewport = page.viewportSize();
  await page.setViewportSize({width: 1280, height: 800});
  await page.evaluate(() => window.scrollTo(0, 0));
  await settled();
  const results = page.locator("#ftq-primary");
  const table = results.getByRole("table").filter({
    has: page.getByRole("columnheader", {name: "Destination", exact: true})
  });
  const label = `${run.subject} / ${run.focal_region}`;
  const sharing = await page.getByRole("contentinfo").innerText();
  check(/entire output folder/i.test(sharing) && sharing.includes("fibrtraq-run-index.csv") &&
    /extract/i.test(sharing), `The report identifies the complete folder to zip and extract: ${label}`);
  check(await table.count() === 1, `One RSC table for ${label}`);
  const bounds = await table.boundingBox();
  check(bounds && bounds.y >= 0 && bounds.y + bounds.height <= 800,
    `All synthetic RSC rows are visible on a laptop's first screen: ${label}`, bounds);
  await screenshot(`first-screen-${run.subject}-${run.focal_region}`);
  const heatmaps = results.getByRole("heading", {name: /density heatmaps/i});
  const heatmapBounds = await heatmaps.boundingBox();
  check(heatmapBounds && heatmapBounds.y >= bounds.y + bounds.height,
    `Density heatmaps follow the RSC table: ${label}`);
  const local = table.getByRole("row").filter({has: page.getByRole("cell", {name: "local", exact: true})});
  check(await local.count() === 1 && await local.locator("a").count() === 0,
    `Local classification is distinct from the target map downloads: ${label}`);
  check(await table.locator("a").count() === run.targets.length,
    `Every named destination has one map link: ${label}`);
  for (const target of run.targets) {
    const row = table.getByRole("row").filter({
      has: page.getByRole("link", {name: target, exact: true})
    });
    const link = row.getByRole("link", {name: target, exact: true});
    check(await link.count() === 1, `Reader can locate the ${target} map from its RSC row: ${label}`);
    if (isKnownExample) {
      const values = (await row.getByRole("cell").allTextContents()).map(text => text.trim());
      const uniqueCount = {cortex: 6, brainstem: 4}[target];
      check(+values[1] === uniqueCount && +values[2] === 2 && +values[3] === 100,
        `Known example row shows unique streamlines, focal coverage and RSC for ${target}`, values);
    }
    const url = await link.evaluate(element => element.href);
    const file = localFile(url, bundleRoot);
    check(path.dirname(file) === path.join(bundleRoot, run.subject, run.focal_region, "maps") &&
      path.basename(file).endsWith(`_per-voxel-hits_target-${target}.nii.gz`),
      `RSC destination links to its own run's count map: ${label} / ${target}`);
    const {download} = await observeFileClick(link, url, false);
    const checksum = hash(fs.readFileSync(file));
    check(await download.failure() === null && hash(fs.readFileSync(await download.path())) === checksum,
      `Downloaded destination NIfTI matches the relocated file: ${label} / ${target}`);
    result.downloads.push({extension: "nii.gz", behavior: "download", target: target, run: label,
      filename: download.suggestedFilename(), target_url: url, sha256: checksum});
  }
  await page.setViewportSize(viewport);
}
async function canvasPixels() {
  return viewer.locator("#scene").evaluate(canvas => {
    // Repeated reads on the live canvas can make Chromium switch rasterizers
    // and change antialiasing even when the scene is identical. Read a copy.
    const snapshot = document.createElement("canvas");
    snapshot.width = canvas.width;
    snapshot.height = canvas.height;
    const context = snapshot.getContext("2d", {willReadFrequently: true});
    context.drawImage(canvas, 0, 0);
    const pixels = context.getImageData(0, 0, snapshot.width, snapshot.height).data;
    let centralPainted = 0, centralFingerprint = 2166136261;
    const centralColors = new Set();
    for (let y = Math.floor(canvas.height * .2); y < canvas.height * .8; y++) {
      for (let x = Math.floor(canvas.width * .2); x < canvas.width * .8; x++) {
        const i = (y * canvas.width + x) * 4;
        for (let channel = 0; channel < 4; channel++) {
          centralFingerprint = Math.imul(centralFingerprint ^ pixels[i + channel], 16777619) >>> 0;
        }
        if (pixels[i + 3] > 0) {
          centralPainted++;
          centralColors.add(`${pixels[i]},${pixels[i + 1]},${pixels[i + 2]},${pixels[i + 3]}`);
        }
      }
    }
    return {width: canvas.width, height: canvas.height, centralPainted,
      centralColors: centralColors.size, centralFingerprint, image: snapshot.toDataURL("image/png")};
  });
}
async function projectedSceneBounds() {
  // Inspect the coordinates actually projected by the generated viewer,
  // independently of whether its canvas happens to contain some painted pixels.
  return viewer.evaluate(() => {
    const points = [...data.brain, ...data.focal];
    for (const target of data.targets) {
      points.push(target.point, ...(target.boundary || []));
    }
    for (const streamline of data.streamlines) {
      points.push(...streamline.points);
      for (const endpoint of streamline.endpoints || []) points.push(endpoint.point);
    }
    const projected = points.map(point => project(point));
    const rectangle = canvas.getBoundingClientRect();
    return {point_count: projected.length, canvas_width: rectangle.width, canvas_height: rectangle.height,
      finite: projected.every(point => Number.isFinite(point.x) && Number.isFinite(point.y)),
      min_x: Math.min(...projected.map(point => point.x)), max_x: Math.max(...projected.map(point => point.x)),
      min_y: Math.min(...projected.map(point => point.y)), max_y: Math.max(...projected.map(point => point.y))};
  });
}

function bundleFilePaths(directory) {
  return fs.readdirSync(directory, {withFileTypes: true}).flatMap(entry => {
    const filename = path.join(directory, entry.name);
    return entry.isDirectory() ? bundleFilePaths(filename) : [filename];
  });
}

(async () => {
  try {
    if (!localSmoke) {
      check(process.platform === "win32", "Acceptance browser runs on native Windows");
      check(manifest.os === "Windows", "Synthetic reports were generated on Windows");
    }
    check(manifest.copy_checksums_identical && manifest.original_outputs_preserved,
      "Relocated bundle preserves the generated original and identical file bytes");
    check(manifest.scientific_tables_passed === 7 && manifest.declared_values_passed === 18,
      "Installed-package fixture passed all scientific known answers");
    const bundlePaths = bundleFilePaths(root);
    const longestPath = bundlePaths.reduce((a, b) => a.length >= b.length ? a : b, "");
    const originalRoot = fs.realpathSync(manifest.original_root);
    const originalPaths = bundleFilePaths(originalRoot);
    const longestOriginalPath = originalPaths.reduce((a, b) => a.length >= b.length ? a : b, "");
    const pathDetails = (directory, longest) => ({
      absolute_root: directory, root_utf16_units: directory.length,
      directory_components: directory.split(path.sep).filter(Boolean).length,
      maximum_utf16_units: longest.length,
      longest_relative_path: path.relative(directory, longest)
    });
    result.path_lengths = {maximum_characters: longestPath.length,
      longest_relative_path: path.relative(root, longestPath),
      original: pathDetails(originalRoot, longestOriginalPath),
      relocated: pathDetails(root, longestPath)};
    check(longestPath.length < 260,
      "Windows browser fixture uses file paths below the 260-character compatibility limit; choose a shorter output root if needed",
      result.path_lengths);
    if (manifest.path_mode === "near-limit") {
      const boundary = manifest.path_boundary;
      result.path_boundary = boundary;
      check(boundary.generated_directly_at_boundary,
        "Near-limit reports were generated at their deep original root");
      check(longestOriginalPath.length === 255 && boundary.original_expected_maximum_utf16_units === 255,
        "Longest directly generated file path is 255 UTF-16 code units", result.path_lengths.original);
      check(longestPath.length === 259 && boundary.relocated_expected_maximum_utf16_units === 259,
        "Longest relocated file path is 259 UTF-16 code units", result.path_lengths.relocated);
      check(boundary.original_planned_maximum_utf16_units <= 259,
        "Successful generation also satisfies the package's planned-path preflight",
        {planned_maximum_utf16_units: boundary.original_planned_maximum_utf16_units});
      // CRAN's failed check used a 44-character temporary prefix. This records
      // comparable absolute-root depth without claiming every CRAN root fits.
      check(originalRoot.length >= 44,
        "Direct generation uses an absolute root at least as long as the observed CRAN temporary prefix");
      const inventory = (directory, paths) => paths.map(filename => path.relative(directory, filename)).sort();
      const baselineRoot = fs.realpathSync(boundary.baseline_root);
      check(boundary.relative_inventory_identical &&
        JSON.stringify(inventory(originalRoot, originalPaths)) ===
          JSON.stringify(inventory(baselineRoot, bundleFilePaths(baselineRoot))),
      "Deep generation preserves the short-root file inventory");
      const rejected = boundary.rejected_request;
      check(rejected.planned_maximum_utf16_units > 259 &&
        rejected.message.startsWith("Windows output path exceeds the supported limit of 259 UTF-16 code units"),
      "Over-limit public workflow rejects the planned path with the specific actionable error", rejected);
      check(rejected.output_root_absent && !fs.existsSync(rejected.output_root),
        "Rejected workflow leaves no partial output directory");
      const sentinelMD5 = crypto.createHash("md5")
        .update(fs.readFileSync(rejected.existing_output_sentinel)).digest("hex");
      check(rejected.existing_output_sentinel_preserved &&
        sentinelMD5 === rejected.existing_output_sentinel_md5 && rejected.baseline_outputs_preserved,
      "Rejected workflow preserves pre-existing outputs and sentinel bytes");
    }
    browser = await chromium.launch({headless: true,
      executablePath: localSmoke ? process.env.FTQ_LOCAL_CHROMIUM : undefined});
    result.browser_version = browser.version();
    result.playwright_version = require("playwright/package.json").version;
    context = await browser.newContext({viewport: {width: 1440, height: 1000},
      deviceScaleFactor: 1, acceptDownloads: true});
    context.setDefaultTimeout(15000);
    page = await context.newPage();
    page.on("pageerror", error => result.page_errors.push(error.message));
    // Reports are self-contained. Do not contact external sites during review.
    await context.route(/^https?:/, route => route.abort());

    // Check what readers can actually see after choosing a run. Assets must
    // resolve within this relocated bundle, not the still-existing original.
    async function inspectRun(bundleRoot, label, expectedRun) {
      await exposeReportImages(label);
      const content = await page.evaluate(() => ({
        title: document.title,
        identifiers: [...document.querySelectorAll("[id], a[name]")].flatMap(e => [e.id, e.name].filter(Boolean)),
        links: [...document.querySelectorAll("a[href]")].map(e => ({raw: e.getAttribute("href"), url: e.href})),
        images: [...document.images].map(e => ({url: e.src, loaded: e.complete && e.naturalWidth > 0 && e.naturalHeight > 0}))
      }));
      check(Boolean(content.title.trim()), `HTML title is present: ${label}`);
      check(content.images.length > 0 && content.images.every(i => i.loaded), `All images load: ${label}`);
      const shown = [];
      for (const img of content.images) {
        if (new URL(img.url).protocol === "file:") shown.push(localFile(img.url, bundleRoot));
      }
      if (expectedRun) {
        const dir = path.join(bundleRoot, expectedRun.subject, expectedRun.focal_region, "figures");
        check(shown.length > 0 && shown.every(name => path.dirname(name) === dir),
          `Images correspond to selected subject and region: ${label}`);
      }
      const linkFiles = [];
      for (const link of content.links) {
        check(!link.raw.includes("\\"), `HTML link uses URL separators: ${label}`);
        const target = new URL(link.url);
        if (target.protocol !== "file:") continue;
        const dest = localFile(link.url, bundleRoot);
        linkFiles.push(dest);
        if (target.hash && /[.]html$/i.test(target.pathname)) {
          check(content.identifiers.includes(decodeURIComponent(target.hash.slice(1))),
            `Local HTML fragment resolves: ${label} ${target.hash}`);
        }
      }
      if (expectedRun) {
        const here = path.join(bundleRoot, expectedRun.subject, expectedRun.focal_region, "tables");
        check(linkFiles.some(name => path.dirname(name) === here && name.endsWith("_per-voxel-hits.csv")),
          `Scientific CSV available for selected run: ${label}`);
      }
      check(linkFiles.length > 0, `Report includes working local links: ${label}`);
      result.local_references_checked = (result.local_references_checked || 0) + linkFiles.length;
      result.pages.push({view: label, image_count: shown.length, local_link_count: linkFiles.length});
    }
    async function openViewer() {
      await page.locator("#ftq-open-viewer").click();
      const iframe = await page.locator("#ftq-viewer").elementHandle();
      viewer = await iframe.contentFrame();
      assert(viewer, "Viewer frame is attached");
      await viewer.locator("#scene").waitFor({state:"visible"});
      await viewer.locator("#scene").scrollIntoViewIfNeeded();
      await settled();
    }
    check(manifest.html_files.length === 1 && manifest.html_files[0] === manifest.report,
      "Single run writes one HTML document in its focal-region folder");
    await open(singleURL);
    check(await page.locator("#ftq-viewer").count() === 0, "Viewer initializes on request");
    await screenshot("results", false);
    await inspectResultsFirst(root, manifest.single_run, true);
    await inspectRun(root, "single", manifest.single_run);
    check(await page.locator("#ftq-primary img").count() > 0 &&
      await page.locator("#ftq-targets img").count() === 5,
    "Results and all five target figures are accessible on the same page");
    await screenshot("results-expanded", true);
    for (const extension of ["csv", "png"]) {
      const link = page.locator(`#ftq-targets a[download][href$=".${extension}"]`).first();
      check(await link.count() === 1, `Target report offers a ${extension.toUpperCase()} download`);
      const url = await link.evaluate(element => element.href);
      const expected = fs.readFileSync(localFile(url, root));
      const observed = await observeFileClick(link, url, extension === "png");
      if (observed.behavior === "download") {
        const download = observed.download;
        check(await download.failure() === null, `${extension.toUpperCase()} download completes`);
        const received = fs.readFileSync(await download.path());
        check(hash(received) === hash(expected), `${extension.toUpperCase()} download matches the linked file`);
        result.downloads.push({extension, behavior: "download", target_url: url,
          page_url: page.url(), filename: download.suggestedFilename(),
          target_sha256: hash(expected), sha256: hash(received)});
      } else {
        // Chromium can display a local PNG despite an HTML download attribute.
        // Accept only the observed exact-target navigation and a decoded image
        // matching the dimensions declared by that linked PNG's IHDR header.
        check(extension === "png", "Only PNG links may open as browser images");
        check(expected.length >= 24 && expected.subarray(0, 8).equals(Buffer.from([137, 80, 78, 71, 13, 10, 26, 10])) &&
          expected.toString("ascii", 12, 16) === "IHDR", "Linked image has a valid PNG signature and dimension header");
        const expectedWidth = expected.readUInt32BE(16), expectedHeight = expected.readUInt32BE(20);
        await page.waitForLoadState("load");
        await page.waitForFunction(() => document.images.length === 1 && document.images[0].complete &&
          document.images[0].naturalWidth > 0 && document.images[0].naturalHeight > 0);
        const image = await page.evaluate(() => ({url: location.href, content_type: document.contentType,
          natural_width: document.images[0].naturalWidth, natural_height: document.images[0].naturalHeight}));
        check(image.url === url && image.content_type === "image/png", "PNG opens at the exact linked URL with image/png content type");
        check(image.natural_width === expectedWidth && image.natural_height === expectedHeight,
          "Browser-decoded PNG dimensions match the linked file header");
        await screenshot("png-opened-in-browser", true);
        result.downloads.push({extension, behavior: "image_navigation", target_url: url,
          target_sha256: hash(expected), page_url: image.url, content_type: image.content_type,
          natural_width: image.natural_width, natural_height: image.natural_height});
      }
    }

    await open(singleURL);
    await openViewer();
    const initial = await canvasPixels();
    check(initial.width > 0 && initial.height > 0 && initial.centralPainted > 20 && initial.centralColors > 8,
      "3-D canvas paints nonuniform scene content away from the compass",
      {central_painted_pixels: initial.centralPainted, central_colors: initial.centralColors});
    await screenshot("viewer-initial");
    const bounds = await viewer.locator("#scene").boundingBox();
    const x = bounds.x + bounds.width * .55, y = bounds.y + bounds.height * .45;
    await page.mouse.move(x, y);
    await page.mouse.down();
    await page.mouse.move(x + 120, y + 70, {steps: 12});
    await page.mouse.up();
    await settled();
    const rotated = await canvasPixels();
    check(rotated.centralFingerprint !== initial.centralFingerprint,
      "Dragging changes central scene pixels, not only the compass");
    await screenshot("viewer-rotated");
    await page.mouse.wheel(0, -300);
    await settled();
    const zoomed = await canvasPixels();
    check(zoomed.centralFingerprint !== rotated.centralFingerprint,
      "Mouse wheel changes rendered zoom in the scene centre");
    await screenshot("viewer-zoomed");
    result.canvas = {width: initial.width, height: initial.height,
      central_painted_pixels: initial.centralPainted, central_colors: initial.centralColors,
      initial_fingerprint: initial.centralFingerprint, rotated_fingerprint: rotated.centralFingerprint,
      zoomed_fingerprint: zoomed.centralFingerprint};
    await viewer.locator("#reset").click();
    await settled();
    check((await canvasPixels()).image === initial.image, "Reset restores the initial rendered scene exactly");
    await screenshot("viewer-reset");

    // Drag near the right edge to exercise trackball roll with real pointer events.
    const edge = await viewer.locator("#scene").boundingBox();
    await page.mouse.move(edge.x + edge.width - 25, edge.y + edge.height * .45);
    await page.mouse.down();
    await page.mouse.move(edge.x + edge.width - 25, edge.y + edge.height * .65, {steps: 12});
    await page.mouse.up();
    await settled();
    check((await canvasPixels()).centralFingerprint !== initial.centralFingerprint,
      "Edge dragging rolls the rendered scene");
    await screenshot("viewer-rolled");
    await viewer.locator("#reset").click();
    await settled();
    check((await canvasPixels()).image === initial.image, "Reset restores the scene after edge roll");

    // Hide endpoint markers so they cannot completely occlude the tiny
    // synthetic focal boundary while testing that layer's visibility control.
    await viewer.locator("#toggleEndpoints").click();
    await settled();
    const layerBaseline = await canvasPixels();
    await viewer.locator("#toggleFocal").click();
    await settled();
    check(await viewer.locator("#toggleFocal").getAttribute("aria-pressed") === "false",
      "Focal-boundary toggle exposes its disabled state");
    check((await canvasPixels()).centralFingerprint !== layerBaseline.centralFingerprint,
      "Focal-boundary toggle changes central rendered pixels");
    await viewer.locator("#toggleFocal").click();
    await viewer.locator("#toggleEndpoints").click();
    await viewer.locator("#advancedControls > summary").click();
    await viewer.locator('.segment[data-mode="focal"]').click();
    await settled();
    check(await viewer.locator('.segment[data-mode="focal"]').getAttribute("aria-pressed") === "true",
      "Advanced path display selects the focal segment");
    check((await canvasPixels()).centralFingerprint !== initial.centralFingerprint,
      "Focal-segment selection changes the rendered paths");
    await viewer.locator('.segment[data-mode="points"]').click();
    await viewer.locator('.endpoint-mode[data-endpoint-mode="Dilated1"]').click();
    await settled();
    check(await viewer.locator('.endpoint-mode[data-endpoint-mode="Dilated1"]').getAttribute("aria-pressed") === "true",
      "Endpoint sensitivity control changes its selected state");
    check((await viewer.locator("#endpointStatus").innerText()).includes("1-voxel"),
      "Endpoint status reports the selected sensitivity mode");
    await screenshot("viewer-advanced");
    await viewer.locator('.endpoint-mode[data-endpoint-mode="Exact"]').click();
    await viewer.locator("#advancedControls > summary").click();
    await page.setViewportSize({width: 500, height: 900});
    await settled();
    await viewer.locator("#reset").click();
    await settled();
    check(await viewer.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth + 2),
      "Narrow viewer fits the viewport without horizontal overflow");
    const narrow = await canvasPixels();
    check(narrow.centralPainted > 20 && narrow.centralColors > 8, "Narrow viewer still paints nonuniform scene content");
    const narrowBounds = await projectedSceneBounds();
    result.narrow_scene_bounds = narrowBounds;
    const markerPadding = 9; // Largest default marker: radius 7 plus half its 3-pixel stroke.
    check(narrowBounds.point_count > 0 && narrowBounds.finite &&
      narrowBounds.min_x >= markerPadding && narrowBounds.min_y >= markerPadding &&
      narrowBounds.max_x <= narrowBounds.canvas_width - markerPadding &&
      narrowBounds.max_y <= narrowBounds.canvas_height - markerPadding,
    "After reset, the complete projected scene fits inside the narrow canvas with marker padding", narrowBounds);
    await screenshot("viewer-narrow", true);
    const cohortRoot = fs.realpathSync(manifest.cohort.relocated_root);
    const cohortURL = pathToFileURL(path.join(cohortRoot, manifest.cohort.report)).href;
    const cohortPaths = bundleFilePaths(cohortRoot);
    check(cohortPaths.filter(p => /[.]html?$/i.test(p)).length === 1 && path.dirname(manifest.cohort.report) === ".",
      "Cohort uses one HTML report at the study root");
    if (manifest.path_mode === "near-limit") {
      check(Math.max(...cohortPaths.map(p => p.length)) === 259, "Relocated cohort reaches the Windows path boundary");
    }
    await page.setViewportSize({width: 1440, height: 1000});
    await open(cohortURL);
    const runs = manifest.cohort.runs;
    check(await page.locator("#ftq-run option").count() === runs.length && runs.length === 4,
      "Cohort offers all four subject/focal combinations");
    for (const [i, run] of runs.entries()) {
      const label = `${run.subject} / ${run.focal_region}`;
      await page.locator("#ftq-run").selectOption(String(i));
      const title = await page.locator("#ftq-current-run").innerText();
      check(title.includes(run.subject) && title.includes(run.focal_region), `Selected run heading: ${label}`);
      const denominator = new RegExp(`^Total number of fNode voxels: ${run.focal_voxels}\\.(?: |$)`);
      check(await page.locator("#ftq-primary").getByText(denominator).isVisible(),
        `RSC denominator matches the selected input mask: ${label}`);
      check(await page.locator("#ftq-viewer").count() === 0, `Selection clears previous viewer: ${label}`);
      await inspectResultsFirst(cohortRoot, run, false);
      await inspectRun(cohortRoot, label, run);
      await openViewer();
      const scene = await canvasPixels();
      check(scene.centralPainted > 20 && scene.centralColors > 8, `Selected viewer draws a new scene: ${label}`);
      const description = await viewer.evaluate(() => data.metadata.title);
      check(description.includes(run.subject) && description.includes(run.focal_region), `Viewer uses selected run's scene: ${label}`, description);
      check(page.url().split("#")[0] === cohortURL, `Subject switching stays inside the cohort document: ${label}`);
      await screenshot(`cohort-viewer-${i + 1}`);
    }
    check(result.page_errors.length === 0, "No report or viewer JavaScript errors", result.page_errors);
    result.passed = true;
  } catch (error) {
    result.error = error.stack || String(error);
    if (page) {
      result.failure_url = page.url();
      try { await screenshot("failure"); } catch (_) { /* Preserve the primary failure. */ }
    }
    process.exitCode = 1;
  } finally {
    fs.writeFileSync(path.join(evidenceDir, "browser-check.json"), JSON.stringify(result, null, 2) + "\n");
    if (browser) await browser.close();
    process.stdout.write(JSON.stringify({passed: result.passed, checks: result.checks.length,
      pages: result.pages.length, local_references: result.local_references_checked,
      error: result.error || null}, null, 2) + "\n");
  }
})();
